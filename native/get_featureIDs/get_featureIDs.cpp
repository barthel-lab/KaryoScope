// KaryoScope: get_featureIDs
//
// Reads a FASTA/FASTQ input (plain or gzipped, or stdin) and queries every
// k-mer against a KMC database, emitting a single combined BED file with one
// record per run of consecutive k-mer positions sharing the same feature id
// (KMC's "counter" field repurposed to hold the feature id).
//
// Output filename pattern:
//   <output_dir>/<prefix>.combined.presmoothed.featureIDs.bed
//
// The default prefix is "<fasta-basename-without-ext>.<kmc-basename>", and
// the default output_dir is the input FASTA's parent directory (or "." for
// stdin). The "presmoothed" tag in the name flags that this BED is the raw
// per-position output; downstream Python code applies hierarchy-aware
// smoothing before splitting per feature set.
//
// This file is adapted from the version in the KaryoScope archive repo.
// Only the include paths changed: header references to the vendored KMC
// API and cxxopts now rely on the Makefile's -I flags rather than relative
// "../external/..." paths.

#include "kmc_api/kmc_file.h"
#include "kmc_api/kmer_api.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <condition_variable>
#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <regex>
#include <span>
#include <sstream>
#include <string_view>
#include <string>
#include <thread>
#include <vector>
#include <zlib.h>

namespace fs = std::filesystem;

/**
 * @brief Resolves the KMC database base path from a directory, full filename, or base name.
 * Validates that both .kmc_pre and .kmc_suf exist.
 */
std::string resolve_kmc_db_path(const std::string& input_path) {
    fs::path p(input_path);

    // Case 1: Input is a directory. Search for valid pairs.
    if (fs::exists(p) && fs::is_directory(p)) {
        std::vector<fs::path> valid_dbs;
        for (const auto& entry : fs::directory_iterator(p)) {
            // Look for .kmc_pre files
            if (entry.is_regular_file() && entry.path().extension() == ".kmc_pre") {
                // Check for corresponding .kmc_suf
                fs::path base_path = entry.path().parent_path() / entry.path().stem();
                fs::path suf_path = base_path;
                suf_path += ".kmc_suf";

                if (fs::exists(suf_path)) {
                    valid_dbs.push_back(base_path);
                }
            }
        }

        if (valid_dbs.empty()) {
            throw std::runtime_error("No valid KMC database (pair of .kmc_pre and .kmc_suf) found in directory: " + input_path);
        } else if (valid_dbs.size() > 1) {
            std::string msg = "Multiple KMC databases found in directory " + input_path + ". Please specify one explicitly:\n";
            for(const auto& db : valid_dbs) msg += "  " + db.string() + "\n";
            throw std::runtime_error(msg);
        }
        return valid_dbs[0].string();
    }

    // Case 2: Input ends in .kmc_pre or .kmc_suf. Strip and verify pair.
    std::string input_str = p.string();
    if (input_str.size() > 8 && input_str.substr(input_str.size() - 8) == ".kmc_pre") {
        std::string base = input_str.substr(0, input_str.size() - 8);
        if (fs::exists(base + ".kmc_suf")) return base;
        throw std::runtime_error("Found .kmc_pre but missing corresponding .kmc_suf for: " + base);
    }
    if (input_str.size() > 8 && input_str.substr(input_str.size() - 8) == ".kmc_suf") {
        std::string base = input_str.substr(0, input_str.size() - 8);
        if (fs::exists(base + ".kmc_pre")) return base;
        throw std::runtime_error("Found .kmc_suf but missing corresponding .kmc_pre for: " + base);
    }

    // Case 3: Treat as base name. Check if extensions exist.
    fs::path pre_path = p; pre_path += ".kmc_pre";
    fs::path suf_path = p; suf_path += ".kmc_suf";

    if (fs::exists(pre_path) && fs::exists(suf_path)) {
        return p.string();
    }

    throw std::runtime_error("KMC database not found. Checked for directory search, specific .kmc_pre/.kmc_suf extensions, and base name path at: " + input_path);
}

/**
 * @brief Returns the filename (without path) from a given path.
 */
std::string_view get_basename(std::string_view path) {
    auto pos = path.find_last_of("/\\");
    if (pos == std::string_view::npos) return path;
    return path.substr(pos + 1);
}

/**
 * @brief Removes known FASTA/FASTQ extensions from a filename.
 */
std::string_view remove_fasta_extensions(std::string_view filename) {
    static constexpr std::array<std::string_view, 10> extensions = {
        ".fasta.gz", ".fa.gz", ".fna.gz", ".fastq.gz", ".fq.gz",
        ".fasta", ".fa", ".fna", ".fastq", ".fq"
    };
    for (auto ext : extensions) {
        if (filename.size() >= ext.size() &&
            std::equal(ext.rbegin(), ext.rend(), filename.rbegin(),
                [](char a, char b) { return std::tolower(a) == std::tolower(b); })) {
            return filename.substr(0, filename.size() - ext.size());
        }
    }
    return filename;
}

/**
 * @brief Returns the prefix of a FASTA/FASTQ file (filename without extension).
 */
std::string get_fasta_prefix(std::string_view path) {
    auto filename = get_basename(path);
    auto prefix = remove_fasta_extensions(filename);
    return std::string(prefix);
}

/**
 * @brief Structure representing a sequence processing task.
 */
struct SequenceTask {
    size_t index;
    std::string name;
    std::string sequence;
};

/**
 * @brief Structure representing the result of processing a sequence.
 */
struct SequenceResult {
    size_t index;
    std::string seq_name;
    std::vector<std::string> original_bed_lines;

    SequenceResult(size_t idx = 0, std::string name = "")
        : index(idx), seq_name(std::move(name)) {}

    SequenceResult(SequenceResult&& other) noexcept
        : index(other.index),
          seq_name(std::move(other.seq_name)),
          original_bed_lines(std::move(other.original_bed_lines)) {}

    SequenceResult& operator=(SequenceResult&& other) noexcept {
        if (this != &other) {
            index = other.index;
            seq_name = std::move(other.seq_name);
            original_bed_lines = std::move(other.original_bed_lines);
        }
        return *this;
    }

    bool operator<(const SequenceResult& other) const {
        return index < other.index;
    }
};

/**
 * @brief Thread-safe queue for inter-thread communication.
 */
template<typename T>
class ThreadSafeQueue {
private:
    std::queue<T> queue;
    mutable std::mutex mutex;
    std::condition_variable cond_var;
    std::atomic_bool finished{false};

public:
    void push(T item) {
        {
            std::lock_guard<std::mutex> lock(mutex);
            queue.push(std::move(item));
        }
        cond_var.notify_one();
    }

    bool pop(T& item) {
        std::unique_lock<std::mutex> lock(mutex);
        cond_var.wait(lock, [this] { return !queue.empty() || finished.load(); });
        if (queue.empty() && finished.load()) return false;
        if (queue.empty()) return false;
        item = std::move(queue.front());
        queue.pop();
        return true;
    }

    void signal_finished() {
        finished.store(true);
        cond_var.notify_all();
    }

    bool is_finished() const { return finished.load(); }
    bool empty_unsafe() { return queue.empty(); }
    bool empty() const {
        std::lock_guard<std::mutex> lock(mutex);
        return queue.empty();
    }
};

/**
 * @brief Clusters and formats original per-position feature IDs as BED lines.
 */
std::vector<std::string> cluster_and_format_original_bed(std::string_view seq_name, std::span<const uint32_t> feature_ids) {
    std::vector<std::string> bed_lines;
    size_t n_features = feature_ids.size();
    if (n_features == 0) return bed_lines;
    size_t block_start_pos = 0;
    uint32_t previous_feature_id = feature_ids[0];
    for (size_t i = 1; i < n_features; ++i) {
        uint32_t current_feature_id = feature_ids[i];
        size_t current_pos = i;
        if (current_feature_id != previous_feature_id) {
            bed_lines.push_back(std::string(seq_name) + "\t" + std::to_string(block_start_pos) + "\t" +
                                std::to_string(current_pos) + "\t" + std::to_string(previous_feature_id));
            block_start_pos = current_pos;
            previous_feature_id = current_feature_id;
        }
    }
    bed_lines.push_back(std::string(seq_name) + "\t" + std::to_string(block_start_pos) + "\t" +
                        std::to_string(n_features) + "\t" + std::to_string(previous_feature_id));
    return bed_lines;
}

/**
 * @brief Processes a single sequence and pushes the result to the results queue.
 */
void process_sequence_bed_worker(
    size_t seq_index,
    std::string_view seq_name,
    std::string_view sequence,
    CKMCFile& kmc_db,
    ThreadSafeQueue<SequenceResult>& results_queue)
{
    SequenceResult result(seq_index, std::string(seq_name));
    if (sequence.empty()) {
        results_queue.push(std::move(result));
        return;
    }

    std::vector<uint32_t> feature_ids_output;
    kmc_db.GetCountersForRead(std::string(sequence), feature_ids_output);

    // Perform original per-position clustering
    result.original_bed_lines = cluster_and_format_original_bed(seq_name, feature_ids_output);

    results_queue.push(std::move(result));
}

/**
 * @brief Worker thread main loop for processing sequence tasks.
 */
void worker_thread_function(
    ThreadSafeQueue<SequenceTask>& task_queue,
    ThreadSafeQueue<SequenceResult>& results_queue,
    CKMCFile& kmc_db)
{
    SequenceTask task;
    while (task_queue.pop(task)) {
        process_sequence_bed_worker(task.index, task.name, task.sequence, kmc_db, results_queue);
    }
}

/**
 * @brief Reads a line from a gzipped file.
 */
bool getline_gz(gzFile gz, std::string& line) {
    static constexpr int buffer_size = 4096;
    char buffer[buffer_size];
    line.clear();
    while (true) {
        char* res = gzgets(gz, buffer, buffer_size);
        if (!res) return !line.empty();
        char* newline = strchr(buffer, '\n');
        if (newline) {
            *newline = '\0';
            line += buffer;
            return true;
        } else {
            line += buffer;
        }
    }
}

/**
 * @brief RAII wrapper for gzFile.
 */
class GzFileRAII {
public:
    explicit GzFileRAII(const char* path, const char* mode) : file_(gzopen(path, mode)) {}
    ~GzFileRAII() { if (file_) gzclose(file_); }
    GzFileRAII(const GzFileRAII&) = delete;
    GzFileRAII& operator=(const GzFileRAII&) = delete;
    GzFileRAII(GzFileRAII&& other) noexcept : file_(other.file_) { other.file_ = nullptr; }
    GzFileRAII& operator=(GzFileRAII&& other) noexcept {
        if (this != &other) {
            if (file_) gzclose(file_);
            file_ = other.file_;
            other.file_ = nullptr;
        }
        return *this;
    }
    gzFile get() const { return file_; }
    explicit operator bool() const { return file_ != nullptr; }
private:
    gzFile file_;
};


int main(int argc, char** argv) {
    // --- Argument Parsing ---
    cxxopts::Options options("get_featureIDs", "K-mer feature analysis tool");
    options.add_options()
        ("d,db", "KMC database path (Name, Name.kmc_pre, Name.kmc_suf, or Directory)", cxxopts::value<std::string>())
        ("i,input", "Input FASTA/FASTQ file, or '-' for stdin", cxxopts::value<std::string>())
        ("o,output", "Output directory (optional, defaults to input file directory)", cxxopts::value<std::string>())
        ("t,threads", "Number of threads", cxxopts::value<int>()->default_value("0"))
        ("p,prefix", "Output filename prefix (optional, defaults to Input.KMC_DB)", cxxopts::value<std::string>())
        ("input-format", "Format of stdin ('fasta' or 'fastq'). Required when using '-i -'", cxxopts::value<std::string>()->default_value(""))
        ("h,help", "Print usage");
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (!result.count("db") || !result.count("input")) {
        std::cerr << "Error: --db and --input are required." << std::endl;
        std::cout << options.help() << std::endl;
        return 1;
    }

    std::string raw_kmc_path = result["db"].as<std::string>();
    std::string fasta_filename_path = result["input"].as<std::string>();

    // --- Resolve KMC Database Path ---
    std::string kmc_db_path;
    try {
        kmc_db_path = resolve_kmc_db_path(raw_kmc_path);
    } catch (const std::exception& e) {
        std::cerr << "Error resolving KMC database: " << e.what() << std::endl;
        return 1;
    }

    int requested_threads = result["threads"].as<int>();
    std::string input_format = result["input-format"].as<std::string>();
    const bool is_stdin = (fasta_filename_path == "-");

    if (is_stdin && input_format.empty()) {
        std::cerr << "Error: --input-format must be set to 'fasta' or 'fastq' when reading from stdin ('-i -')." << std::endl;
        return 1;
    }
    if (!input_format.empty() && input_format != "fasta" && input_format != "fastq") {
        std::cerr << "Error: --input-format must be 'fasta' or 'fastq'." << std::endl;
        return 1;
    }

    // --- Determine Thread Count ---
    int num_threads;
    unsigned int hardware_threads = std::thread::hardware_concurrency();
    if (requested_threads <= 0) {
        num_threads = hardware_threads > 0 ? hardware_threads : 4;
    } else {
        num_threads = (hardware_threads > 0 && static_cast<unsigned int>(requested_threads) > hardware_threads)
            ? hardware_threads : requested_threads;
    }

    // --- Determine Output Directory ---
    std::string output_dir;
    if (result.count("output")) {
        output_dir = result["output"].as<std::string>();
    } else {
        if (is_stdin) {
            output_dir = ".";
        } else {
            fs::path p(fasta_filename_path);
            if (p.has_parent_path()) {
                output_dir = p.parent_path().string();
            } else {
                output_dir = ".";
            }
        }
        std::cerr << "INFO: Output directory defaulting to: " << output_dir << std::endl;
    }

    std::cerr << "INFO: Launching " << num_threads << " worker threads." << std::endl;

    // --- Create Output Directory ---
    try {
        fs::create_directories(output_dir);
    } catch (const std::exception& e) {
        std::cerr << "Error creating output directory '" << output_dir << "': " << e.what() << std::endl;
        return 1;
    }

    // --- Generate Base Filename ---
    // Use resolved KMC path for basename generation to be accurate
    std::string kmc_basename(get_basename(kmc_db_path));
    std::string base_output_name;

    if (result.count("prefix")) {
        base_output_name = result["prefix"].as<std::string>();
        std::cerr << "INFO: Using user-provided output prefix: " << base_output_name << std::endl;
    } else {
        std::string fasta_prefix;
        if (is_stdin) {
            fasta_prefix = "stdin";
        } else {
            fasta_prefix = get_fasta_prefix(fasta_filename_path);
        }
        if (kmc_basename.empty() || fasta_prefix.empty()) {
             std::cerr << "Error: Could not determine base names for output files." << std::endl;
             return 1;
        }
        base_output_name = fasta_prefix + "." + kmc_basename;
        std::cerr << "INFO: Using generated output prefix: " << base_output_name << std::endl;
    }

    // --- Construct Output Paths and Open Files ---
    fs::path dir_path = output_dir;

    // Holder for combined file paths and streams
    std::string combined_original_full_path;
    std::ofstream combined_original_output_stream;

    // Set paths for combined files (Implied default: combined_only=true)
    combined_original_full_path = (dir_path / (base_output_name + ".combined.presmoothed.featureIDs.bed")).string();

    // Open combined stream
    combined_original_output_stream.open(combined_original_full_path);
    if (!combined_original_output_stream.is_open()) {
        std::cerr << "Error opening file: " << combined_original_full_path << std::endl;
        return 1;
    }

    // --- Open KMC DB ---
    CKMCFile kmc_db;
    if (!kmc_db.OpenForRA(kmc_db_path)) {
        std::cerr << "Error opening KMC database: " << kmc_db_path << std::endl;
        return 1;
    }

    // --- Threading Setup ---
    ThreadSafeQueue<SequenceTask> task_queue;
    ThreadSafeQueue<SequenceResult> results_queue;
    std::vector<std::thread> worker_threads;

    // --- Output Thread ---
    std::thread output_thread(
        [&](ThreadSafeQueue<SequenceResult>& results_queue) {
            size_t next_expected_index = 0;
            std::map<size_t, SequenceResult> buffered_results;

            auto process_and_write_result = [&](SequenceResult& result) {
                // Write combined/original BED lines directly
                for (const auto& line : result.original_bed_lines) {
                    combined_original_output_stream << line << "\n";
                }
            };

            SequenceResult result;
            while (true) {
                if (!results_queue.pop(result)) {
                    if (results_queue.is_finished() && results_queue.empty_unsafe()) break;
                    continue;
                }
                buffered_results[result.index] = std::move(result);
                while (true) {
                    auto it = buffered_results.find(next_expected_index);
                    if (it == buffered_results.end()) break;
                    process_and_write_result(it->second);
                    buffered_results.erase(it);
                    ++next_expected_index;
                }
            }
            while (!buffered_results.empty()) {
                 auto it = buffered_results.find(next_expected_index);
                 if (it == buffered_results.end()) {
                     std::cerr << "FATAL: Missing result for sequence index " << next_expected_index << std::endl;
                     break;
                 }
                 process_and_write_result(it->second);
                 buffered_results.erase(it);
                 ++next_expected_index;
            }
        },
        std::ref(results_queue)
    );

    // --- Launch Worker Threads ---
    for (int i = 0; i < num_threads; ++i) {
        worker_threads.emplace_back(
            worker_thread_function,
            std::ref(task_queue),
            std::ref(results_queue),
            std::ref(kmc_db)
        );
    }

    // --- Main Thread: Read FASTA/FASTQ and Dispatch Tasks ---
    bool is_gzipped = false;
    bool is_fastq = false;

    if (is_stdin) {
        std::cerr << "INFO: Reading from standard input..." << std::endl;
        is_fastq = (input_format == "fastq");
    } else {
        auto lower_ext = [](std::string_view s) {
            std::string ext(s);
            std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
            return ext;
        };
        std::string lower_filename = lower_ext(fasta_filename_path);
        if (lower_filename.size() >= 3 && lower_filename.substr(lower_filename.size() - 3) == ".gz") {
            is_gzipped = true;
            lower_filename = lower_filename.substr(0, lower_filename.size() - 3);
        }
        if (lower_filename.size() >= 6 &&
            (lower_filename.substr(lower_filename.size() - 6) == ".fastq" ||
             lower_filename.substr(lower_filename.size() - 3) == ".fq")) {
            is_fastq = true;
        }
    }

    std::ifstream plain_file_stream;
    GzFileRAII gz_file_handle(nullptr, "rb");
    std::istream* input_stream = nullptr;

    if (is_stdin) {
        input_stream = &std::cin;
    } else if (is_gzipped) {
        gz_file_handle = GzFileRAII(fasta_filename_path.c_str(), "rb");
        if (!gz_file_handle) {
            std::cerr << "Error opening gzipped file: " << fasta_filename_path << std::endl;
            return 1;
        }
    } else {
        plain_file_stream.open(fasta_filename_path);
        if (!plain_file_stream.is_open()) {
            std::cerr << "Error opening plain file: " << fasta_filename_path << std::endl;
            return 1;
        }
        input_stream = &plain_file_stream;
    }

    std::function<bool(std::string&)> get_next_line;
    if (is_gzipped && !is_stdin) {
        get_next_line = [&](std::string& l) { return getline_gz(gz_file_handle.get(), l); };
    } else {
        get_next_line = [&](std::string& l) { return !!std::getline(*input_stream, l); };
    }

    std::string line;
    std::string current_sequence;
    std::string current_seq_name;
    bool first_sequence = true;
    size_t current_index = 0;

    if (!is_fastq) { // FASTA processing
        while (get_next_line(line)) {
            line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
            line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!first_sequence) {
                    task_queue.push({current_index++, current_seq_name, current_sequence});
                }
                first_sequence = false;
                current_sequence.clear();
                current_seq_name = line.substr(1);
                size_t first_space = current_seq_name.find_first_of(" \t");
                if (first_space != std::string::npos)
                    current_seq_name = current_seq_name.substr(0, first_space);
            } else {
                current_sequence += line;
            }
        }
        if (!current_seq_name.empty() && !current_sequence.empty()) {
            task_queue.push({current_index++, current_seq_name, current_sequence});
        }
    } else { // FASTQ processing
        while (true) {
            std::array<std::string, 4> record_lines;
            bool record_complete = true;
            for(int i = 0; i < 4; ++i) {
                if (!get_next_line(record_lines[i])) {
                    record_complete = false;
                    break;
                }
            }
            if (!record_complete) break;

            std::string& name_line = record_lines[0];
            std::string& seq_line = record_lines[1];

            name_line.erase(0, name_line.find_first_not_of(" \t\n\r\f\v"));
            name_line.erase(name_line.find_last_not_of(" \t\n\r\f\v") + 1);
            seq_line.erase(0, seq_line.find_first_not_of(" \t\n\r\f\v"));
            seq_line.erase(seq_line.find_last_not_of(" \t\n\r\f\v") + 1);

            if (!name_line.empty() && !seq_line.empty() && name_line[0] == '@') {
                std::string seq_name = name_line.substr(1);
                size_t first_space = seq_name.find_first_of(" \t");
                if (first_space != std::string::npos)
                    seq_name = seq_name.substr(0, first_space);
                task_queue.push({current_index++, seq_name, seq_line});
            } else {
                std::cerr << "Warning: Skipping malformed FASTQ record near index " << current_index << std::endl;
            }
        }
    }

    // --- Signal & Join ---
    task_queue.signal_finished();
    std::cerr << "INFO: Finished reading input. Waiting for workers (" << current_index << " sequences dispatched)..." << std::endl;
    for (auto& t : worker_threads) if (t.joinable()) t.join();

    std::cerr << "INFO: All workers finished." << std::endl;
    results_queue.signal_finished();

    std::cerr << "INFO: Waiting for output thread to finish writing..." << std::endl;
    if (output_thread.joinable()) output_thread.join();

    // --- Cleanup ---
    if (combined_original_output_stream.is_open()) combined_original_output_stream.close();

    std::cerr << "INFO: Processing complete. Output written to directory: " << output_dir << std::endl;
    return 0;
}
