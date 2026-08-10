"""Unit + integration tests for :mod:`karyoscope.core.karyotype`."""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import ClassVar

import pytest
from click.testing import CliRunner

import karyoscope.core.karyotype as karyotype_mod
from karyoscope.cli import main
from karyoscope.core.io.scaffold_map import MapRow
from karyoscope.core.karyotype import (
    PREDEFINED_SEX_SYSTEMS,
    RenderInput,
    _effective_hap,
    _haps_natural_sort_key,
    _legend_sort_key,
    convert_svg,
    get_expected_haps,
    render_karyotype,
)
from karyoscope.exceptions import KaryotypeError

# --- get_expected_haps ----------------------------------------------


class TestGetExpectedHaps:
    haps: ClassVar[list[str]] = ["hap1", "hap2"]

    def test_autosome_returns_all_haps(self) -> None:
        assert get_expected_haps("chr1", "male", self.haps, "XY") == self.haps
        assert get_expected_haps("chr1", "female", self.haps, "XY") == self.haps

    def test_xy_male_chrY_first_hap_only(self) -> None:
        assert get_expected_haps("chrY", "male", self.haps, "XY") == ["hap1"]

    def test_xy_male_chrX_second_hap_only(self) -> None:
        # chrX is at index 1 in ['chrY', 'chrX']; with two haps that maps
        # to ['hap2'].
        assert get_expected_haps("chrX", "male", self.haps, "XY") == ["hap2"]

    def test_xy_female_chrY_empty(self) -> None:
        assert get_expected_haps("chrY", "female", self.haps, "XY") == []

    def test_xy_female_chrX_both_haps(self) -> None:
        assert get_expected_haps("chrX", "female", self.haps, "XY") == self.haps

    def test_unknown_sex_no_expectation_for_sex_chroms(self) -> None:
        assert get_expected_haps("chrY", None, self.haps, "XY") == []
        assert get_expected_haps("chr1", None, self.haps, "XY") == self.haps

    def test_reference_uses_first_hap(self) -> None:
        assert get_expected_haps("chr1", "reference", self.haps, "XY") == ["hap1"]

    def test_zw_female_chrW_first_hap(self) -> None:
        # ZW: heterogametic_sex='female', sex_chromosomes=['chrZ', 'chrW'].
        # Index of chrW is 1, so first hap.
        assert get_expected_haps("chrW", "female", self.haps, "ZW") == ["hap2"]
        assert get_expected_haps("chrZ", "female", self.haps, "ZW") == ["hap1"]
        assert get_expected_haps("chrZ", "male", self.haps, "ZW") == self.haps
        assert get_expected_haps("chrW", "male", self.haps, "ZW") == []

    def test_unknown_system_raises(self) -> None:
        with pytest.raises(KaryotypeError, match="unknown sex-determination"):
            get_expected_haps("chr1", "male", self.haps, "ABC")

    def test_dict_system_accepted(self) -> None:
        # Custom dict matching the predefined schema.
        custom = PREDEFINED_SEX_SYSTEMS["XY"]
        assert get_expected_haps("chrY", "male", self.haps, custom) == ["hap1"]


class TestLegendSortKey:
    """The legend in the karyotype SVG sorts feature names with a
    specific layout: chromosomes (chr*) at the very top in natural
    order, then "categorized" (hierarchy root), then hierarchy-order
    features, then any unranked features alphabetical, then "novel"
    at the very bottom.
    """

    @staticmethod
    def _sorted(names: list[str], feature_order: list[str] | None = None) -> list[str]:
        return sorted(names, key=lambda f: _legend_sort_key(f, feature_order))

    def test_chromosome_fs_production_layout(self) -> None:
        """Production CHM13 chromosome feature set: chromosomes first
        in natural order, then categorized, then the categorical
        groupings in hierarchy.tsv order, then novel."""
        # The order the hierarchy file gives for the chromosome FS's
        # non-leaf categories.
        hierarchy_order = [
            "autosome",
            "acrocentric",
            "metacentric",
            "submetacentric",
            "sex",
        ]
        # Simulated set of features seen in data for chromosome FS.
        seen = [
            "novel",
            "chr1",
            "chr22",
            "chr2",
            "chrX",
            "chrY",
            "autosome",
            "acrocentric",
            "metacentric",
            "submetacentric",
            "sex",
            "categorized",
        ]
        result = self._sorted(seen, feature_order=hierarchy_order)
        assert result == [
            # Chromosomes natural-ordered at the very top
            "chr1",
            "chr2",
            "chr22",
            "chrX",
            "chrY",
            # Then "categorized" (hierarchy root, pinned)
            "categorized",
            # Then categorical groupings in hierarchy.tsv order
            "autosome",
            "acrocentric",
            "metacentric",
            "submetacentric",
            "sex",
            # Novel at the very bottom
            "novel",
        ]

    def test_unknown_features_alphabetical_after_hierarchy_entries(self) -> None:
        """Future-database fallback: features in data but not in
        hierarchy_order end up between known hierarchy entries and
        novel, sorted alphabetically."""
        hierarchy_order = ["alpha", "beta"]
        seen = ["zeta", "alpha", "novel", "beta", "gamma", "delta"]
        result = self._sorted(seen, feature_order=hierarchy_order)
        # alpha, beta first (hierarchy order); delta/gamma/zeta
        # alphabetical at the end (unknown bucket); novel last.
        assert result == ["alpha", "beta", "delta", "gamma", "zeta", "novel"]

    def test_no_feature_order_alphabetical(self) -> None:
        """When feature_order is None (e.g. caller doesn't have a
        hierarchy), non-special features fall into one alphabetical
        bucket between chromosomes-and-categorized and novel."""
        seen = ["zeta", "alpha", "novel", "categorized", "chrX", "chr1", "beta"]
        result = self._sorted(seen)
        assert result == [
            "chr1",
            "chrX",
            "categorized",
            "alpha",
            "beta",
            "zeta",
            "novel",
        ]

    def test_chromosomes_natural_chrM_after_numeric(self) -> None:
        # Non-numeric chr* suffixes (chrM, chrW, chrZ etc.) sort
        # alphabetically AFTER chr1..chr22.
        assert self._sorted(["chr1", "chrM", "chr2", "chrY", "chrX"]) == [
            "chr1",
            "chr2",
            "chrM",
            "chrX",
            "chrY",
        ]

    def test_categorized_above_other_features(self) -> None:
        """``"categorized"`` (hierarchy root) goes immediately after
        chromosomes, before any other category."""
        seen = ["acrocentric", "categorized", "autosome", "chr1"]
        result = self._sorted(seen, feature_order=["autosome", "acrocentric"])
        assert result == ["chr1", "categorized", "autosome", "acrocentric"]

    def test_novel_always_last(self) -> None:
        # Even with a feature_order that doesn't mention novel, it
        # pins to the bottom.
        assert self._sorted(["novel", "alpha", "chr1"], feature_order=["alpha"]) == [
            "chr1",
            "alpha",
            "novel",
        ]


class TestHapsNaturalSortKey:
    """Hap column ordering follows the HPRC convention: paternal first
    (= hap1), maternal second (= hap2). Numeric hapN forms come before
    biological labels; ``unassigned`` is always last.
    """

    @staticmethod
    def _sorted(haps: list[str]) -> list[str]:
        return sorted(haps, key=_haps_natural_sort_key)

    def test_paternal_before_maternal(self) -> None:
        # HPRC convention: paternal = hap1 (first column).
        assert self._sorted(["maternal", "paternal"]) == ["paternal", "maternal"]

    def test_hap1_before_hap2(self) -> None:
        assert self._sorted(["hap2", "hap1"]) == ["hap1", "hap2"]

    def test_hap_numbers_ordered_numerically(self) -> None:
        # Lexical sort would put "hap10" before "hap2"; we want numeric.
        assert self._sorted(["hap2", "hap10", "hap1"]) == ["hap1", "hap2", "hap10"]

    def test_unassigned_always_last(self) -> None:
        assert self._sorted(["unassigned", "paternal", "maternal"]) == [
            "paternal",
            "maternal",
            "unassigned",
        ]

    def test_hapN_comes_before_biological_labels(self) -> None:
        assert self._sorted(["paternal", "hap1"]) == ["hap1", "paternal"]

    def test_other_labels_alphabetical(self) -> None:
        assert self._sorted(["zeta", "alpha", "beta"]) == ["alpha", "beta", "zeta"]


class TestGetExpectedHapsWithDataDrivenInference:
    """When per-chrom-hap data is passed, the heterogametic chromosome's
    expected hap is inferred from where chrY actually lives, not from
    haplotypes[0]. Fixes the maternal/paternal labelling bug where
    sort order puts ``maternal`` first but chrY is actually paternal.
    """

    # HG002-style: maternal/paternal labels, chrY lives in paternal,
    # chrX lives in maternal (XY male). Hap order follows the HPRC
    # convention (paternal = hap1 first, maternal = hap2 second);
    # matches what ``_haps_natural_sort_key`` produces for this set.
    hg002_haps: ClassVar[list[str]] = ["paternal", "maternal"]
    hg002_seqs: ClassVar[dict[str, dict[str, list[str]]]] = {
        "chr1": {"paternal": ["chr1_PATERNAL"], "maternal": ["chr1_MATERNAL"]},
        "chrX": {"paternal": [], "maternal": ["chrX_MATERNAL"]},
        "chrY": {"paternal": ["chrY_PATERNAL"], "maternal": []},
    }

    def test_male_hg002_chrY_inferred_to_paternal(self) -> None:
        result = get_expected_haps(
            "chrY",
            "male",
            self.hg002_haps,
            "XY",
            sequences_per_chrom_hap=self.hg002_seqs,
        )
        assert result == ["paternal"], (
            "data-driven inference should put chrY on paternal, "
            f"not the alphabetical-first hap; got {result}"
        )

    def test_male_hg002_chrX_inferred_to_maternal(self) -> None:
        result = get_expected_haps(
            "chrX",
            "male",
            self.hg002_haps,
            "XY",
            sequences_per_chrom_hap=self.hg002_seqs,
        )
        assert result == ["maternal"], (
            f"data-driven inference should put chrX on maternal; got {result}"
        )

    def test_autosomes_unaffected_by_inference(self) -> None:
        # chr1 should still return both haps regardless of inference.
        result = get_expected_haps(
            "chr1",
            "male",
            self.hg002_haps,
            "XY",
            sequences_per_chrom_hap=self.hg002_seqs,
        )
        assert result == self.hg002_haps

    def test_unknown_sex_unaffected_by_inference(self) -> None:
        # With sex=None, sex chromosomes still get [] (data-driven only).
        # The inference parameter is irrelevant.
        for chrom in ("chrX", "chrY"):
            result = get_expected_haps(
                chrom,
                None,
                self.hg002_haps,
                "XY",
                sequences_per_chrom_hap=self.hg002_seqs,
            )
            assert result == [], (
                f"sex=None should always return [] for sex chromosomes; got {result} for {chrom}"
            )

    def test_falls_back_to_sort_order_when_no_chrY_data(self) -> None:
        """Cancer case: --sex male but chrY is lost. Inference can't
        determine the heterogametic hap, so fall back to haplotypes[0].
        With the HPRC sort order (paternal first), the fallback lands
        on the biologically-correct hap -- the empty chrY column is
        labelled paternal as a real chrY would be.
        """
        cancer_seqs = {
            "chrX": {"paternal": ["chrX_PATERNAL"], "maternal": ["chrX_MATERNAL"]},
            "chrY": {"paternal": [], "maternal": []},  # chrY lost
        }
        result_y = get_expected_haps(
            "chrY",
            "male",
            self.hg002_haps,
            "XY",
            sequences_per_chrom_hap=cancer_seqs,
        )
        # Fall-back to haplotypes[0] = "paternal" (HPRC convention).
        assert result_y == ["paternal"]

    def test_falls_back_to_sort_order_when_chrY_in_both_haps(self) -> None:
        """Ambiguous case: chrY data appears in both haps (unusual but
        possible with mis-labelled contigs). Inference returns None;
        fall back to haplotypes[0]."""
        weird_seqs = {
            "chrY": {"paternal": ["chrY_PATERNAL"], "maternal": ["chrY_MATERNAL"]},
        }
        result = get_expected_haps(
            "chrY",
            "male",
            self.hg002_haps,
            "XY",
            sequences_per_chrom_hap=weird_seqs,
        )
        # Ambiguous → fall back to haplotypes[0] = "paternal".
        assert result == ["paternal"]

    def test_hap1_hap2_convention_still_works(self) -> None:
        """For the conventional ``hap1``/``hap2`` labelling, chrY
        should still end up on hap1 -- via inference when data shows
        it there, via the fallback otherwise."""
        haps = ["hap1", "hap2"]
        seqs = {
            "chrX": {"hap1": [], "hap2": ["chrX_seq"]},
            "chrY": {"hap1": ["chrY_seq"], "hap2": []},
        }
        assert get_expected_haps("chrY", "male", haps, "XY", sequences_per_chrom_hap=seqs) == [
            "hap1"
        ]
        assert get_expected_haps("chrX", "male", haps, "XY", sequences_per_chrom_hap=seqs) == [
            "hap2"
        ]

    def test_no_per_chrom_hap_passed_keeps_legacy_behaviour(self) -> None:
        """Calls without the new parameter should behave identically
        to the original archive logic (backward compat)."""
        assert get_expected_haps("chrY", "male", ["hap1", "hap2"], "XY") == ["hap1"]
        assert get_expected_haps("chrX", "male", ["hap1", "hap2"], "XY") == ["hap2"]


# --- render_karyotype: unit-level (no external tools) ---------------


def _row(
    name: str,
    *,
    chrom: str,
    hap: str,
    stats: str = "PCQ",
    flipped: bool = False,
    length: int = 1000,
) -> MapRow:
    return MapRow(
        new_name=name,
        original_name=name.split("_")[-1],
        input_file="x.fa",
        hap=hap,
        chromosome=chrom,
        flipped=flipped,
        length=length,
        stats=stats,
    )


def _read_svg_text(p: Path) -> str:
    return p.read_text()


def _label_x(svg: str, label: str) -> float:
    """Return the x of the first ``<text ...>label</text>`` column header."""
    m = re.search(rf"<text\b([^>]*)>{re.escape(label)}</text>", svg)
    assert m is not None, f"column label {label!r} not found in SVG"
    xm = re.search(r'x="([-\d.]+)"', m.group(1))
    assert xm is not None, f"no x on label {label!r}"
    return float(xm.group(1))


class TestEffectiveHap:
    def test_infers_hap_from_contig_name_over_file_label(self) -> None:
        # Combined-FASTA case: scaffold labelled every contig with the
        # sample stem, but the contig name carries the true haplotype.
        row = _row(
            "chr16_GM00392_haplotype2-0000061",
            chrom="chr16",
            hap="GM00392",
        )
        assert _effective_hap(row) == "hap2"

    def test_detects_unassigned_from_contig_name(self) -> None:
        row = _row("chr2_GM00392_unassigned-0000409", chrom="chr2", hap="GM00392")
        assert _effective_hap(row) == "unassigned"

    def test_falls_back_to_file_label_when_no_marker(self) -> None:
        # One-file-per-haplotype convention: contig name has no marker, so
        # the file-level hap label remains authoritative.
        row = _row("chr1_hap1_ptg000001l", chrom="chr1", hap="hap1")
        assert _effective_hap(row) == "hap1"


class TestHaplotypeColumns:
    def test_columns_ordered_by_true_hap_not_size(self, tmp_path: Path) -> None:
        # GM00392 chr16: the dup is on hap2 (full chr16, large) while hap1
        # holds only 16p-ter (small). Both contigs carry the file-level
        # label "GM00392"; ordering by that label (then size) would draw
        # the large hap2 contig first/left. Deriving the hap from the
        # contig name must put hap1 left of hap2 regardless of size.
        ri = RenderInput(
            map_rows=[
                _row(
                    "chr16_GM00392_haplotype2-0000061",
                    chrom="chr16",
                    hap="GM00392",
                    stats="TPCQT",
                    length=174_000_000,
                ),
                _row(
                    "chr16_GM00392_haplotype1-0000013",
                    chrom="chr16",
                    hap="GM00392",
                    stats="TP",
                    length=15_000_000,
                ),
            ],
            binned_bed={
                "chr16_GM00392_haplotype2-0000061": [(0, 500, "rA"), (500, 1000, "rB")],
                "chr16_GM00392_haplotype1-0000013": [(0, 500, "rA")],
            },
        )
        out = tmp_path / "chr16.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#2ca02c", "rB": "#d62728"},
            mode="genome",
            seed_human_chromosomes=False,
            output_path=out,
        )
        text = _read_svg_text(out)
        # Both haplotype columns are labelled ("h1"/"h2")...
        assert "<text" in text and ">h1</text>" in text and ">h2</text>" in text
        # ...and hap1 is drawn to the left of hap2 (true-haplotype order,
        # not size order).
        assert _label_x(text, "h1") < _label_x(text, "h2")

    def test_unassigned_segregated_and_labelled(self, tmp_path: Path) -> None:
        # GM00392 chr2: both haplotypes fully assembled plus a tiny
        # unassigned fragment. The unassigned contig must land in its own
        # labelled column, to the right of the real haplotypes.
        ri = RenderInput(
            map_rows=[
                _row("chr2_GM00392_haplotype1-0000012", chrom="chr2", hap="GM00392", length=243),
                _row("chr2_GM00392_haplotype2-0000056", chrom="chr2", hap="GM00392", length=242),
                _row("chr2_GM00392_unassigned-0000409", chrom="chr2", hap="GM00392", length=59),
            ],
            binned_bed={
                "chr2_GM00392_haplotype1-0000012": [(0, 500, "rA")],
                "chr2_GM00392_haplotype2-0000056": [(0, 500, "rA")],
                "chr2_GM00392_unassigned-0000409": [(0, 500, "rA")],
            },
        )
        out = tmp_path / "chr2.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#2ca02c"},
            mode="genome",
            seed_human_chromosomes=False,
            output_path=out,
        )
        text = _read_svg_text(out)
        # Unassigned column is flagged with a compact "u" tag...
        assert ">u</text>" in text
        # ...and sits to the right of both real haplotypes.
        assert _label_x(text, "u") > _label_x(text, "h2")
        assert _label_x(text, "h2") > _label_x(text, "h1")


class TestRenderKaryotypeUnit:
    def test_renders_genome_mode_minimal(self, tmp_path: Path) -> None:
        # One contig, chr1, hap1, genome mode.
        ri = RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", stats="TPCQT", length=1000)],
            binned_bed={
                "chr1_hap1_a": [(0, 500, "rA"), (500, 1000, "rB")],
            },
        )
        out = tmp_path / "out.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#2ca02c", "rB": "#d62728"},
            mode="genome",
            output_path=out,
        )
        assert out.is_file()
        text = _read_svg_text(out)
        assert text.startswith("<?xml") or text.startswith("<svg")
        assert "<svg" in text
        # Should contain at least one rectangle (for the rendered bins
        # plus the background plus the sequence outline).
        assert text.count("<rect") >= 2
        # Colours appear in the output.
        assert "#2ca02c" in text or "#d62728" in text

    def test_renders_unknown_mode_raises(self, tmp_path: Path) -> None:
        ri = RenderInput(map_rows=[], binned_bed={})
        with pytest.raises(KaryotypeError, match="unknown mode"):
            render_karyotype(
                [ri],
                colors={},
                mode="wat",  # type: ignore[arg-type]
                output_path=tmp_path / "x.svg",
            )

    def test_centromere_mode_requires_centromere_ranges(self, tmp_path: Path) -> None:
        ri = RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_hap1_a": [(0, 500, "rA")]},
            centromere_ranges=None,
        )
        with pytest.raises(KaryotypeError, match="centromere mode requires"):
            render_karyotype(
                [ri],
                colors={"rA": "#000"},
                mode="centromere",
                output_path=tmp_path / "x.svg",
            )

    def test_centromere_mode_with_ranges_renders(self, tmp_path: Path) -> None:
        ri = RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", length=2000)],
            binned_bed={"chr1_hap1_a": [(800, 1200, "rA")]},
            centromere_ranges={"chr1_hap1_a": (800, 1200)},
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#abcdef"},
            mode="centromere",
            output_path=out,
        )
        assert out.is_file()
        assert "#abcdef" in out.read_text()

    def test_subtelomere_mode_skips_non_telomere_contigs(self, tmp_path: Path) -> None:
        ri = RenderInput(
            map_rows=[
                # Has telomere
                _row("chr1_hap1_a", chrom="chr1", hap="hap1", stats="TPCQ", length=1_000_000),
                # No telomere
                _row("chr1_hap1_b", chrom="chr1", hap="hap1", stats="PCQ", length=1_000_000),
            ],
            binned_bed={
                "chr1_hap1_a": [(0, 100_000, "rA")],
                "chr1_hap1_b": [(0, 100_000, "rB")],
            },
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#abcdef", "rB": "#fedcba"},
            mode="subtelomere",
            output_path=out,
        )
        text = out.read_text()
        # The contig with a telomere gets rendered; the one without doesn't.
        assert "#abcdef" in text
        assert "#fedcba" not in text

    def test_unknown_feature_is_hard_error(self, tmp_path: Path) -> None:
        # Stage 6d: missing colours used to silently fall back to
        # white (conflating real categorisation with the novel
        # sentinel). They now raise. The caller is expected to have
        # validated colours upstream via validate_colors.
        ri = RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_hap1_a": [(0, 1000, "mystery_feature")]},
        )
        out = tmp_path / "x.svg"
        with pytest.raises(KaryotypeError, match="mystery_feature"):
            render_karyotype([ri], colors={}, mode="genome", output_path=out)


# --- title + legend ------------------------------------------------


class TestTitleBand:
    def _basic_input(self) -> RenderInput:
        return RenderInput(
            map_rows=[_row("chr1_hap1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_hap1_a": [(0, 1000, "rA")]},
        )

    def test_title_text_appears_in_svg(self, tmp_path: Path) -> None:
        out = tmp_path / "x.svg"
        render_karyotype(
            [self._basic_input()],
            colors={"rA": "#2ca02c"},
            mode="genome",
            output_path=out,
            sample_label="HG002.maternal",
            database_id="KS_human_CHM13_v2",
            feature_set_label="region",
            smoothed=True,
        )
        text = out.read_text()
        # Title fields should all appear in the SVG (drawsvg writes
        # Text as literal text nodes).
        assert "HG002.maternal" in text
        assert "KS_human_CHM13_v2" in text
        assert "genome" in text
        assert "region" in text
        assert "smoothed" in text

    def test_no_title_suppresses_title(self, tmp_path: Path) -> None:
        out = tmp_path / "x.svg"
        render_karyotype(
            [self._basic_input()],
            colors={"rA": "#000"},
            mode="genome",
            output_path=out,
            sample_label="HG002",
            database_id="db1",
            feature_set_label="region",
            show_title=False,
        )
        text = out.read_text()
        # The sample label should NOT be present when show_title=False.
        assert "HG002" not in text


class TestLegend:
    def test_legend_lists_drawn_features_only(self, tmp_path: Path) -> None:
        # Two features in the colors map, only one in the binned data.
        # The legend should list only the one actually drawn.
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_h1_a": [(0, 1000, "rA")]},
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#2ca02c", "rB": "#d62728"},  # rB unused
            mode="genome",
            output_path=out,
            show_title=False,
        )
        text = out.read_text()
        # rA should appear as a legend label; rB should not.
        assert text.count("rA") >= 1
        assert "rB" not in text

    def test_title_fits_the_canvas_over_few_chromosomes(self, tmp_path: Path) -> None:
        # A long title over a narrow karyotype must widen the canvas to hold it.
        # SVG viewers let text overflow the canvas, so the SVG looks fine while
        # the rasterised PNG loses both ends -- which is exactly how this shipped:
        # the Arabidopsis genome plot (5 chromosomes, ~104-char title) rendered
        # correctly on one node and clipped on another whose default font was
        # wider. The canvas must fit the title even for a font appreciably wider
        # than the one used here.
        rows = [
            _row(f"Chr{i}_h1_c{i}", chrom=f"Chr{i}", hap="hap1", length=30_000_000)
            for i in range(1, 6)
        ]
        ri = RenderInput(
            map_rows=rows,
            binned_bed={f"Chr{i}_h1_c{i}": [(0, 30_000_000, f"Chr{i}")] for i in range(1, 6)},
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={f"Chr{i}": "#2ca02c" for i in range(1, 6)},
            mode="genome",
            output_path=out,
            sample_label="Col-CEN_v1.2",
            database_id="HKS_arabidopsis_ColCEN",
            feature_set_label="chromosome",
            smoothed=True,
        )
        svg = out.read_text()
        width = float(re.search(r'width="([\d.]+)"', svg).group(1))
        title = re.search(r">([^<]*database[^<]*)<", svg).group(1)
        # 7.9 px/char is above the ~7.79 that clipped on the node where this was
        # found, so the canvas has real headroom rather than only just fitting.
        assert width >= len(title) * 7.9

    def test_every_text_element_declares_a_font_family(self, tmp_path: Path) -> None:
        # Without font-family the renderer picks its own default sans-serif, so
        # the same SVG rasterises to different text widths on different machines
        # -- and the canvas is sized to hold the title, so a wider default clips
        # it. Pinning the family is what makes the output reproducible.
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_h1_a": [(0, 1000, "rA")]},
        )
        out = tmp_path / "x.svg"
        render_karyotype([ri], colors={"rA": "#2ca02c"}, mode="genome", output_path=out)
        svg = out.read_text()
        n_text = svg.count("<text")
        assert n_text > 0
        assert svg.count("font-family") == n_text

    def test_legend_collapses_features_sharing_a_group(self, tmp_path: Path) -> None:
        # The motivating case: the CHM13 cytoband database has 833 features in
        # a handful of colours, and the legend silently truncated to the ~51
        # that fit the canvas. Grouping collapses them to one row per group.
        length = 3_000_000
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=length)],
            binned_bed={
                "chr1_h1_a": [
                    (0, 1_000_000, "p11.1"),
                    (1_000_000, 2_000_000, "q21.3"),
                    (2_000_000, length, "p13.2"),
                ]
            },
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"p11.1": "#000000", "q21.3": "#000000", "p13.2": "#ffffff"},
            legend_groups={"p11.1": "gpos100", "q21.3": "gpos100", "p13.2": "gneg"},
            mode="genome",
            output_path=out,
            show_title=False,
        )
        text = out.read_text()
        # Group labels replace the per-feature labels entirely.
        assert "gpos100" in text
        assert "gneg" in text
        assert "p11.1" not in text
        assert "q21.3" not in text
        assert "p13.2" not in text

    def test_legend_ungrouped_without_groups(self, tmp_path: Path) -> None:
        # A database with no legend_group column (every one shipped before the
        # option existed) must keep its per-feature legend unchanged.
        length = 2_000_000
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=length)],
            binned_bed={"chr1_h1_a": [(0, 1_000_000, "rA"), (1_000_000, length, "rB")]},
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#2ca02c", "rB": "#d62728"},
            mode="genome",
            output_path=out,
            show_title=False,
        )
        text = out.read_text()
        assert "rA" in text
        assert "rB" in text

    def test_outlines_drawn_on_dark_background(self, tmp_path: Path) -> None:
        # Outlines were white-background-only and hardcoded black, so a dark
        # plot got no sequence border -- and a black-filled cytoband then
        # merged into the backdrop. The outline must contrast with the page.
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_h1_a": [(0, 1000, "gpos100")]},
        )
        outs = {}
        for bg in ("white", "black"):
            p = tmp_path / f"{bg}.svg"
            render_karyotype(
                [ri],
                colors={"gpos100": "#000000"},
                mode="genome",
                output_path=p,
                show_title=False,
                background_color=bg,
            )
            outs[bg] = p.read_text()
        # Dark theme strokes white, light theme strokes black -- and neither
        # leaves the sequence unbordered.
        assert 'stroke="#FFFFFF"' in outs["black"]
        assert 'stroke="#000000"' in outs["white"]

    def test_no_hardcoded_black_strokes_remain(self) -> None:
        # The regression this guards: strokes that ignore the theme. Every
        # stroke in the module must derive from outline_color/text_color.
        src = Path(karyotype_mod.__file__).read_text()
        assert 'stroke="black"' not in src

    def test_legend_omits_subpixel_features(self, tmp_path: Path) -> None:
        # A feature whose entire drawn extent is a fraction of a pixel cannot
        # be found in the figure, so a legend row for it only sends the reader
        # hunting for a colour that was never visibly rendered. Mirrors the
        # observed CHM13 case: a 48 bp `categorized` tail on a 248 Mb chr1.
        length = 248_387_298
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=length)],
            binned_bed={
                "chr1_h1_a": [(0, length - 48, "chr1"), (length - 48, length, "categorized")]
            },
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"chr1": "#2ca02c", "categorized": "#B0C4DE"},
            mode="genome",
            output_path=out,
            show_title=False,
        )
        text = out.read_text()
        assert "chr1" in text
        assert "categorized" not in text

    def test_legend_keeps_small_but_visible_features(self, tmp_path: Path) -> None:
        # The filter is on visibility, not smallness: a feature occupying a few
        # percent of a chromosome is easily visible and must keep its row.
        length = 1_000_000
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=length)],
            binned_bed={"chr1_h1_a": [(0, 950_000, "rA"), (950_000, length, "rB")]},
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#2ca02c", "rB": "#d62728"},
            mode="genome",
            output_path=out,
            show_title=False,
        )
        text = out.read_text()
        assert "rA" in text
        assert "rB" in text

    def test_no_legend_suppresses_legend(self, tmp_path: Path) -> None:
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_h1_a": [(0, 1000, "uniquely_named_feature_xyz")]},
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"uniquely_named_feature_xyz": "#000"},
            mode="genome",
            output_path=out,
            show_title=False,
            show_legend=False,
        )
        text = out.read_text()
        # The feature name only appears as a legend label; the
        # rendered rectangles use the colour, not the name. With
        # legend off, the name should be entirely absent.
        assert "uniquely_named_feature_xyz" not in text

    def test_legend_respects_feature_order(self, tmp_path: Path) -> None:
        # Three features drawn; the requested order puts them in
        # reverse alphabetical, which is different from any natural
        # fallback so we can verify the order was honoured. Length
        # is generous so the karyotype is tall enough for all three
        # legend rows to fit (short contigs trigger the truncation
        # safety in real-world too-small renders, irrelevant here).
        ri = RenderInput(
            map_rows=[
                _row("chr1_h1_a", chrom="chr1", hap="hap1", length=30_000_000),
            ],
            binned_bed={
                "chr1_h1_a": [
                    (0, 10_000_000, "aaa"),
                    (10_000_000, 20_000_000, "bbb"),
                    (20_000_000, 30_000_000, "ccc"),
                ],
            },
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"aaa": "#100", "bbb": "#020", "ccc": "#003"},
            mode="genome",
            output_path=out,
            show_title=False,
            feature_order=["ccc", "bbb", "aaa"],  # reversed
        )
        text = out.read_text()
        # The legend labels appear in the SVG in the order they were
        # drawn, which is the feature_order order. We just need
        # ``ccc`` to appear before ``bbb`` before ``aaa`` somewhere
        # after the karyotype draws.
        pos_a = text.find("aaa")
        pos_b = text.find("bbb")
        pos_c = text.find("ccc")
        assert pos_a > 0 and pos_b > 0 and pos_c > 0, "all three labels should be present"
        assert pos_c < pos_b < pos_a

    def test_categorized_pins_to_top_novel_pins_to_bottom(self, tmp_path: Path) -> None:
        # "categorized" is the hierarchy root -- only in the parent
        # column, never in feature_order. "novel" is the k-mer-not-
        # in-index sentinel, never in the hierarchy. Both deserve
        # fixed positions in the legend regardless of the rest of
        # the sort. This test puts a real feature (rA) in the
        # middle and verifies the order is: categorized, rA, novel.
        ri = RenderInput(
            map_rows=[
                _row("chr1_h1_a", chrom="chr1", hap="hap1", length=30_000_000),
            ],
            binned_bed={
                "chr1_h1_a": [
                    (0, 10_000_000, "categorized"),
                    (10_000_000, 20_000_000, "rA"),
                    (20_000_000, 30_000_000, "novel"),
                ],
            },
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"categorized": "#aaa", "rA": "#2ca02c"},  # novel handled by sentinel
            mode="genome",
            output_path=out,
            show_title=False,
            feature_order=["rA"],  # only rA in hierarchy order
        )
        text = out.read_text()
        pos_categorized = text.find("categorized")
        pos_ra = text.find("rA")
        pos_novel = text.find("novel")
        assert pos_categorized > 0 and pos_ra > 0 and pos_novel > 0
        # Visual order in SVG matches draw order: categorized first
        # (top of legend), then rA, then novel (bottom).
        assert pos_categorized < pos_ra < pos_novel

    def test_fallback_sort_pins_chromosomes_top_then_categorized_then_novel(
        self, tmp_path: Path
    ) -> None:
        # When feature_order is None the chr-then-alpha sort applies.
        # Under the HPRC-aligned layout chromosomes (chr*) pin to the
        # very top, then "categorized" (hierarchy root), then any
        # other features alphabetical, then "novel" at the bottom.
        ri = RenderInput(
            map_rows=[
                _row("chr1_h1_a", chrom="chr1", hap="hap1", length=30_000_000),
            ],
            binned_bed={
                "chr1_h1_a": [
                    (0, 10_000_000, "categorized"),
                    (10_000_000, 20_000_000, "chr1"),
                    (20_000_000, 30_000_000, "novel"),
                ],
            },
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"categorized": "#aaa", "chr1": "#1f77b4"},
            mode="genome",
            output_path=out,
            show_title=False,
            # No feature_order -- fallback sort applies.
        )
        text = out.read_text()
        pos_categorized = text.find("categorized")
        # The last "chr1" occurrence is the legend label; the earlier
        # ones are chromosome-column labels above the karyotype.
        pos_chr1_legend = text.rfind("chr1")
        pos_novel = text.find("novel")
        # chr1 (legend) appears before categorized, which appears
        # before novel.
        assert pos_chr1_legend < pos_categorized < pos_novel

    def test_title_uses_database_and_feature_set_nouns(self, tmp_path: Path) -> None:
        ri = RenderInput(
            map_rows=[_row("chr1_h1_a", chrom="chr1", hap="hap1", length=1000)],
            binned_bed={"chr1_h1_a": [(0, 1000, "rA")]},
        )
        out = tmp_path / "x.svg"
        render_karyotype(
            [ri],
            colors={"rA": "#000"},
            mode="genome",
            output_path=out,
            sample_label="HG002",
            database_id="KS_human_CHM13_v2",
            feature_set_label="gene",
        )
        text = out.read_text()
        # Nouns appear next to their tags so the title reads as a
        # sentence rather than a string of unlabelled IDs.
        assert "KS_human_CHM13_v2 database" in text
        assert "genome view" in text
        assert "gene feature set" in text


# --- format conversion --------------------------------------------


def _has_cairo() -> bool:
    """True if cairosvg can import (its native libcairo is present)."""
    try:
        import cairosvg  # noqa: F401

        return True
    except OSError:
        return False


class TestConvertSvg:
    def _make_svg(self, tmp_path: Path) -> Path:
        # A minimal valid SVG cairosvg will accept.
        p = tmp_path / "src.svg"
        p.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">'
            '<rect x="10" y="10" width="80" height="80" fill="#2ca02c"/></svg>'
        )
        return p

    @pytest.mark.skipif(not _has_cairo(), reason="needs native libcairo")
    def test_svg_to_pdf(self, tmp_path: Path) -> None:
        src = self._make_svg(tmp_path)
        dst = tmp_path / "out.pdf"
        convert_svg(src, dst)
        # PDF magic bytes
        assert dst.is_file()
        assert dst.read_bytes()[:4] == b"%PDF"

    @pytest.mark.skipif(not _has_cairo(), reason="needs native libcairo")
    def test_svg_to_png(self, tmp_path: Path) -> None:
        src = self._make_svg(tmp_path)
        dst = tmp_path / "out.png"
        convert_svg(src, dst)
        # PNG magic bytes
        assert dst.is_file()
        assert dst.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"

    def test_svg_to_svg_copies(self, tmp_path: Path) -> None:
        src = self._make_svg(tmp_path)
        dst = tmp_path / "out.svg"
        convert_svg(src, dst)
        assert dst.read_text() == src.read_text()

    def test_unsupported_format_raises(self, tmp_path: Path) -> None:
        src = self._make_svg(tmp_path)
        dst = tmp_path / "out.bmp"
        with pytest.raises(KaryotypeError, match="unsupported output format"):
            convert_svg(src, dst)

    def test_pdf_without_cairo_gives_actionable_error(self, tmp_path: Path) -> None:
        if _has_cairo():
            pytest.skip("libcairo is installed; can't exercise the missing-cairo branch")
        src = self._make_svg(tmp_path)
        dst = tmp_path / "out.pdf"
        with pytest.raises(KaryotypeError, match="libcairo"):
            convert_svg(src, dst)


# --- CLI parsing ---------------------------------------------------


class TestCliSurface:
    def test_help_runs(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["karyotype", "--help"])
        assert result.exit_code == 0
        assert "render" in result.output.lower()

    def test_no_input_fails(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["karyotype"])
        assert result.exit_code != 0

    def test_outdir_and_output_conflict(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n")
        result = cli_runner.invoke(
            main,
            [
                "karyotype",
                "-i",
                str(fa),
                "--outdir",
                str(tmp_path / "a"),
                "--output",
                str(tmp_path / "b.svg"),
            ],
        )
        assert result.exit_code != 0
        assert "--outdir and --output" in result.output

    def test_help_lists_combine_chromosomes(self, cli_runner: CliRunner) -> None:
        result = cli_runner.invoke(main, ["karyotype", "--help"])
        assert result.exit_code == 0
        assert "--combine-chromosomes" in result.output

    def test_negative_scaffold_gap_size_rejected(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
    ) -> None:
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n")
        result = cli_runner.invoke(
            main,
            ["karyotype", "-i", str(fa), "--combine-chromosomes", "--scaffold-gap-size", "-1"],
        )
        assert result.exit_code != 0
        assert "scaffold-gap-size" in result.output

    def test_scaffold_db_with_combine_chromosomes_rejected(
        self,
        cli_runner: CliRunner,
        populated_db_root: Path,
        tmp_path: Path,
    ) -> None:
        """The two flags are mutually exclusive, and the combine path relies on it.

        Combined layout reads the plotted feature set's scaffolded BEDs from
        the layout database, which a separate --scaffold-db does not have. The
        combine path also reads scaffold_map.tsv under the layout database's
        id, so it only agrees with what scaffold_run wrote while the two ids
        are the same. Nothing pinned this rejection; if it were relaxed, the
        first symptom would be "scaffold map not found".
        """
        fa = tmp_path / "x.fa"
        fa.write_text(">a\nACGT\n")
        result = cli_runner.invoke(
            main,
            [
                "karyotype",
                "-i",
                str(fa),
                "--db-root",
                str(populated_db_root),
                "--scaffold-db",
                "some_other_db",
                "--combine-chromosomes",
            ],
        )
        assert result.exit_code != 0
        assert "--scaffold-db cannot be combined with --combine-chromosomes" in result.output


# --- integration tests against the dummy DB ------------------------


def _has(name: str) -> bool:
    return shutil.which(name) is not None


_required = pytest.mark.skipif(
    not (_has("seqtk") and _has("bgzip")),
    reason="needs seqtk and bgzip on PATH",
)


@pytest.fixture
def dummy_assembly_fasta(tmp_path: Path) -> Path:
    seed = "ACGTGCTAGCTAGGCTATCGTAC"
    fa = tmp_path / "assembly.fa"
    fa.write_text(f">seq_for_chr1\n{seed[:21]}\n>seq_for_chr2\n{seed[2:23]}\n")
    return fa


@_required
@pytest.mark.integration
def test_karyotype_genome_mode_against_dummy_db(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "karyotype",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "genome",
            "--bin-size",
            "10",  # tiny -- sequences are 21 bp
            "--min-scaffold-length",
            "1",
            "--feature-set",
            "region",
            "--no-human-chroms",  # dummy db has only chr1, chr2
        ],
    )
    assert result.exit_code == 0, result.output

    svgs = list(out_dir.glob("*.karyotype.svg"))
    assert svgs, f"no SVGs produced; got: {list(out_dir.iterdir())}"
    text = svgs[0].read_text()
    assert "<svg" in text
    assert "</svg>" in text
    # Colour from the dummy db's region set.
    assert "#2ca02c" in text or "#d62728" in text or "#9467bd" in text


@_required
@pytest.mark.integration
def test_karyotype_combine_chromosomes_against_dummy_db(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """--combine-chromosomes cascades scaffold's combine path and renders
    from the combined-chromosome BEDs.

    Verifies the karyotype outputs carry the ``combined_chromosomes`` tag
    and that the combined scaffolded FASTA + AGP side artifacts are
    produced.
    """
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "karyotype",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "genome",
            "--bin-size",
            "10",
            "--min-scaffold-length",
            "1",
            "--feature-set",
            "region",
            "--no-human-chroms",
            "--combine-chromosomes",
        ],
    )
    assert result.exit_code == 0, result.output

    tagged_svgs = list(out_dir.glob("*.combined_chromosomes.karyotype.svg"))
    assert tagged_svgs, f"no combined SVGs; got: {list(out_dir.iterdir())}"
    text = tagged_svgs[0].read_text()
    assert "<svg" in text and "</svg>" in text

    # The combine cascade leaves the combined FASTA + AGP as side artifacts.
    assert list(out_dir.glob("*.scaffolded.combined_chromosomes.fa*"))
    assert list(out_dir.glob("*.scaffolded.combined_chromosomes.agp"))
    # The combined binned BED (keyed by <chrom>_<hap>) was materialised.
    assert list(out_dir.glob("*.scaffolded.combined_chromosomes.binned*.bed.gz"))


@_required
@pytest.mark.integration
def test_karyotype_combine_centromere_mode_against_dummy_db(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """Combine + centromere mode detects centromeres in combined coords
    and writes a combined-tagged centromeres BED."""
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "karyotype",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "centromere",
            "--bin-size",
            "10",
            "--min-scaffold-length",
            "1",
            "--feature-set",
            "region",
            "--no-human-chroms",
            "--combine-chromosomes",
        ],
    )
    assert result.exit_code == 0, result.output
    assert list(out_dir.glob("*.centromere.region.smoothed.combined_chromosomes.karyotype.svg"))
    assert list(out_dir.glob("*.centromeres.combined_chromosomes.bed*"))


@_required
@pytest.mark.integration
def test_karyotype_centromere_mode_against_dummy_db(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "karyotype",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "centromere",
            "--bin-size",
            "10",
            "--min-scaffold-length",
            "1",
            "--feature-set",
            "region",
            "--no-human-chroms",
        ],
    )
    # Centromere mode requires non-empty centromere ranges across all
    # inputs. The dummy db's rA/rB/rC region features all classify as
    # centromere via get_simple_region, so we should get a non-empty SVG.
    assert result.exit_code == 0, result.output
    svgs = list(out_dir.glob("*.centromere.region.smoothed.karyotype.svg"))
    assert svgs


@_required
@pytest.mark.integration
def test_karyotype_no_scaffolding_path_renders(
    cli_runner: CliRunner,
    populated_db_root: Path,
    dummy_assembly_fasta: Path,
    tmp_path: Path,
) -> None:
    """``--no-scaffolding`` skips the full-resolution scaffolded BED but
    still renders a valid SVG via the bin-time map-application path.

    Verifies:
    1. The run exits cleanly.
    2. The SVG is produced and well-formed.
    3. No ``*.smoothed.scaffolded.bed[.gz]`` file is left on disk (the
       whole point of the flag).
    4. The binned-scaffolded BED still exists (the new path produces
       it via bin-then-rewrite).
    """
    out_dir = tmp_path / "out"
    result = cli_runner.invoke(
        main,
        [
            "karyotype",
            "-i",
            f"hap1={dummy_assembly_fasta}",
            "--outdir",
            str(out_dir),
            "--mode",
            "genome",
            "--bin-size",
            "10",
            "--min-scaffold-length",
            "1",
            "--feature-set",
            "region",
            "--no-human-chroms",
            "--no-scaffolding",
        ],
    )
    assert result.exit_code == 0, result.output

    svgs = list(out_dir.glob("*.karyotype.svg"))
    assert svgs, f"no SVGs produced; got: {list(out_dir.iterdir())}"
    text = svgs[0].read_text()
    assert "<svg" in text
    assert "</svg>" in text

    # The full-resolution scaffolded BED must NOT exist.
    scaffolded = list(out_dir.glob("*.region.smoothed.scaffolded.bed*"))
    assert not scaffolded, (
        f"--no-scaffolding should skip writing full-resolution scaffolded BEDs, "
        f"but found: {scaffolded}"
    )

    # The binned-scaffolded BED MUST exist (the bin-time map application
    # path still produces it for the renderer to consume).
    binned_scaffolded = list(out_dir.glob("*.region.smoothed.scaffolded.binned*.bed.gz"))
    assert binned_scaffolded, f"binned scaffolded BED missing; got: {list(out_dir.iterdir())}"

    # And the scaffold_map.tsv must exist (it's how the binner applied
    # rename + flip).
    maps = list(out_dir.glob("*.scaffold_map.tsv"))
    assert maps, "scaffold_map.tsv should still be written when --no-scaffolding"


class TestColorsFilenameTag:
    """The --colors output tag keeps custom-colour renders from clobbering defaults."""

    def test_default_colors_no_tag(self) -> None:
        from karyoscope.core.karyotype_run import _colors_filename_tag

        assert _colors_filename_tag(None) == ""

    def test_custom_colors_tagged_by_stem(self) -> None:
        from pathlib import Path

        from karyoscope.core.karyotype_run import _colors_filename_tag

        assert _colors_filename_tag(Path("/x/colors_chromosome.tsv")) == ".colors_chromosome"
        # A custom file and the default produce different tags -> different filenames.
        assert _colors_filename_tag(Path("/x/my_palette.tsv")) != _colors_filename_tag(None)


# --- scaffolding prerequisite note --------------------------------------


def test_prereq_note_names_the_role_set_pulled_in_by_scaffolding() -> None:
    """Asking for 2 feature sets and seeing annotate report 3 is baffling.

    Scaffolding needs the chromosome- and region-assignment sets to place
    and orient contigs regardless of what is being plotted, so the count
    in the progress output exceeds --feature-set. The note explains which
    extra set appeared and why.
    """
    from karyoscope.core.karyotype_run import _scaffolding_prereq_note

    note = _scaffolding_prereq_note(
        requested=["chromosome", "repeat"],
        scaffold_manifest_roles={
            "chromosome_assignment": "chromosome",
            "region_assignment": "region",
        },
        scaffold_available=["chromosome", "region", "repeat"],
        centromere_fs=None,
        scaffold_db_id=None,
    )
    assert "region" in note
    assert "not rendered" in note
    # The already-requested role set must not be listed as an extra.
    assert "chromosome" not in note


def test_prereq_note_is_empty_when_the_roles_were_requested() -> None:
    """No note when nothing extra is pulled in -- keep the common case clean."""
    from karyoscope.core.karyotype_run import _scaffolding_prereq_note

    assert (
        _scaffolding_prereq_note(
            requested=["chromosome", "region"],
            scaffold_manifest_roles={
                "chromosome_assignment": "chromosome",
                "region_assignment": "region",
            },
            scaffold_available=["chromosome", "region"],
            centromere_fs=None,
            scaffold_db_id=None,
        )
        == ""
    )


def test_prereq_note_includes_the_centromere_set_and_names_a_layout_db() -> None:
    from karyoscope.core.karyotype_run import _scaffolding_prereq_note

    note = _scaffolding_prereq_note(
        requested=["cytoband"],
        scaffold_manifest_roles={
            "chromosome_assignment": "chromosome",
            "region_assignment": "region",
        },
        scaffold_available=["chromosome", "region"],
        centromere_fs="region",
        scaffold_db_id="KS_human_CHM13_v2",
    )
    assert "chromosome, region" in note
    assert "KS_human_CHM13_v2" in note
    # centromere_fs duplicates region_assignment here; it must not repeat.
    assert note.count("region") == 1


def test_prereq_note_stays_silent_when_roles_cannot_be_resolved() -> None:
    """A plot-only database has no roles; that is the cascade's error to report."""
    from karyoscope.core.karyotype_run import _scaffolding_prereq_note

    assert (
        _scaffolding_prereq_note(
            requested=["cytoband"],
            scaffold_manifest_roles={},
            scaffold_available=["cytoband"],
            centromere_fs=None,
            scaffold_db_id=None,
        )
        == ""
    )
