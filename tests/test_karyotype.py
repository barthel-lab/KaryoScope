"""Unit + integration tests for :mod:`karyoscope.core.karyotype`."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import ClassVar

import pytest
from click.testing import CliRunner

from karyoscope.cli import main
from karyoscope.core.io.scaffold_map import MapRow
from karyoscope.core.karyotype import (
    PREDEFINED_SEX_SYSTEMS,
    RenderInput,
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
