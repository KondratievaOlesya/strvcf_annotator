"""Unit tests for STR reference management."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

from strvcf_annotator.core.str_reference import (
    find_overlapping_str,
    get_str_at_position,
    load_str_reference,
)


class TestLoadSTRReference:
    """Test suite for load_str_reference."""

    @pytest.fixture
    def temp_bed_file(self):
        """Create temporary BED file."""
        content = """chr1\t100\t115\t3\tCAG
chr1\t200\t212\t4\tATCG
chr2\t300\t318\t3\tGAT"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            temp_path = f.name

        yield temp_path

        # Cleanup
        Path(temp_path).unlink()

    def test_load_basic(self, temp_bed_file):
        """Test basic BED file loading."""
        df = load_str_reference(temp_bed_file)

        assert len(df) == 3
        assert list(df.columns) == ["CHROM", "START", "END", "PERIOD", "RU", "COUNT"]

    def test_coordinate_conversion(self, temp_bed_file):
        """Test BED to VCF coordinate conversion."""
        df = load_str_reference(temp_bed_file)

        # BED START 100 should become VCF START 101
        assert df.iloc[0]["START"] == 101
        # BED END should remain the same
        assert df.iloc[0]["END"] == 115

    def test_count_calculation(self, temp_bed_file):
        """Test repeat count calculation."""
        df = load_str_reference(temp_bed_file)

        # (115 - 101 + 1) / 3 = 5
        assert df.iloc[0]["COUNT"] == 5.0

    def test_sorting(self, temp_bed_file):
        """Test that output is sorted."""
        df = load_str_reference(temp_bed_file)

        # Check chromosomes are sorted
        chroms = df["CHROM"].tolist()
        assert chroms == sorted(chroms)

        # Check positions within chromosome are sorted
        for chrom in df["CHROM"].unique():
            chrom_df = df[df["CHROM"] == chrom]
            positions = chrom_df["START"].tolist()
            assert positions == sorted(positions)


class TestFindOverlappingSTR:
    """Test suite for find_overlapping_str."""

    @pytest.fixture
    def str_df(self):
        """Create sample STR DataFrame."""
        data = {
            "CHROM": ["chr1", "chr1", "chr2"],
            "START": [101, 201, 301],
            "END": [115, 212, 318],
            "PERIOD": [3, 4, 3],
            "RU": ["CAG", "ATCG", "GAT"],
            "COUNT": [5.0, 3.0, 6.0],
        }
        return pd.DataFrame(data)

    def test_exact_overlap(self, str_df):
        """Test exact position overlap."""
        result = find_overlapping_str(str_df, "chr1", 101, 115)

        assert result is not None
        assert result["START"] == 101
        assert result["RU"] == "CAG"

    def test_partial_overlap(self, str_df):
        """Test partial overlap."""
        result = find_overlapping_str(str_df, "chr1", 105, 110)

        assert result is not None
        assert result["START"] == 101

    def test_no_overlap(self, str_df):
        """Test no overlap."""
        result = find_overlapping_str(str_df, "chr1", 150, 160)

        assert result is None

    def test_wrong_chromosome(self, str_df):
        """Test wrong chromosome."""
        result = find_overlapping_str(str_df, "chr3", 101, 115)

        assert result is None

    def test_variant_extends_beyond(self, str_df):
        """Test variant extending beyond STR."""
        result = find_overlapping_str(str_df, "chr1", 110, 120)

        assert result is not None
        assert result["START"] == 101


class TestGetSTRAtPosition:
    """Test suite for get_str_at_position."""

    @pytest.fixture
    def str_df(self):
        """Create sample STR DataFrame."""
        data = {
            "CHROM": ["chr1", "chr1", "chr2"],
            "START": [101, 201, 301],
            "END": [115, 212, 318],
            "PERIOD": [3, 4, 3],
            "RU": ["CAG", "ATCG", "GAT"],
            "COUNT": [5.0, 3.0, 6.0],
        }
        return pd.DataFrame(data)

    def test_position_in_str(self, str_df):
        """Test position within STR."""
        result = get_str_at_position(str_df, "chr1", 105)

        assert result is not None
        assert result["RU"] == "CAG"

    def test_position_outside_str(self, str_df):
        """Test position outside STR."""
        result = get_str_at_position(str_df, "chr1", 150)

        assert result is None

    def test_position_at_boundary(self, str_df):
        """Test position at STR boundary."""
        result = get_str_at_position(str_df, "chr1", 101)

        assert result is not None
        assert result["RU"] == "CAG"
