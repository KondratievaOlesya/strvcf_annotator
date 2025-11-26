"""Unit tests for validation utilities."""

import tempfile
from pathlib import Path

import pytest

from strvcf_annotator.utils.validation import (
    ValidationError,
    validate_bed_file,
    validate_directory_path,
    validate_file_path,
    validate_str_bed_file,
    validate_vcf_file,
)


class TestValidateFilePath:
    """Test suite for validate_file_path."""

    def test_valid_existing_file(self):
        """Test validation of existing file."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_path = f.name

        try:
            result = validate_file_path(temp_path, must_exist=True)
            assert isinstance(result, Path)
            assert result.exists()
        finally:
            Path(temp_path).unlink()

    def test_nonexistent_file_required(self):
        """Test validation fails for nonexistent file when required."""
        with pytest.raises(ValidationError, match="does not exist"):
            validate_file_path('/nonexistent/file.txt', must_exist=True)

    def test_nonexistent_file_not_required(self):
        """Test validation passes for nonexistent file when not required."""
        result = validate_file_path('/nonexistent/file.txt', must_exist=False)
        assert isinstance(result, Path)

    def test_empty_path(self):
        """Test validation fails for empty path."""
        with pytest.raises(ValidationError, match="cannot be empty"):
            validate_file_path('', must_exist=False)


class TestValidateDirectoryPath:
    """Test suite for validate_directory_path."""

    def test_valid_existing_directory(self):
        """Test validation of existing directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            result = validate_directory_path(temp_dir, must_exist=True)
            assert isinstance(result, Path)
            assert result.is_dir()

    def test_create_directory(self):
        """Test directory creation."""
        with tempfile.TemporaryDirectory() as temp_dir:
            new_dir = Path(temp_dir) / 'new_subdir'
            result = validate_directory_path(str(new_dir), must_exist=False, create=True)

            assert result.exists()
            assert result.is_dir()

    def test_nonexistent_directory_required(self):
        """Test validation fails for nonexistent directory when required."""
        with pytest.raises(ValidationError, match="does not exist"):
            validate_directory_path('/nonexistent/dir', must_exist=True, create=False)

    def test_empty_path(self):
        """Test validation fails for empty path."""
        with pytest.raises(ValidationError, match="cannot be empty"):
            validate_directory_path('', must_exist=False)


class TestValidateVCFFile:
    """Test suite for validate_vcf_file."""

    def test_valid_vcf(self):
        """Test validation of valid VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            temp_path = f.name

        try:
            result = validate_vcf_file(temp_path)
            assert result is True
        finally:
            Path(temp_path).unlink()

    def test_vcf_no_samples(self):
        """Test validation fails for VCF without samples."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\t.\t.
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            temp_path = f.name

        try:
            with pytest.raises(ValidationError, match="no samples"):
                validate_vcf_file(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_invalid_vcf(self):
        """Test validation fails for invalid VCF."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write("This is not a VCF file")
            temp_path = f.name

        try:
            with pytest.raises(ValidationError, match="Invalid VCF"):
                validate_vcf_file(temp_path)
        finally:
            Path(temp_path).unlink()


class TestValidateBEDFile:
    """Test suite for validate_bed_file."""

    def test_valid_bed(self):
        """Test validation of valid BED file."""
        bed_content = "chr1\t100\t200\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write(bed_content)
            temp_path = f.name

        try:
            result = validate_bed_file(temp_path)
            assert result is True
        finally:
            Path(temp_path).unlink()

    def test_bed_too_few_columns(self):
        """Test validation fails for BED with too few columns."""
        bed_content = "chr1\t100\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write(bed_content)
            temp_path = f.name

        try:
            with pytest.raises(ValidationError, match="fewer than 3 columns"):
                validate_bed_file(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_bed_non_numeric_coordinates(self):
        """Test validation fails for non-numeric coordinates."""
        bed_content = "chr1\tabc\t200\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write(bed_content)
            temp_path = f.name

        try:
            with pytest.raises(ValidationError, match="not numeric"):
                validate_bed_file(temp_path)
        finally:
            Path(temp_path).unlink()


class TestValidateSTRBEDFile:
    """Test suite for validate_str_bed_file."""

    def test_valid_str_bed(self):
        """Test validation of valid STR BED file."""
        bed_content = "chr1\t100\t115\t3\tCAG\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write(bed_content)
            temp_path = f.name

        try:
            result = validate_str_bed_file(temp_path)
            assert result is True
        finally:
            Path(temp_path).unlink()

    def test_str_bed_too_few_columns(self):
        """Test validation fails for STR BED with too few columns."""
        bed_content = "chr1\t100\t115\t3\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write(bed_content)
            temp_path = f.name

        try:
            with pytest.raises(ValidationError, match="at least 5 columns"):
                validate_str_bed_file(temp_path)
        finally:
            Path(temp_path).unlink()

    def test_str_bed_non_numeric_period(self):
        """Test validation fails for non-numeric PERIOD."""
        bed_content = "chr1\t100\t115\tabc\tCAG\n"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write(bed_content)
            temp_path = f.name

        try:
            with pytest.raises(ValidationError, match="must be numeric"):
                validate_str_bed_file(temp_path)
        finally:
            Path(temp_path).unlink()
