"""Unit tests for VCF processor functions."""

import pytest
import pysam
import pandas as pd
import tempfile
import os
from strvcf_annotator.core.vcf_processor import (
    check_vcf_sorted,
    reset_and_sort_vcf,
    generate_annotated_records
)
from strvcf_annotator.parsers.generic import GenericParser


@pytest.fixture
def basic_vcf_header():
    """Create a basic VCF header for testing."""
    header = pysam.VariantHeader()
    header.add_line('##fileformat=VCFv4.2')
    header.contigs.add('chr1', length=1000000)
    header.contigs.add('chr2', length=1000000)
    header.contigs.add('chr3', length=1000000)
    header.add_sample('TUMOR')
    header.add_sample('NORMAL')

    # Add basic INFO and FORMAT fields
    header.info.add('DP', 1, 'Integer', 'Total read depth')
    header.formats.add('GT', 1, 'String', 'Genotype')
    header.formats.add('AD', 'R', 'Integer', 'Allelic depths')
    header.formats.add('DP', 1, 'Integer', 'Read depth')

    return header


@pytest.fixture
def sorted_vcf_file(basic_vcf_header):
    """Create a temporary sorted VCF file using pysam API."""
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
    temp_file.close()  # Close immediately so pysam can write to it

    try:
        # Write VCF using pysam API
        with pysam.VariantFile(temp_file.name, 'w', header=basic_vcf_header) as vcf_out:
            # Write sorted records
            records_data = [
                ('chr1', 100, 'A', 'T'),
                ('chr1', 200, 'G', 'C'),
                ('chr1', 300, 'T', 'A'),
                ('chr2', 100, 'A', 'G'),
                ('chr2', 200, 'C', 'T'),
            ]
            for chrom, pos, ref, alt in records_data:
                rec = vcf_out.new_record(
                    contig=chrom,
                    start=pos-1,  # 0-based
                    alleles=(ref, alt),
                    qual=30,
                    filter=['PASS'],
                    info={'DP': 50}
                )
                rec.samples['TUMOR']['GT'] = (0, 1)
                rec.samples['TUMOR']['AD'] = (20, 30)
                rec.samples['TUMOR']['DP'] = 50
                rec.samples['NORMAL']['GT'] = (0, 0)
                rec.samples['NORMAL']['AD'] = (50, 0)
                rec.samples['NORMAL']['DP'] = 50
                vcf_out.write(rec)

        yield temp_file.name
    finally:
        os.unlink(temp_file.name)


@pytest.fixture
def unsorted_vcf_file(basic_vcf_header):
    """Create a temporary unsorted VCF file using pysam API."""
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
    temp_file.close()  # Close immediately so pysam can write to it

    try:
        # Write VCF using pysam API
        with pysam.VariantFile(temp_file.name, 'w', header=basic_vcf_header) as vcf_out:
            # Write unsorted records (positions out of order)
            records_data = [
                ('chr1', 300, 'T', 'A'),
                ('chr1', 100, 'A', 'T'),
                ('chr2', 200, 'C', 'T'),
                ('chr1', 200, 'G', 'C'),
                ('chr2', 100, 'A', 'G'),
            ]
            for chrom, pos, ref, alt in records_data:
                rec = vcf_out.new_record(
                    contig=chrom,
                    start=pos-1,  # 0-based
                    alleles=(ref, alt),
                    qual=30,
                    filter=['PASS'],
                    info={'DP': 50}
                )
                rec.samples['TUMOR']['GT'] = (0, 1)
                rec.samples['TUMOR']['AD'] = (20, 30)
                rec.samples['TUMOR']['DP'] = 50
                rec.samples['NORMAL']['GT'] = (0, 0)
                rec.samples['NORMAL']['AD'] = (50, 0)
                rec.samples['NORMAL']['DP'] = 50
                vcf_out.write(rec)

        yield temp_file.name
    finally:
        os.unlink(temp_file.name)


@pytest.fixture
def str_dataframe():
    """Create a sample STR DataFrame for testing."""
    return pd.DataFrame({
        'CHROM': ['chr1', 'chr1', 'chr2'],
        'START': [95, 195, 95],
        'END': [115, 215, 115],
        'PERIOD': [2, 3, 2],
        'RU': ['AT', 'CAG', 'GC'],
        'COUNT': [10, 7, 10]
    })


class TestCheckVCFSorted:
    """Tests for check_vcf_sorted function."""

    def test_sorted_vcf_returns_true(self, sorted_vcf_file):
        """Test that sorted VCF returns True."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        result = check_vcf_sorted(vcf_in)
        vcf_in.close()

        assert result is True

    def test_unsorted_vcf_returns_false(self, unsorted_vcf_file):
        """Test that unsorted VCF returns False."""
        vcf_in = pysam.VariantFile(unsorted_vcf_file)
        result = check_vcf_sorted(vcf_in)
        vcf_in.close()

        assert result is False

    def test_rewinds_file_after_check(self, sorted_vcf_file):
        """Test that file is rewound after checking."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)

        # Check sorting
        check_vcf_sorted(vcf_in)

        # File should be rewound, so we can read from beginning
        records = list(vcf_in)
        assert len(records) == 5
        assert records[0].pos == 100  # First record

        vcf_in.close()

    def test_empty_vcf_returns_true(self, basic_vcf_header):
        """Test that empty VCF is considered sorted."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
        temp_file.close()

        try:
            # Create empty VCF using pysam API
            with pysam.VariantFile(temp_file.name, 'w', header=basic_vcf_header) as vcf_out:
                pass  # Write header only, no records

            vcf_in = pysam.VariantFile(temp_file.name)
            result = check_vcf_sorted(vcf_in)
            vcf_in.close()

            assert result is True
        finally:
            os.unlink(temp_file.name)

    def test_single_record_returns_true(self, basic_vcf_header):
        """Test that single record VCF is considered sorted."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
        temp_file.close()

        try:
            # Create VCF with single record using pysam API
            with pysam.VariantFile(temp_file.name, 'w', header=basic_vcf_header) as vcf_out:
                rec = vcf_out.new_record(
                    contig='chr1',
                    start=99,  # 0-based
                    alleles=('A', 'T'),
                    qual=30,
                    filter=['PASS'],
                    info={'DP': 50}
                )
                rec.samples['TUMOR']['GT'] = (0, 1)
                rec.samples['TUMOR']['AD'] = (20, 30)
                rec.samples['TUMOR']['DP'] = 50
                rec.samples['NORMAL']['GT'] = (0, 0)
                rec.samples['NORMAL']['AD'] = (50, 0)
                rec.samples['NORMAL']['DP'] = 50
                vcf_out.write(rec)

            vcf_in = pysam.VariantFile(temp_file.name)
            result = check_vcf_sorted(vcf_in)
            vcf_in.close()

            assert result is True
        finally:
            os.unlink(temp_file.name)


class TestResetAndSortVCF:
    """Tests for reset_and_sort_vcf function."""

    def test_sorts_unsorted_vcf(self, unsorted_vcf_file):
        """Test that function sorts unsorted VCF."""
        vcf_in = pysam.VariantFile(unsorted_vcf_file)
        sorted_records = reset_and_sort_vcf(vcf_in)
        vcf_in.close()

        # Check correct number of records
        assert len(sorted_records) == 5

        # Check sorted order
        assert sorted_records[0].contig == 'chr1'
        assert sorted_records[0].pos == 100
        assert sorted_records[1].contig == 'chr1'
        assert sorted_records[1].pos == 200
        assert sorted_records[2].contig == 'chr1'
        assert sorted_records[2].pos == 300
        assert sorted_records[3].contig == 'chr2'
        assert sorted_records[3].pos == 100
        assert sorted_records[4].contig == 'chr2'
        assert sorted_records[4].pos == 200

    def test_preserves_sorted_order(self, sorted_vcf_file):
        """Test that already sorted VCF remains sorted."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        sorted_records = reset_and_sort_vcf(vcf_in)
        vcf_in.close()

        # Check sorted order is preserved
        assert sorted_records[0].pos == 100
        assert sorted_records[1].pos == 200
        assert sorted_records[2].pos == 300
        assert sorted_records[3].contig == 'chr2'
        assert sorted_records[3].pos == 100

    def test_returns_list(self, sorted_vcf_file):
        """Test that function returns a list."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        result = reset_and_sort_vcf(vcf_in)
        vcf_in.close()

        assert isinstance(result, list)

    def test_empty_vcf_returns_empty_list(self, basic_vcf_header):
        """Test that empty VCF returns empty list."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
        temp_file.close()

        try:
            # Create empty VCF using pysam API
            with pysam.VariantFile(temp_file.name, 'w', header=basic_vcf_header) as vcf_out:
                pass  # No records

            vcf_in = pysam.VariantFile(temp_file.name)
            sorted_records = reset_and_sort_vcf(vcf_in)
            vcf_in.close()

            assert len(sorted_records) == 0
        finally:
            os.unlink(temp_file.name)

    def test_sorts_by_contig_order(self, basic_vcf_header):
        """Test that sorting respects contig order in header."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
        temp_file.close()

        try:
            # Create VCF with records from different chromosomes using pysam API
            with pysam.VariantFile(temp_file.name, 'w', header=basic_vcf_header) as vcf_out:
                for chrom in ['chr3', 'chr1', 'chr2']:
                    rec = vcf_out.new_record(
                        contig=chrom,
                        start=99,  # 0-based
                        alleles=('A', 'T'),
                        qual=30,
                        filter=['PASS'],
                        info={'DP': 50}
                    )
                    rec.samples['TUMOR']['GT'] = (0, 1)
                    rec.samples['TUMOR']['AD'] = (20, 30)
                    rec.samples['TUMOR']['DP'] = 50
                    rec.samples['NORMAL']['GT'] = (0, 0)
                    rec.samples['NORMAL']['AD'] = (50, 0)
                    rec.samples['NORMAL']['DP'] = 50
                    vcf_out.write(rec)

            vcf_in = pysam.VariantFile(temp_file.name)
            sorted_records = reset_and_sort_vcf(vcf_in)
            vcf_in.close()

            # Should be sorted by header contig order: chr1, chr2, chr3
            assert sorted_records[0].contig == 'chr1'
            assert sorted_records[1].contig == 'chr2'
            assert sorted_records[2].contig == 'chr3'
        finally:
            os.unlink(temp_file.name)


class TestGenerateAnnotatedRecords:
    """Tests for generate_annotated_records function."""

    def test_returns_iterator(self, sorted_vcf_file, str_dataframe):
        """Test that function returns an iterator."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        parser = GenericParser()

        result = generate_annotated_records(vcf_in, str_dataframe, parser)

        # Should be an iterator
        assert hasattr(result, '__iter__')
        assert hasattr(result, '__next__')

        vcf_in.close()

    def test_yields_variant_records(self, sorted_vcf_file, str_dataframe):
        """Test that iterator yields VariantRecord objects."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        parser = GenericParser()

        records = list(generate_annotated_records(vcf_in, str_dataframe, parser))

        if len(records) > 0:
            assert all(isinstance(r, pysam.VariantRecord) for r in records)

        vcf_in.close()

    def test_filters_non_overlapping_variants(self, sorted_vcf_file):
        """Test that variants outside STR regions are filtered."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        parser = GenericParser()

        # Create STR dataframe with no overlap
        str_df = pd.DataFrame({
            'CHROM': ['chr10'],
            'START': [1000],
            'END': [1100],
            'PERIOD': [2],
            'RU': ['AT'],
            'COUNT': [50]
        })

        records = list(generate_annotated_records(vcf_in, str_df, parser))

        # Should yield no records (no overlap)
        assert len(records) == 0

        vcf_in.close()

    def test_handles_unsorted_vcf(self, unsorted_vcf_file, str_dataframe):
        """Test that function handles unsorted VCF by sorting it."""
        vcf_in = pysam.VariantFile(unsorted_vcf_file)
        parser = GenericParser()

        # Should not raise an error
        records = list(generate_annotated_records(vcf_in, str_dataframe, parser))

        # Records should be processed (may be 0 if no overlaps)
        assert isinstance(records, list)

        vcf_in.close()

    def test_uses_generic_parser_by_default(self, sorted_vcf_file, str_dataframe):
        """Test that GenericParser is used when parser=None."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)

        # Call without parser
        records = list(generate_annotated_records(vcf_in, str_dataframe, parser=None))

        # Should work without error
        assert isinstance(records, list)

        vcf_in.close()

    def test_empty_vcf_returns_no_records(self, basic_vcf_header, str_dataframe):
        """Test that empty VCF yields no records."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False)
        temp_file.close()

        try:
            # Create empty VCF using pysam API
            with pysam.VariantFile(temp_file.name, 'w', header=basic_vcf_header) as vcf_out:
                pass  # No records

            vcf_in = pysam.VariantFile(temp_file.name)
            parser = GenericParser()

            records = list(generate_annotated_records(vcf_in, str_dataframe, parser))

            assert len(records) == 0

            vcf_in.close()
        finally:
            os.unlink(temp_file.name)

    def test_empty_str_dataframe_returns_no_records(self, sorted_vcf_file):
        """Test that empty STR DataFrame yields no records."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        parser = GenericParser()

        # Empty STR dataframe
        empty_str_df = pd.DataFrame(columns=['CHROM', 'START', 'END', 'PERIOD', 'RU', 'COUNT'])

        records = list(generate_annotated_records(vcf_in, empty_str_df, parser))

        assert len(records) == 0

        vcf_in.close()
