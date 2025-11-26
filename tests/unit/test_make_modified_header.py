import os

import pysam
import pytest

from strvcf_annotator.core.annotation import make_modified_header


@pytest.fixture(
    params=[
        "test.vcf.gz",
        "pindel_header.vcf",  # pindel
        "mutec2_indel.vcf.gz",  # with index file
        "TCGA-DC-6682.vcf",  # TCGA tumor only with GangSTR annotation
    ]
)
def vcf_in(request, data_dir):
    vcf_path = os.path.abspath(os.path.join(data_dir, request.param))
    return pysam.VariantFile(vcf_path)


class TestMakeModifiedHeader:
    def test_returns_variant_header(self, vcf_in):
        new_header = make_modified_header(vcf_in)
        assert isinstance(new_header, pysam.VariantHeader)

    def test_new_info_fields_present(self, vcf_in):
        new_header = make_modified_header(vcf_in)
        for key in ["RU", "PERIOD", "REF", "PERFECT"]:
            assert key in new_header.info

    def test_new_format_fields_present(self, vcf_in):
        new_header = make_modified_header(vcf_in)
        assert "REPCN" in new_header.formats

    def test_field_definitions_correct(self, vcf_in):
        header = make_modified_header(vcf_in)

        # Check INFO fields
        assert "RU" in header.info
        assert header.info["RU"].number == 1
        assert header.info["RU"].type == "String"

        assert "PERIOD" in header.info
        assert header.info["PERIOD"].number == 1
        assert header.info["PERIOD"].type == "Integer"

        assert "REF" in header.info
        assert header.info["REF"].number == 1
        assert header.info["REF"].type == "Integer"

        assert "PERFECT" in header.info
        assert header.info["PERFECT"].number == 1
        assert header.info["PERFECT"].type == "String"

        # Check FORMAT field
        assert "REPCN" in header.formats
        assert header.formats["REPCN"].number == 2
        assert header.formats["REPCN"].type == "Integer"

    def test_modified_header_removes_no_unexpected_lines(self, vcf_in):
        original_header = str(vcf_in.header).strip().splitlines()
        modified_header = str(make_modified_header(vcf_in)).strip().splitlines()

        removed_lines = set(original_header) - set(modified_header)

        # Only flag lines that are NOT ones we explicitly replaced
        replaced_prefixes = {
            "##INFO=<ID=RU",
            "##INFO=<ID=PERIOD",
            "##INFO=<ID=REF",
            "##INFO=<ID=PERFECT",
            "##FORMAT=<ID=REPCN",
        }

        # Keep only removed lines that are not expected to be removed
        unexpected_removed = [
            line
            for line in removed_lines
            if not any(line.startswith(prefix) for prefix in replaced_prefixes)
        ]

        if unexpected_removed:
            removed_text = "\n".join(unexpected_removed)
            raise AssertionError(f"Unexpected header lines removed:\n{removed_text}")

    def test_modified_header_only_adds_expected_lines(self, vcf_in):
        original_header = str(vcf_in.header).strip().splitlines()
        modified_header = str(make_modified_header(vcf_in)).strip().splitlines()

        expected_new_info = {
            "##INFO=<ID=RU",
            "##INFO=<ID=PERIOD",
            "##INFO=<ID=REF",
            "##INFO=<ID=PERFECT",
        }
        expected_new_format = {"##FORMAT=<ID=REPCN"}

        # Structural VCF header lines that are mandatory per VCF specification
        # or may be added by pysam when creating a new VariantHeader
        structural_prefixes = {
            "##fileformat",  # Mandatory first line in VCF spec
            "##contig",  # Contig definitions
            "##FILTER",  # Filter definitions
            "##reference",  # Reference genome
            "#CHROM",  # Column header line
        }

        added_lines = set(modified_header) - set(original_header)
        expected_prefixes = expected_new_info.union(expected_new_format).union(structural_prefixes)

        unexpected_lines = [
            line
            for line in added_lines
            if not any(line.startswith(prefix) for prefix in expected_prefixes)
        ]

        assert not unexpected_lines, f"Unexpected new lines in header: {unexpected_lines}"

    def test_modified_header_includes_all_expected_fields(self, vcf_in):
        modified_header = str(make_modified_header(vcf_in)).strip().splitlines()

        expected_new_lines = {
            "##INFO=<ID=RU",
            "##INFO=<ID=PERIOD",
            "##INFO=<ID=REF",
            "##INFO=<ID=PERFECT",
            "##FORMAT=<ID=REPCN",
        }

        for prefix in expected_new_lines:
            assert any(line.startswith(prefix) for line in modified_header), (
                f"Missing expected header line: {prefix}"
            )
