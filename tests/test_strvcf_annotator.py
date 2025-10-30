# #!/usr/bin/env python

# """Tests for `strvcf_annotator` package."""

# import pytest
# import hashlib
# import os
# import strvcf_annotator.strvcf_annotator as strvcf_annotator
# import pysam


# # def test_intersect_vcf_with_str(tmp_path):
# #     from strvcf_annotator import run_annotation_pipeline  # hypothetical main function

# #     input_vcf = "tests/data/input.vcf"
# #     str_bed = "tests/data/str_regions.bed"
# #     expected_hash = "b3d3c0f4a9...abcdef"  # precomputed

# #     output_vcf = tmp_path / "out.vcf"
# #     run_annotation_pipeline(input_vcf, str_bed, output_vcf)

# #     assert file_hash(output_vcf) == expected_hash

# def file_hash(file_path):
#     return hashlib.md5(open(output_vcf,'rb').read()).hexdigest()


  
# @pytest.fixture(
#     params=[
#         "MUTO-INTL.DO253554.SA615942.wgs.20240808.gatk-mutect2.somatic.indel.vcf.gz", # mutect2
#         "MUTO-INTL.DO253205.SA615836.wgs.20240922.sanger-wgs.somatic.indel.annotated.vcf",  # pindel
#         "APGI-AU.DO32825.SA407790.wgs.20210623.gatk-mutect2.somatic.indel.vcf.gz", # with index file
#         "TCGA-DC-6682.vcf" # TCGA tumor only with GangSTR annotation
#     ],
#     hash=[
#         ""
#     ]
# )

# def vcf_in(request):
#     base_dir = os.path.dirname(__file__)
#     vcf_path = os.path.abspath(os.path.join(base_dir, "data", request.param))
#     return pysam.VariantFile(vcf_path)

# class TestMakeModifiedHeader:
#     def test_returns_variant_header(self, vcf_in):
#         new_header = strvcf_annotator.make_modified_header(vcf_in)
#         assert isinstance(new_header, pysam.VariantHeader)

#     def test_new_info_fields_present(self, vcf_in):
#         new_header = strvcf_annotator.make_modified_header(vcf_in)
#         for key in ['RU', 'PERIOD', 'REF', 'PERFECT']:
#             assert key in new_header.info

#     def test_new_format_fields_present(self, vcf_in):
#         new_header = strvcf_annotator.make_modified_header(vcf_in)
#         assert 'REPCN' in new_header.formats

#     def test_field_definitions_correct(self, vcf_in):
#         header = strvcf_annotator.make_modified_header(vcf_in)
        
#         # Check INFO fields
#         assert 'RU' in header.info
#         assert header.info['RU'].number == 1
#         assert header.info['RU'].type == 'String'

#         assert 'PERIOD' in header.info
#         assert header.info['PERIOD'].number == 1
#         assert header.info['PERIOD'].type == 'Integer'

#         assert 'REF' in header.info
#         assert header.info['REF'].number == 1
#         assert header.info['REF'].type == 'Integer'

#         assert 'PERFECT' in header.info
#         assert header.info['PERFECT'].number == 1
#         assert header.info['PERFECT'].type == 'String'
        
#         # Check FORMAT field
#         assert 'REPCN' in header.formats
#         assert header.formats['REPCN'].number == 2
#         assert header.formats['REPCN'].type == 'Integer'
    
#     def test_modified_header_removes_no_unexpected_lines(self, vcf_in):
#         original_header = str(vcf_in.header).strip().splitlines()
#         modified_header = str(strvcf_annotator.make_modified_header(vcf_in)).strip().splitlines()

#         removed_lines = set(original_header) - set(modified_header)

#         # Only flag lines that are NOT ones we explicitly replaced
#         replaced_prefixes = {
#             '##INFO=<ID=RU',
#             '##INFO=<ID=PERIOD',
#             '##INFO=<ID=REF',
#             '##INFO=<ID=PERFECT',
#             '##FORMAT=<ID=REPCN',
#         }

#         # Keep only removed lines that are not expected to be removed
#         unexpected_removed = [
#             line for line in removed_lines
#             if not any(line.startswith(prefix) for prefix in replaced_prefixes)
#         ]

#         if unexpected_removed:
#             removed_text = "\n".join(unexpected_removed)
#             assert False, f"Unexpected header lines removed:\n{removed_text}"



#     def test_modified_header_only_adds_expected_lines(self, vcf_in):
#         original_header = str(vcf_in.header).strip().splitlines()
#         modified_header = str(strvcf_annotator.make_modified_header(vcf_in)).strip().splitlines()

#         expected_new_info = {
#             '##INFO=<ID=RU',
#             '##INFO=<ID=PERIOD',
#             '##INFO=<ID=REF',
#             '##INFO=<ID=PERFECT'
#         }
#         expected_new_format = {
#             '##FORMAT=<ID=REPCN'
#         }

#         added_lines = set(modified_header) - set(original_header)
#         expected_prefixes = expected_new_info.union(expected_new_format)

#         unexpected_lines = [
#             line for line in added_lines
#             if not any(line.startswith(prefix) for prefix in expected_prefixes)
#         ]

#         assert not unexpected_lines, f"Unexpected new lines in header: {unexpected_lines}"

#     def test_modified_header_includes_all_expected_fields(self, vcf_in):
#         modified_header = str(strvcf_annotator.make_modified_header(vcf_in)).strip().splitlines()

#         expected_new_lines = {
#             '##INFO=<ID=RU',
#             '##INFO=<ID=PERIOD',
#             '##INFO=<ID=REF',
#             '##INFO=<ID=PERFECT',
#             '##FORMAT=<ID=REPCN'
#         }

#         for prefix in expected_new_lines:
#             assert any(line.startswith(prefix) for line in modified_header), f"Missing expected header line: {prefix}"

# class TestProcessVcf:
#     def test_process_vcf(self, request):
#     base_dir = os.path.dirname(__file__)  # directory where this test file lives

#     input_vcf = os.path.abspath(os.path.join("data", "MUTO-INTL.DO253554.SA615942.wgs.20240808.gatk-mutect2.somatic.indel.vcf.gz"))
#     str_bed = os.path.abspath(os.path.join("data", "GRCh38_repeats.bed"))
#     output_vcf = os.path.abspath(os.path.join("output", "MUTO-INTL.DO253554.SA615942.wgs.20240808.gatk-mutect2.somatic.indel.vcf"))
#     expected_hash = "269f36179979776f46c2d23b77be512b"

#     os.makedirs(os.path.dirname(output_vcf), exist_ok=True)  # ensure output dir exists
#     str_df = strvcf_annotator.load_str_reference(str_bed)
#     strvcf_annotator.process_vcf(input_vcf, str_df, output_vcf)

#     assert file_hash(output_vcf) == expected_hash 
# # class TestBuildNewRecord:

# #     @classmethod
# #     def setup_class(cls):
# #         # Create a simple header with the required fields
# #         base_header = pysam.VariantHeader()
# #         base_header.contigs.add("1", length=1000)
# #         base_header.info.add("DP", 1, "Integer", "Read depth")
# #         base_header.formats.add("GT", 1, "String", "Genotype")
# #         base_header.add_sample("S1")

# #         cls.header = make_modified_header(pysam.VariantFile(header=base_header))

# #         # Define a minimal STR region
# #         cls.str_row = pd.Series({
# #             'CHROM': '1',
# #             'START': 100,
# #             'END': 111,
# #             'RU': 'CA',
# #             'PERIOD': 2
# #         })

# #         # Create a variant: CA deletion inside STR region
# #         variant_vcf = pysam.VariantFile(header=base_header, mode='w')
# #         rec = variant_vcf.new_record(
# #             contig="1", start=104, stop=106,
# #             alleles=("CA", "C"),
# #             filter=["PASS"]
# #         )
# #         rec.samples["S1"]["GT"] = (0, 1)
# #         cls.original_record = rec

# #     def test_returns_variant_record(self, vcf_in):
# #         new_rec = build_new_record(self.original_record, self.str_row, self.header)
# #         assert isinstance(new_rec, pysam.VariantRecord)

# #     def test_alleles_are_repeat_sequences(self, vcf_in):
# #         new_rec = build_new_record(self.original_record, self.str_row, self.header)
# #         assert new_rec.ref.startswith("CA")
# #         assert new_rec.alts[0] != "C"  # Should be a repeat string, not just raw ALT from VCF

# #     def test_info_fields_are_set_correctly(self, vcf_in):
# #         new_rec = build_new_record(self.original_record, self.str_row, self.header)
# #         assert new_rec.info["RU"] == "CA"
# #         assert new_rec.info["PERIOD"] == 2
# #         assert isinstance(new_rec.info["REF"], int)
# #         assert new_rec.info["PERFECT"] in {"TRUE", "FALSE"}

# #     def test_format_fields_are_set_correctly(self, vcf_in):
# #         new_rec = build_new_record(self.original_record, self.str_row, self.header)
# #         assert "REPCN" in new_rec.samples["S1"]
# #         alleles = new_rec.samples["S1"]["REPCN"]
# #         assert isinstance(alleles, tuple)
# #         assert all(isinstance(a, int) for a in alleles)

# #     def test_position_and_span_correct(self, vcf_in):
# #         new_rec = build_new_record(self.original_record, self.str_row, self.header)
# #         assert new_rec.start == self.str_row["START"] - 1
# #         assert new_rec.stop == self.str_row["END"]