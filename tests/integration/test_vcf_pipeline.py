import hashlib
import os
import pytest

def file_hash(path):
    """Calculate MD5 hash of a file."""
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

@pytest.fixture(
    params=[
        (
            "MUTO-INTL.DO253554.SA615942.wgs.20240808.gatk-mutect2.somatic.indel.vcf.gz",  
            "269f36179979776f46c2d23b77be512b"
        ),
        (
            "MUTO-INTL.DO253205.SA615836.wgs.20240922.sanger-wgs.somatic.indel.annotated.vcf",  
            "abcdef1234567890abcdef1234567890"
        ),
        (
            "APGI-AU.DO32825.SA407790.wgs.20210623.gatk-mutect2.somatic.indel.vcf.gz",  
            "fedcba0987654321fedcba0987654321"
        ),
        (
            "TCGA-DC-6682.vcf",  
            "11112222333344445555666677778888"
        )
    ]
)
def vcf_case(request, data_dir):
    input_vcf = os.path.abspath(os.path.join(data_dir, request.param[0]))
    expected_hash = request.param[1]
    return input_vcf, expected_hash, request


class TestProcessVcf:
    def test_process_vcf(self, vcf_case, output_dir):
        input_vcf, expected_hash, request = vcf_case
        update_hashes = request.config.getoption("--update-vcf-hashes")

      
        str_bed = os.path.abspath(os.path.join(base_dir, "data", "GRCh38_repeats.bed"))
        output_filename = os.path.basename(input_vcf).replace(".vcf.gz", "").replace(".vcf", "") + ".processed.vcf"
        output_path = os.path.abspath(os.path.join(base_dir, "output", output_filename))
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        str_df = strvcf_annotator.load_str_reference(str_bed)
        strvcf_annotator.process_vcf(input_vcf, str_df, output_path)

        actual_hash = file_hash(output_path)

        if update_hashes:
            print(f"New hash for {input_vcf}: {actual_hash}")
            # optionally save to a JSON/yaml file here
        else:
            assert actual_hash == expected_hash, f"{input_vcf} hash mismatch: {actual_hash}"
