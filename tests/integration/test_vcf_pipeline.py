import hashlib
import os
import pytest
from strvcf_annotator import STRAnnotator

# Get base directory for test data
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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
            "test.vcf.gz",
            "0cb0cac184724cf9ea3ea27c037704cc"
        ),
        (
            "pindel_header.vcf",
            "9bd195a201d6b3317645ce5d44d40a2e"
        ),
        (
            "mutec2_indel.vcf.gz",
            "9ce2671c373aead335e0f2eba02bf537"
        ),
        (
            "TCGA-DC-6682.vcf",
            "0321d7668390582f4a3ee0f538c60b0c"
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

        # Use new API
        annotator = STRAnnotator(str_bed)
        annotator.annotate_vcf_file(input_vcf, output_path)

        actual_hash = file_hash(output_path)

        if update_hashes:
            print(f"New hash for {input_vcf}: {actual_hash}")
            # optionally save to a JSON/yaml file here
        else:
            assert actual_hash == expected_hash, f"{input_vcf} hash mismatch: {actual_hash}"
