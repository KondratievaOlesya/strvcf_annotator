"""Console script for strvcf_annotator."""

import argparse
from strvcf_annotator import full_pipeline

def main():
    parser = argparse.ArgumentParser(description="Annotate STR regions in VCFs using a BED file")
    parser.add_argument("--input_dir", required=True, help="Directory with input VCF files")
    parser.add_argument("--str_bed", required=True, help="Path to BED file with STR regions")
    parser.add_argument("--output_dir", required=True, help="Directory to store annotated VCFs")
    parser.add_argument("--temp_dir", required=True, help="Directory to store intermediate filtered VCFs")

    args = parser.parse_args()
    
    full_pipeline(
        input_dir=args.input_dir,
        str_bed_path=args.str_bed,
        output_dir=args.output_dir,
        temp_dir=args.temp_dir
    )

if __name__ == "__main__":
    app()
