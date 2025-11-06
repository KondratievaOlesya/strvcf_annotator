# Usage Examples

This page contains practical examples of how to use `strvcf_annotator` for various scenarios.

## Contents

1. [Basic usage](#basic-usage)
2. [Streaming](#streaming)
3. [Batch processing](#batch-processing)
4. [Filtering and analysis](#filtering-and-analysis)
5. [Custom parsers](#custom-parsers)
6. [Pipeline integration](#pipeline-integration)

## Basic usage

### Example 1: Simple annotation of a single file

```python
from strvcf_annotator import annotate_vcf

# The simplest way
annotate_vcf('input.vcf', 'repeats.bed', 'output.vcf')
```

### Example 2: Using the STRAnnotator class

```python
from strvcf_annotator import STRAnnotator

# Create annotator
annotator = STRAnnotator('repeats.bed')

# Annotate file
annotator.annotate_vcf_file('input.vcf', 'output.vcf')

# Get statistics
stats = annotator.get_statistics()
print(f"Processed {stats['total_regions']} STR regions")
```

### Example 3: CLI usage

```bash
# Annotate a single file
strvcf-annotator --input input.vcf --str-bed repeats.bed --output output.vcf

# With verbose logging
strvcf-annotator --input input.vcf --str-bed repeats.bed --output output.vcf --verbose
```

## Streaming

### Example 4: Processing large files

```python
import pysam
from strvcf_annotator import STRAnnotator

# Create annotator
annotator = STRAnnotator('repeats.bed')

# Open input and output files
vcf_in = pysam.VariantFile('large_input.vcf')
vcf_out = pysam.VariantFile('large_output.vcf', 'w', header=vcf_in.header)

# Streaming processing
count = 0
for record in annotator.annotate_vcf_stream(vcf_in):
    vcf_out.write(record)
    count += 1
    
    # Progress every 1000 records
    if count % 1000 == 0:
        print(f"Processed {count} records...")

vcf_out.close()
vcf_in.close()

print(f"Total processed records: {count}")
```

### Example 5: Streaming with filtering

```python
import pysam
from strvcf_annotator import STRAnnotator

annotator = STRAnnotator('repeats.bed')
vcf_in = pysam.VariantFile('input.vcf')

# Create a new header with STR fields
from strvcf_annotator.core.annotation import make_modified_header
new_header = make_modified_header(vcf_in)
vcf_out = pysam.VariantFile('filtered_output.vcf', 'w', header=new_header)

# Keep only perfect repeats
for record in annotator.annotate_vcf_stream(vcf_in):
    if record.info.get('PERFECT') == 'TRUE':
        vcf_out.write(record)

vcf_out.close()
vcf_in.close()
```

## Batch processing

### Example 6: Processing a directory

```python
from strvcf_annotator import STRAnnotator
import logging

# Configure logging to track progress
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Create annotator
annotator = STRAnnotator('repeats.bed')

# Batch processing
annotator.process_directory('vcf_files/', 'annotated_files/')
```

### Example 7: CLI batch processing

```bash
# Process all VCF files in a directory
strvcf-annotator --input-dir vcf_files/ --str-bed repeats.bed --output-dir annotated/
```

### Example 8: Skip already processed files

```python
from strvcf_annotator import STRAnnotator
from pathlib import Path

annotator = STRAnnotator('repeats.bed')

input_dir = Path('vcf_files')
output_dir = Path('annotated')
output_dir.mkdir(exist_ok=True)

for vcf_file in input_dir.glob('*.vcf'):
    output_file = output_dir / f"{vcf_file.stem}.annotated.vcf"
    
    # Skip already processed files
    if output_file.exists():
        print(f"Skipping {vcf_file.name} - already processed")
        continue
    
    print(f"Processing {vcf_file.name}...")
    annotator.annotate_vcf_file(str(vcf_file), str(output_file))
```

## Filtering and analysis

### Example 9: STR region statistics

```python
from strvcf_annotator import STRAnnotator

annotator = STRAnnotator('repeats.bed')

# Get statistics
stats = annotator.get_statistics()

print("=== STR region statistics ===")
print(f"Total regions: {stats['total_regions']}")
print(f"Chromosomes: {stats['chromosomes']}")
print(f"Unique repeat units: {stats['unique_repeat_units']}")
print(f"Mean repeat count: {stats['mean_repeat_count']:.2f}")
print(f"Median repeat count: {stats['median_repeat_count']:.2f}")

print("\nDistribution by period:")
for period, count in sorted(stats['period_distribution'].items()):
    print(f"  Period {period}: {count} regions")
```

### Example 10: Finding STRs at specific positions

```python
from strvcf_annotator import STRAnnotator

annotator = STRAnnotator('repeats.bed')

# Search for STRs at specific positions
positions = [
    ('chr1', 1000000),
    ('chr2', 2500000),
    ('chr3', 500000)
]

for chrom, pos in positions:
    str_region = annotator.get_str_at_position(chrom, pos)
    
    if str_region:
        print(f"{chrom}:{pos}")
        print(f"  Repeat unit: {str_region['RU']}")
        print(f"  Period: {str_region['PERIOD']}")
        print(f"  Count: {str_region['COUNT']:.1f}")
        print(f"  Region: {str_region['START']}-{str_region['END']}")
    else:
        print(f"{chrom}:{pos} - no STR found")
```

### Example 11: Filtering by repeat type

```python
import pysam
from strvcf_annotator import STRAnnotator

annotator = STRAnnotator('repeats.bed')
vcf_in = pysam.VariantFile('input.vcf')

# Create separate files for different repeat types
from strvcf_annotator.core.annotation import make_modified_header
new_header = make_modified_header(vcf_in)

trinucleotide_vcf = pysam.VariantFile('trinucleotide.vcf', 'w', header=new_header)
tetranucleotide_vcf = pysam.VariantFile('tetranucleotide.vcf', 'w', header=new_header)

for record in annotator.annotate_vcf_stream(vcf_in):
    period = record.info.get('PERIOD')
    
    if period == 3:
        trinucleotide_vcf.write(record)
    elif period == 4:
        tetranucleotide_vcf.write(record)

trinucleotide_vcf.close()
tetranucleotide_vcf.close()
vcf_in.close()
```

### Example 12: Exporting to CSV for downstream analysis

```python
import pysam
import csv
from strvcf_annotator import STRAnnotator

annotator = STRAnnotator('repeats.bed')
vcf_in = pysam.VariantFile('input.vcf')

# Export annotations to CSV
with open('str_annotations.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['CHROM', 'POS', 'RU', 'PERIOD', 'REF_COPIES', 'PERFECT', 'SAMPLE', 'REPCN'])
    
    for record in annotator.annotate_vcf_stream(vcf_in):
        chrom = record.contig
        pos = record.pos
        ru = record.info.get('RU')
        period = record.info.get('PERIOD')
        ref_copies = record.info.get('REF')
        perfect = record.info.get('PERFECT')
        
        for sample_name, sample_data in record.samples.items():
            repcn = sample_data.get('REPCN', (0, 0))
            writer.writerow([chrom, pos, ru, period, ref_copies, perfect, sample_name, f"{repcn[0]}/{repcn[1]}"])

vcf_in.close()
print("Export complete: str_annotations.csv")
```

## Custom parsers

### Example 13: Simple custom parser

```python
from strvcf_annotator import STRAnnotator
from strvcf_annotator.parsers.base import BaseVCFParser

class SimpleParser(BaseVCFParser):
    """Simple parser that returns the sample's genotype if present."""
    
    def get_genotype(self, record, sample_idx):
        samples = list(record.samples.values())
        if sample_idx >= len(samples):
            return None
        
        gt = samples[sample_idx].get('GT', None)
        if gt and gt != (None, None):
            return gt
        return None
    
    def has_variant(self, record, sample_idx):
        gt = self.get_genotype(record, sample_idx)
        return gt is not None and any(a != 0 for a in gt)
    
    def extract_info(self, record, sample_idx):
        samples = list(record.samples.values())
        if sample_idx >= len(samples):
            return {}
        
        info = {}
        sample = samples[sample_idx]
        
        if 'DP' in sample:
            info['DP'] = sample['DP']
        if 'AD' in sample:
            info['AD'] = sample['AD']
        
        return info
    
    def validate_record(self, record):
        return record is not None and len(record.samples) > 0

# Usage
parser = SimpleParser()
annotator = STRAnnotator('repeats.bed', parser=parser)
annotator.annotate_vcf_file('input.vcf', 'output.vcf')
```

### Example 14: Parser with quality filtering

```python
from strvcf_annotator.parsers.base import BaseVCFParser

class QualityFilterParser(BaseVCFParser):
    """Parser with quality-based filtering."""
    
    def __init__(self, min_gq=20, min_dp=10):
        self.min_gq = min_gq
        self.min_dp = min_dp
    
    def get_genotype(self, record, sample_idx):
        samples = list(record.samples.values())
        if sample_idx >= len(samples):
            return None
        
        sample = samples[sample_idx]
        
        # Quality checks
        gq = sample.get('GQ', 0)
        dp = sample.get('DP', 0)
        
        if gq < self.min_gq or dp < self.min_dp:
            return None
        
        return sample.get('GT', None)
    
    def has_variant(self, record, sample_idx):
        gt = self.get_genotype(record, sample_idx)
        return gt is not None and any(a != 0 for a in gt)
    
    def extract_info(self, record, sample_idx):
        samples = list(record.samples.values())
        if sample_idx >= len(samples):
            return {}
        return {'DP': samples[sample_idx].get('DP', 0)}
    
    def validate_record(self, record):
        return record is not None

# Usage with custom thresholds
parser = QualityFilterParser(min_gq=30, min_dp=15)
annotator = STRAnnotator('repeats.bed', parser=parser)
annotator.annotate_vcf_file('input.vcf', 'high_quality_output.vcf')
```

## Pipeline integration

### Example 15: Integration with pandas

```python
import pysam
import pandas as pd
from strvcf_annotator import STRAnnotator

def vcf_to_dataframe(vcf_path, str_bed_path):
    """Convert an annotated VCF to a pandas DataFrame."""
    
    annotator = STRAnnotator(str_bed_path)
    vcf_in = pysam.VariantFile(vcf_path)
    
    data = []
    for record in annotator.annotate_vcf_stream(vcf_in):
        for sample_name, sample_data in record.samples.items():
            row = {
                'chrom': record.contig,
                'pos': record.pos,
                'ru': record.info.get('RU'),
                'period': record.info.get('PERIOD'),
                'ref_copies': record.info.get('REF'),
                'perfect': record.info.get('PERFECT') == 'TRUE',
                'sample': sample_name,
                'gt': sample_data.get('GT'),
                'repcn': sample_data.get('REPCN')
            }
            data.append(row)
    
    vcf_in.close()
    return pd.DataFrame(data)

# Usage
df = vcf_to_dataframe('input.vcf', 'repeats.bed')
print(df.head())

# Analysis
print(f"\nTotal variants: {len(df)}")
print(f"Perfect repeats: {df['perfect'].sum()}")
print(f"\nDistribution by period:")
print(df['period'].value_counts())
```

### Example 16: Parallel processing

```python
from strvcf_annotator import STRAnnotator
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import logging

logging.basicConfig(level=logging.INFO)

def process_single_file(args):
    """Process a single file."""
    vcf_file, str_bed_path, output_dir = args
    
    try:
        annotator = STRAnnotator(str_bed_path)
        output_file = output_dir / f"{vcf_file.stem}.annotated.vcf"
        
        annotator.annotate_vcf_file(str(vcf_file), str(output_file))
        return f"Success: {vcf_file.name}"
    
    except Exception as e:
        return f"Error {vcf_file.name}: {e}"

def parallel_process_directory(input_dir, str_bed_path, output_dir, max_workers=4):
    """Parallel processing of a directory."""
    
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Prepare file list
    vcf_files = list(input_path.glob('*.vcf'))
    args_list = [(f, str_bed_path, output_path) for f in vcf_files]
    
    # Parallel processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(process_single_file, args_list)
    
    # Print results
    for result in results:
        print(result)

# Usage
parallel_process_directory('vcf_files/', 'repeats.bed', 'annotated/', max_workers=4)
```

### Example 17: Integration with Snakemake

```python
# Snakefile
rule annotate_str:
    input:
        vcf="data/{sample}.vcf",
        bed="resources/repeats.bed"
    output:
        vcf="results/{sample}.annotated.vcf"
    log:
        "logs/{sample}.log"
    script:
        "scripts/annotate_str.py"

# scripts/annotate_str.py
from strvcf_annotator import annotate_vcf
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO
)

annotate_vcf(
    snakemake.input.vcf,
    snakemake.input.bed,
    snakemake.output.vcf
)
```

### Example 18: Processing with logging

```python
from strvcf_annotator import STRAnnotator
import logging
from datetime import datetime

# Configure logging
log_file = f"annotation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

# Processing with logging
try:
    logger.info("Starting annotation")
    
    annotator = STRAnnotator('repeats.bed')
    stats = annotator.get_statistics()
    logger.info(f"Loaded {stats['total_regions']} STR regions")
    
    annotator.annotate_vcf_file('input.vcf', 'output.vcf')
    
    logger.info("Annotation completed successfully")
    
except Exception as e:
    logger.error(f"Error during annotation: {e}", exc_info=True)
    raise
```

## Additional examples

More examples can be found in the `examples/` directory of the repository:

* `basic_usage.py` - Basic usage examples
