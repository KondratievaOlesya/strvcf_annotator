"""Main module."""
import pandas as pd
import pysam
import subprocess
import os
import argparse
import tempfile
import pyranges as pr
# from constrain_dumb import estimate_genotype_from_support
import logging

logger = logging.getLogger(__name__)


def load_str_reference(str_path):
    df = pd.read_csv(str_path, sep='\t', header=None)
    df.columns = ['CHROM', 'START', 'END', 'PERIOD', 'RU']
    df['START'] = df['START'] + 1 # BED are 0-based, but VCFs are 1-based
    # END is not included for bedtools, but we are including, so no need for +1
    df['COUNT'] = (df['END'] - df['START'] + 1) / df['PERIOD']
    df.sort_values(by=['CHROM', 'START'], inplace=True)
    return df

def make_modified_header(vcf_in: pysam.VariantFile) -> pysam.VariantHeader:
    """
    Create a modified VCF header with updated INFO and FORMAT fields 
    for STR-specific annotations.

    This function takes an input VCF header from a `pysam.VariantFile` object,
    copies all contigs, filters, and samples, retains existing INFO and FORMAT 
    fields except a selected list, and adds new STR-related annotations 
    to the INFO and FORMAT sections. Existing fields with the same 
    names as the new ones are **skipped**.

    INFO fields replaced:
        - `RU`: Repeat unit
        - `PERIOD`: Repeat period (length of unit)
        - `REF`: Reference copy number
        - `PERFECT`: Indicates perfect repeats in REF and ALT

    FORMAT field replaced:
        - `REPCN`: Genotype as number of repeat motif copies

    Parameters
    ----------
    vcf_in : pysam.VariantFile
        The input VCF file object from which the header is derived.

    Returns
    -------
    pysam.VariantHeader
        A new VariantHeader object with updated fields suitable 
        for STR annotation.

    Notes
    -----
    The function does not modify the input header in-place. It returns a new 
    `VariantHeader` object with the modified annotations.
    """
    info_to_replace = ['RU', 'PERIOD', 'REF', 'PERFECT']
    format_to_replace = ['REPCN']
    new_header = pysam.VariantHeader()

    # Step 1: Copy raw header lines, except the ones we want to modify
    for rec in vcf_in.header.records:
        if rec.key == "INFO" and rec['ID'] in info_to_replace:
            continue
        if rec.key == "FORMAT" and rec['ID'] in format_to_replace:
            continue
        try:
            new_header.add_record(rec)
        except Exception:
            # fallback: skip unrecognized or malformed lines
            pass

    # Step 2: Copy contigs explicitly (optional if already in records)
    for contig in vcf_in.header.contigs.values():
        if contig.name not in new_header.contigs:
            new_header.contigs.add(contig.name, length=contig.length)

    # Step 3: Add samples
    for sample in vcf_in.header.samples:
        new_header.add_sample(sample)

    # Step 4: Add modified INFO and FORMAT fields
    new_header.info.add('RU', 1, 'String', 'Repeat unit')
    new_header.info.add('PERIOD', 1, 'Integer', 'Repeat period (length of unit)')
    new_header.info.add('REF', 1, 'Integer', 'Reference copy number')
    new_header.info.add('PERFECT', 1, 'String', 'Indicates if the repeat sequence is a perfect RU repeat (both REF and ALT)')
    new_header.formats.add('REPCN', 2, 'Integer', 'Genotype given in number of copies of the repeat motif')

    return new_header

def get_gt(sample_data):
    gt = sample_data.get('GT', (None, None))
    if (gt == (None, None)) or (gt == ('.', '.')):
        is_pindel = any(field in sample_data for field in ['PP', 'NP', 'PB', 'NB', 'PR', 'NR'])
        if is_pindel:
            gt = infer_pindel_genotype(sample_data)
    return gt


def extract_repeat_sequence(str_row):
    return str_row['RU'] * int(str_row['COUNT'])

def apply_variant_to_repeat(pos, ref, alt, repeat_start, repeat_seq):
    """
    Applies a variant (ref → alt) to the STR repeat_seq that starts at repeat_start.
    The variant starts at `pos`. Handles variants that start before the repeat.
    """
    relative_pos = pos - repeat_start

    # Variant starts before the STR
    if relative_pos < 0:
        # Clip how much of the REF applies before the STR
        ref_clip = -relative_pos
        ref = ref[ref_clip:]  # remove the part before STR
        alt = alt[ref_clip:]  # same for ALT
        relative_pos = 0      # mutation starts at beginning of repeat_seq

    # Apply mutation safely inside repeat_seq
    before_mut = repeat_seq[:relative_pos]
    after_mut = repeat_seq[relative_pos + len(ref):] if relative_pos + len(ref) <= len(repeat_seq) else ''
    
    mutated = before_mut + alt + after_mut
    return mutated

def count_repeat_units(seq, motif):
    """Count non-overlapping occurrences of `motif` anywhere in `seq`."""
    count = 0
    i = 0
    mlen = len(motif)
    while i <= len(seq) - mlen:
        if seq[i:i + mlen] == motif:
            count += 1
            i += mlen  # skip forward by motif length (non-overlapping)
        else:
            i += 1
    return count

def normalize_info_fields(record, header):
    fixed_info = {}

    for key, val in record.info.items():
        if key not in header.info:
            continue  # skip unknown

        desc = header.info[key]

        # Handle Flags (should be included as key only if present)
        if desc.type == "Flag":
            if val:
                fixed_info[key] = True  # will be serialized correctly
            continue  # no value assignment needed

        # Handle single String field (but given as tuple or list)
        if desc.type == "String" and desc.number == 1:
            if isinstance(val, (list, tuple)):
                fixed_info[key] = "|".join(map(str, val))
            else:
                fixed_info[key] = val
            continue
        
        # Clip R fields (REF + ALT)
        if desc.number == "R" and isinstance(val, (list, tuple)):
            fixed_info[key] = list(val[:2])
            continue

        # Convert other tuple-like values to list (e.g., for Type=R or A)
        # if isinstance(val, tuple):
        #     fixed_info[key] = list(val)
        # else:
        #     fixed_info[key] = val

    return fixed_info

def check_vcf_sorted(vcf_in):
    """Check if VCF is sorted by chrom and position."""
    last_chrom = None
    last_pos = -1
    for rec in vcf_in:
        chrom = rec.contig
        pos = rec.pos
        if last_chrom is not None:
            if chrom < last_chrom or (chrom == last_chrom and pos < last_pos):
                return False
        last_chrom, last_pos = chrom, pos
    return True


def reset_and_sort_vcf(vcf_in):
    """Rewind and sort VCF records in memory."""
    header = vcf_in.header
    records = list(vcf_in)
    contig_order = {c: i for i, c in enumerate(header.contigs.keys())}
    records.sort(key=lambda r: (contig_order.get(r.contig, float("inf")), r.pos))
    return records


def should_skip_genotype(record):
    """Determine if the record should be skipped based on sample genotypes."""
    samples = list(record.samples.keys())
    if len(samples) != 2:
        return False
    gt_0 = get_gt(record.samples[samples[0]])
    gt_1 = get_gt(record.samples[samples[1]])
    if any(gt is None or len(gt) > 2 for gt in [gt_0, gt_1]):
        return True
    return gt_0 == gt_1


def build_new_record(
    record: pysam.VariantRecord,
    str_row: pd.Series,
    header: pysam.VariantHeader
) -> pysam.VariantRecord:
    """
    Build a new VCF record with STR-based alleles and annotations.

    Given a mutation record (e.g., SNV or indel) and a repeat region (STR) 
    that overlaps the mutation, this function constructs a new VCF record 
    where the alleles represent the full repeat sequences (before and after 
    mutation). It also adds STR-specific annotations to INFO and FORMAT fields.

    The resulting record uses the entire STR repeat as the reference allele 
    and applies the variant to this repeat to generate the alternate allele 
    (e.g., "AAAAAA" > "A" for a deletion). The function then calculates the 
    repeat copy number for both alleles and annotates:
    
    - INFO fields: `RU`, `PERIOD`, `REF`, `PERFECT`
    - FORMAT field: `REPCN` (repeat copy number per allele per sample)

    Parameters
    ----------
    record : pysam.VariantRecord
        The original VCF record with a point mutation or indel overlapping the STR.

    str_row : dict or pandas.Series
        STR metadata corresponding to the region containing the mutation.
        Must include `'CHROM'`, `'START'`, `'END'`, `'RU'`, and `'PERIOD'`.

    header : pysam.VariantHeader
        A modified VCF header containing all required STR-specific INFO and FORMAT fields.

    Returns
    -------
    pysam.VariantRecord
        A new VCF record with alleles defined as full STR sequences and annotated
        with STR-related metadata in INFO and FORMAT fields.

    Notes
    -----
    - The function logs a warning if the reference base in the VCF does not match
      the STR sequence at the mutation position.
    - The reference allele is reconstructed from the STR region; the mutation is
      applied to it to get the alternate allele.
    - The `PERFECT` field is marked TRUE only if both REF and ALT consist of 
      exact multiples of the repeat unit.
    """
    
    repeat_start = str_row['START']
    repeat_seq = extract_repeat_sequence(str_row)
    ref_base = record.ref
    alt_base = record.alts[0]

    reference_seq = repeat_seq
    tmp_seq = apply_variant_to_repeat(record.pos, ref_base, ref_base, repeat_start, repeat_seq)
    if tmp_seq != reference_seq:
        logger.warning(
            f"Reference mismatch: VCF {record.contig}:{record.pos} {record.alleles[0]}>{record.alleles[1]}, "
            f"STR panel {str_row['CHROM']}:{str_row['START']}-{str_row['END']} RU={str_row['RU']}"
        )

    mutated_seq = apply_variant_to_repeat(record.pos, ref_base, alt_base, repeat_start, repeat_seq)

    ru = str_row['RU']
    ref_len = count_repeat_units(reference_seq, ru)
    alt_len = count_repeat_units(mutated_seq, ru)
    perfect = (reference_seq == ru * ref_len) and (mutated_seq == ru * alt_len)

    info = normalize_info_fields(record, header)
    info.update({
        'RU': ru,
        'PERIOD': int(str_row['PERIOD']),
        'REF': int(ref_len),
        'PERFECT': "TRUE" if perfect else "FALSE"
    })

    new_record = header.new_record(
        contig=record.contig,
        start=repeat_start - 1,
        stop=str_row['END'],
        id='.',
        alleles=(reference_seq, mutated_seq),
        info=info,
        filter=record.filter.keys()
    )

    for sample in record.samples:
        new_record.samples[sample].update(record.samples[sample])
        gt = get_gt(record.samples[sample])
        alleles = [
            ref_len if allele == 0 else
            alt_len if allele == 1 else
            0 for allele in gt
        ]
        new_record.samples[sample]['GT'] = gt
        new_record.samples[sample]['REPCN'] = alleles

    return new_record


def generate_annotated_records(vcf_in, str_df):
    """
    Generator that yields annotated VCF records.
    
    Args:
        vcf_in: pysam.VariantFile input VCF (already opened)
        str_df: pandas DataFrame with STR BED info (CHROM, START, END, PERIOD, RU)
    
    Yields:
        Annotated pysam.VariantRecord
    """

    header = make_modified_header(vcf_in)

    # Warn if VCF is unsorted
    if not check_vcf_sorted(vcf_in):
        logger.warning("Input VCF is not sorted - sorting in memory.")
        records = reset_and_sort_vcf(vcf_in)
    else:
        vcf_in.reset()
        records = vcf_in.fetch()

    str_idx = 0
    str_list = str_df.sort_values(by=["CHROM", "START"]).to_dict("records")

    for record in vcf_in.fetch():
        while str_idx < len(str_list) and (
            str_list[str_idx]['CHROM'] != record.chrom or
            (str_list[str_idx]['CHROM'] == record.chrom and str_list[str_idx]['END'] < record.pos)
        ):
            str_idx += 1
        if str_idx >= len(str_list):
            break

        str_row = str_list[str_idx]

        variant_end = record.pos + len(record.ref)
        if str_row['CHROM'] != record.chrom or variant_end < str_row['START'] or variant_end > str_row['END']:
            continue  # no overlap

        if should_skip_genotype(record):
            continue

        yield build_new_record(record, str_row, header)



def annotate_vcf(vcf_in, str_df):
    from tempfile import SpooledTemporaryFile
    header = make_modified_header(vcf_in)
    tmp = SpooledTemporaryFile(mode='w+')
    vcf_out = pysam.VariantFile(tmp, 'w', header=header)

    for rec in generate_annotated_records(vcf_in, str_df):
        vcf_out.write(rec)

    vcf_out.close()
    tmp.seek(0)
    return pysam.VariantFile(tmp)

# Write to file
def process_vcf(vcf_path, str_df, output_path):
    vcf_in = pysam.VariantFile(vcf_path)
    new_header = make_modified_header(vcf_in)
    vcf_out = pysam.VariantFile(output_path, 'w', header=new_header)


    str_idx = 0
    str_list = str_df.to_dict('records')

    # Initialize stats
    perfect_count = 0
    gt_equal_count = 0
    gt_invalid_count = 0
    written_count = 0
    for record in generate_annotated_records(vcf_in, str_df):
        vcf_out.write(record)
        
    vcf_out.close()
    
    total = written_count + gt_equal_count + gt_invalid_count

    print("========== Summary ==========")
    print(f"Total evaluated records: {total}")
    print(f"Written records:          {written_count} ({written_count / total:.1%})")
    print(f"PERFECT STRs:             {perfect_count} ({perfect_count / total:.1%})")
    print(f"Skipped (GT equal):       {gt_equal_count} ({gt_equal_count / total:.1%})")
    print(f"Skipped (GT invalid):     {gt_invalid_count} ({gt_invalid_count / total:.1%})")
    print("=============================")



def full_pipeline(input_dir, str_bed_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    str_df = load_str_reference(str_bed_path)

    for filename in os.listdir(input_dir):
        if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
            input_path = os.path.join(input_dir, filename)
            base = os.path.splitext(os.path.splitext(filename)[0])[0]
            output_vcf = os.path.join(output_dir, f"{base}.annotated.vcf")

            if os.path.exists(output_vcf):
                print(f"Skipping {filename} — already processed.")
                continue

            print(f"Processing {filename}...")
            process_vcf(input_path, str_df, output_vcf)
            print(f" → Output: {output_vcf}")



full_pipeline(
    input_dir="/media/pho/Elements1/muto_intl_vcf",
    str_bed_path="/media/pho/Elements1/resources/h_sapiens/GRCh38_repeats.bed",
    output_dir="/home/pho/muto_intl_str"
)