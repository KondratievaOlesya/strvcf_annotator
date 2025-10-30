import pandas as pd
import pysam
import subprocess
import os
import argparse
import tempfile
from constrain_dumb import estimate_genotype_from_support

def run_bedtools_intersect(vcf_path, bed_path, output_path):
    # Step 1: Normalize with bcftools
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as norm_tmp:
        norm_vcf = norm_tmp.name

    # Path to your reference genome FASTA (required for correct REF/ALT normalization)
    reference_fa = "/home/pho/muto_intl_filtered_vcf/GRCh38_full_analysis_set_plus_decoy_hla.fa"  # <-- Replace this!

    norm_command = [
        'bcftools', 'norm',
        '-m', '-',           # split multiallelics
        '-f', reference_fa,  # normalize REF/ALT
        '-O', 'v',
        '-o', norm_vcf,
        vcf_path
    ]
    subprocess.run(norm_command, check=True)

    # Step 2: bedtools intersect
    bedtools_command = [
        'bedtools', 'intersect',
        '-header',
        '-a', norm_vcf,
        '-b', bed_path
    ]
    with open(output_path, 'w') as out:
        subprocess.run(bedtools_command, stdout=out, check=True)

    # Optional cleanup
    os.remove(norm_vcf)


def load_str_reference(str_path):
    df = pd.read_csv(str_path, sep='\t', header=None)
    df.columns = ['CHROM', 'START', 'END', 'PERIOD', 'RU']
    df['START'] = df['START'] + 1 # BED are 0-based, but VCFs are 1-based
    # END is not included for bedtools, but we are including, so no need for +1
    df['COUNT'] = (df['END'] - df['START'] + 1) / df['PERIOD']
    df.sort_values(by=['CHROM', 'START'], inplace=True)
    return df


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


def infer_pindel_genotype(sample_data, threshold=0.6):
    # Get support values, defaulting to 0 if missing
    pp = sample_data.get('PP', 0)  # ALT on positive strand
    np = sample_data.get('NP', 0)  # ALT on negative strand
    pr = sample_data.get('PR', 0)  # TOTAL on positive strand
    nr = sample_data.get('NR', 0)  # TOTAL on negative strand
    
    fc = sample_data.get('FC', 0)  # Fragment calls
    if fc < 30:
        # dont bother
        return (None, None)
    
    # Avoid division by zero
    pos_ratio = (pp / pr) if pr > 0 else 0
    neg_ratio = (np / nr) if nr > 0 else 0

    allele1 = 1 if pos_ratio > threshold else 0
    allele2 = 1 if neg_ratio > threshold else 0
    
    constrain_alleles = estimate_genotype_from_support(sample_data)    
    
    my_allels = tuple(sorted([allele1, allele2]))
    #if constrain_alleles != my_allels:
        # print(f'Constrain suggestion: {constrain_alleles}\nThreshold filter: {my_allels}')
    return constrain_alleles

def get_gt(sample_data):
    gt = sample_data.get('GT', (None, None))
    if (gt == (None, None)) or (gt == ('.', '.')):
        is_pindel = any(field in sample_data for field in ['PP', 'NP', 'PB', 'NB', 'PR', 'NR'])
        if is_pindel:
            gt = infer_pindel_genotype(sample_data)
    return gt

def make_modified_header(vcf_in):
    info_to_replace=['RU', 'PERIOD', 'REF', 'PERFECT']
    format_to_replace=['REPCN']
    new_header = pysam.VariantHeader()

    # Copy contigs
    for contig in vcf_in.header.contigs.values():
        new_header.contigs.add(contig.name, length=contig.length)

    # Copy all INFO fields except ones we want to overwrite
    for key, val in vcf_in.header.info.items():
        if key not in info_to_replace:
            new_header.info.add(
                key, val.number, val.type, val.description
            )

    # Copy all FORMAT fields except ones we want to overwrite
    for key, val in vcf_in.header.formats.items():
        if key not in format_to_replace:
            new_header.formats.add(
                key, val.number, val.type, val.description
            )

    # Copy filters
    for key, val in vcf_in.header.filters.items():
        if key == "PASS":
            continue  # PASS is default, skip re-adding
        new_header.filters.add(key, val.number, val.type, val.description)

    # Copy samples
    for sample in vcf_in.header.samples:
        new_header.add_sample(sample)

    # Add custom fields (overwriting ones we excluded)
    new_header.info.add('RU', 1, 'String', 'Repeat unit')
    new_header.info.add('PERIOD', 1, 'Integer', 'Repeat period (length of unit)')
    new_header.info.add('REF', 1, 'Integer', 'Reference copy number')
    new_header.info.add('PERFECT', 1, 'String', 'Indicates if the repeat sequence is a perfect RU repeat (both REF and ALT)')
    new_header.formats.add('REPCN', 2, 'Integer', 'Genotype given in number of copies of the repeat motif')

    return new_header

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
    
    for record in vcf_in:
        while str_idx < len(str_list) and (
            str_list[str_idx]['CHROM'] != record.chrom or
            (str_list[str_idx]['CHROM'] == record.chrom and str_list[str_idx]['END'] < record.pos)
        ):
            str_idx += 1
        if str_idx >= len(str_list):
            break

        str_row = str_list[str_idx]
        
        variant_start = record.pos
        variant_end = record.pos + len(record.ref)

        if str_row['CHROM'] != record.chrom or variant_end < str_row['START'] or variant_end > str_row['END']:
            continue  # No overlap
        
        repeat_start = str_row['START']
            
        repeat_seq = extract_repeat_sequence(str_row)
        ref_base = record.ref
        alt_base = record.alts[0]
        
        # It should be the same as repeat_seq
        reference_seq = repeat_seq # apply_variant_to_repeat(record.pos, ref_base, ref_base, repeat_start, repeat_seq)
        tmp_seq = apply_variant_to_repeat(record.pos, ref_base, ref_base, repeat_start, repeat_seq)
        if tmp_seq != reference_seq:
            print(f"[WARNING] reference in vcf ({record.contig}:{record.pos} {record.alleles[0]}>{record.alleles[1]}) is not the same as in STR panel ({str_row['CHROM']}:{str_row['START']}-{str_row['END']} with RU={str_row['RU']})")
            print(f"Assumed STR: {repeat_seq}")
            print(f"Assumed reference: {tmp_seq}")
            print('Assume STR panel is right\n')
            
        mutated_seq = apply_variant_to_repeat(record.pos, ref_base, alt_base, repeat_start, repeat_seq)

        ru = str_row['RU']
        ref_len = count_repeat_units(reference_seq, ru)
        alt_len = count_repeat_units(mutated_seq, ru)
        perfect = (reference_seq == ru * ref_len) and (mutated_seq == ru * alt_len)    
        info = normalize_info_fields(record, vcf_in.header)

        info.update(
            {
                'RU': ru,
                'PERIOD': int(str_row['PERIOD']),
                'REF': int(ref_len),
                'PERFECT':  "TRUE" if perfect else "FALSE"
            }
        )
        filter_keys = record.filter.keys()
        new_record = vcf_out.new_record(
            contig=record.contig,
            start=repeat_start - 1, # It automatically adds 1 
            stop=str_row['END'],
            id='.',
            alleles=(reference_seq, mutated_seq),
            info=info,
            filter=filter_keys,
        )
        
        # Compare GTs of both samples
        sample_names = list(record.samples.keys())
        if len(sample_names) == 2:
            gt_0 = get_gt(record.samples[sample_names[0]])
            gt_1 = get_gt(record.samples[sample_names[1]])
            if any(gt is None or len(gt) > 2 for gt in [gt_0, gt_1]):
                gt_invalid_count += 1
                continue
            if gt_0 == gt_1:
                gt_equal_count += 1
                continue
            
    
        # Copy original FORMAT fields and add REPCN
        for sample in record.samples:
            new_record.samples[sample].update(record.samples[sample])
            gt = get_gt(record.samples[sample])
                                
            alleles = []
            for allele in gt:
                if allele == 0:
                    alleles.append(ref_len)
                elif allele == 1:
                    alleles.append(alt_len)
                else:
                    # No info
                    alleles.append(0)
            new_record.samples[sample]['GT'] = gt
            new_record.samples[sample]['REPCN'] = alleles

        vcf_out.write(new_record)
        written_count += 1
        if perfect:
            perfect_count += 1    

    vcf_out.close()
    
    total = written_count + gt_equal_count + gt_invalid_count

    print("========== Summary ==========")
    print(f"Total evaluated records: {total}")
    print(f"Written records:          {written_count} ({written_count / total:.1%})")
    print(f"PERFECT STRs:             {perfect_count} ({perfect_count / total:.1%})")
    print(f"Skipped (GT equal):       {gt_equal_count} ({gt_equal_count / total:.1%})")
    print(f"Skipped (GT invalid):     {gt_invalid_count} ({gt_invalid_count / total:.1%})")
    print("=============================")



def full_pipeline(input_dir, str_bed_path, output_dir, temp_dir):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    str_df = load_str_reference(str_bed_path)

    for filename in os.listdir(input_dir):
        if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
            input_path = os.path.join(input_dir, filename)
            base = os.path.splitext(os.path.splitext(filename)[0])[0]
            temp_filtered_vcf = os.path.join(temp_dir, f"{base}.filtered.vcf")
            output_vcf = os.path.join(output_dir, f"{base}.annotated.vcf")

            if os.path.exists(output_vcf):
                print(f"Skipping {filename} — already processed.")
                continue

            print(f"Processing {filename}...")
            run_bedtools_intersect(input_path, str_bed_path, temp_filtered_vcf)
            process_vcf(temp_filtered_vcf, str_df, output_vcf)
            print(f" → Output: {output_vcf}")


# full_pipeline(
#     input_dir="/media/pho/Elements1/muto_intl_vcf",
#     str_bed_path="/media/pho/Elements1/resources/h_sapiens/GRCh38_repeats.bed",
#     output_dir="/home/pho/muto_intl_str",
#     temp_dir="/home/pho/muto_intl_filtered_vcf"
# )

if __name__ == "__main__":
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


