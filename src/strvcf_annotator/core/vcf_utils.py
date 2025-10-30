# VCF parsing, normalization, genotype handling
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
   