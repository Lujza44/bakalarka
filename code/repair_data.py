import json

wrong_strand = ['D1S1656', 'D2S1338', 'FGA', 'D5S818', 'CSF1PO', 'D7S820', 'vWA', 'PentaE', 'D19S433'] # done
wrong_after_flank = ['D1S1656', 'D5S818', 'D7S820', 'vWA', 'D13S317', 'D18S51', 'D21S11'] # done
wrong_before_flank = 'PentaD'
wrong_both = 'D19S433'
wrong_ref_allele = 'PentaE'

def complement_dna(s):
    trans_table = str.maketrans('ATCG', 'TAGC')
    return s.translate(trans_table)

def reverse_string(s):
    return s[::-1]

json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)

# unify strand
for marker_name in wrong_strand:
    marker = data["markers"].get(marker_name)

    ref_allele = marker['referenceAllele']
 
    seq = ref_allele['sequence']
    bef = ref_allele['before']
    aft = ref_allele['after']

    ref_allele['sequence'] = reverse_string(complement_dna(seq))
    ref_allele['before'] = reverse_string(complement_dna(aft))
    ref_allele['after'] = reverse_string(complement_dna(bef))

# repair wrong after flanking regions
for marker_name in wrong_after_flank:
    marker = data["markers"].get(marker_name)

    STRsize = marker["referenceAllele"]["STRsize"]
        
    for variant in marker["alleleVariants"]:
        allele_str = str(variant["allele"])
        if '.' in allele_str:
            integer_part, decimal_part = allele_str.split('.')
            integer_part = int(integer_part)
            decimal_part = int(decimal_part[0])  # Assuming only one digit after decimal
        else:
            integer_part = int(allele_str)
            decimal_part = 0

        if integer_part > 0:
            theSize = (integer_part * STRsize) + decimal_part
            
            for sequenceVariant in variant["sequenceVariants"]:
                sequence = sequenceVariant["sequence"]
                if len(sequence) > theSize:
                    remaining_sequence = sequence[theSize:]
                    sequenceVariant["sequence"] = sequence[:int(theSize)]
                        
                    for flankingVariant in sequenceVariant["flankingRegionsVariants"]:
                        flankingVariant["after"] = remaining_sequence + flankingVariant["after"]

# repair wrong before flanking regions (PentaD only)
marker = data["markers"].get(wrong_before_flank)
for allele_variant in marker['alleleVariants']:
    if allele_variant['allele'] > 0.0:
        for seq_variant in allele_variant['sequenceVariants']:
            first_5_letters = seq_variant['sequence'][:5]
            seq_variant['sequence'] = seq_variant['sequence'][5:]
            for flank_var in seq_variant['flankingRegionsVariants']:
                flank_var['before'] += first_5_letters


# repair wrong both flanking regions (D19S433 only)
marker = data["markers"].get(wrong_both)
for allele_variant in marker['alleleVariants']:
    if allele_variant['allele'] > 0.0:
        for seq_variant in allele_variant['sequenceVariants']:
            last_18_letters = seq_variant['sequence'][-18:]
            seq_variant['sequence'] = seq_variant['sequence'][:-18]    
            for flank_var in seq_variant['flankingRegionsVariants']:
                flank_var['after'] = last_18_letters + flank_var['after']

with open(json_file_path, 'w') as file:
    json.dump(data, file, indent=4)