import json

wrong_strand = ['D1S1656', 'D2S1338', 'FGA', 'D5S818', 'CSF1PO', 'D7S820', 'vWA', 'PentaE', 'D19S433']
wrong_after_flank = ['D1S1656', 'D5S818', 'D7S820', 'vWA', 'D13S317', 'D18S51', 'D21S11']
wrong_before_flank = 'PentaD'
wrong_both = 'D19S433'
wrong_ref = 'PentaE'

def complement_dna(s):
    trans_table = str.maketrans('ATCG', 'TAGC')
    return s.translate(trans_table)

def reverse_string(s):
    return s[::-1]

def unify_strand():
    for marker_name in wrong_strand:
        marker = data["markers"].get(marker_name)

        ref_allele = marker['referenceAllele']
    
        seq = ref_allele['sequence']
        bef = ref_allele['before']
        aft = ref_allele['after']

        ref_allele['sequence'] = reverse_string(complement_dna(seq))
        ref_allele['before'] = reverse_string(complement_dna(aft))
        ref_allele['after'] = reverse_string(complement_dna(bef))

def repair_after_flanks():
    for marker_name in wrong_after_flank:
        marker = data["markers"].get(marker_name)

        STRsize = marker["STRsize"]
            
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

def repair_before_flanks():
    marker = data["markers"].get(wrong_before_flank)
    for allele_variant in marker['alleleVariants']:
        if allele_variant['allele'] > 0.0:
            for seq_variant in allele_variant['sequenceVariants']:
                first_5_letters = seq_variant['sequence'][:5]
                seq_variant['sequence'] = seq_variant['sequence'][5:]
                for flank_var in seq_variant['flankingRegionsVariants']:
                    flank_var['before'] += first_5_letters

def repair_both_flanks():
    marker = data["markers"].get(wrong_both)
    for allele_variant in marker['alleleVariants']:
        if allele_variant['allele'] > 0.0:
            for seq_variant in allele_variant['sequenceVariants']:
                last_18_letters = seq_variant['sequence'][-18:]
                seq_variant['sequence'] = seq_variant['sequence'][:-18]    
                for flank_var in seq_variant['flankingRegionsVariants']:
                    flank_var['after'] = last_18_letters + flank_var['after']

def repair_ref_allele():
    marker = data["markers"].get(wrong_ref)
    ref_allele = marker.get('referenceAllele', {})
    after_sequence = ref_allele.get('after', '')

    sequence_to_find = "AAAGA" * 5  # This represents AAAGA repeated 5 times

    index = after_sequence.find(sequence_to_find)

    before_part = after_sequence[:index]  # Everything before the sequence
    after_part = after_sequence[index + len(sequence_to_find):]  # Everything after the sequence

    ref_allele['sequence'] = sequence_to_find
    ref_allele['before'] = before_part
    ref_allele['after'] = after_part
    ref_allele['chromosome'] = "Chr15"
    ref_allele['STRsize'] = 5

def find_repeats(sequence, repeat_length):
    if repeat_length == 0: return []
    repeats = set()
    for i in range(0, len(sequence), repeat_length):
        substring = sequence[i:i+repeat_length]
        if substring not in repeats:
            repeats.add(substring)
    return list(repeats)

def adjust_ref_allele():
    for marker, details in data['markers'].items():
        reference_allele = details.get('referenceAllele', {})
        str_size = details.get('STRsize', 4)  # Default to 4 if not specified
        ref_seq = reference_allele.get('sequence', '')

        repeats = find_repeats(ref_seq, str_size)
        details['repeats'] = repeats

        ref_before = reference_allele.get('before', '')
        ref_after = reference_allele.get('after', '')

        allele_variants = details.get('alleleVariants', [])
        for variant in allele_variants:
            for sequence_variant in variant.get('sequenceVariants', []):
                flanking_variants = sequence_variant.get('flankingRegionsVariants', [])
                if flanking_variants:
                    flanking_variant = flanking_variants[0]
                    bef_size = len(flanking_variant.get('before', ''))
                    aft_size = len(flanking_variant.get('after', ''))
                    reference_allele['before'] = ref_before[-bef_size:] if bef_size > 0 else ""
                    reference_allele['after'] = ref_after[:aft_size] if aft_size > 0 else ""
                    break

json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)

unify_strand()
repair_after_flanks()
repair_before_flanks() # PentaD only
repair_both_flanks() # D19S433 only # TODO
repair_ref_allele() # PentaE only
adjust_ref_allele()

d21 = data['markers'].get("D21S11")
d21['repeats'] = ["TCTA", "TCTG"]

with open(json_file_path, 'w') as file:
    json.dump(data, file, indent=4)