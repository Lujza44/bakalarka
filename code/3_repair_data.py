import json

other_strand = ['D1S1656', 'D2S1338', 'FGA', 'D5S818', 'CSF1PO', 'D7S820', 'vWA', 'PentaE', 'D19S433']
wrong_after_flank = ['D13S317', 'D18S51', 'D19S433', 'D21S11']
wrong_before_flank = ['D1S1656', 'D5S818', 'D7S820', 'PentaD', 'vWA']
wrong_both = 'D19S433'
wrong_ref = 'PentaE'  
long_before_flank = ['PentaD', 'D22S1045']

def complement_dna(s):
    trans_table = str.maketrans('ATCG', 'TAGC')
    return s.translate(trans_table)

def reverse_string(s):
    return s[::-1]

def unify_strand():
    for marker_name in other_strand:
        marker = data["markers"].get(marker_name)
        for allele_var in marker['lengthVariants']:
            for seq_var in allele_var['sequenceVariants']:
                # reverse sequence
                seq = seq_var['sequence']
                seq_var['sequence'] = reverse_string(complement_dna(seq))
                if seq_var['flankingRegionsVariants']:
                    for flank_var in seq_var['flankingRegionsVariants']:
                        # reverse before
                        # reverse after
                        bef = flank_var['before']
                        aft = flank_var['after']
                        flank_var['before'] = reverse_string(complement_dna(aft))
                        flank_var['after'] = reverse_string(complement_dna(bef))

def repair_after_flanks():
    for marker_name in wrong_after_flank:
        marker = data["markers"].get(marker_name)

        STRsize = marker["STRlength"]
            
        for variant in marker["lengthVariants"]:
            allele_str = str(variant["numberOfRepeats"])
            if '.' in allele_str:
                integer_part, decimal_part = allele_str.split('.')
                integer_part = int(integer_part)
                decimal_part = int(decimal_part[0])
            else:
                integer_part = int(allele_str)
                decimal_part = 0

            if integer_part > 0:
                theSize = (integer_part * STRsize) + decimal_part
                
                for sequenceVariant in variant["sequenceVariants"]:
                    sequence = sequenceVariant["sequence"]
                    if len(sequence) > theSize:
                        remaining_sequence = sequence[theSize:] # skipnem theSize pismen
                        sequenceVariant["sequence"] = sequence[:int(theSize)] # v sequence nechavam prvych theSize pismen
                            
                        for flankingVariant in sequenceVariant["flankingRegionsVariants"]:
                            flankingVariant["after"] = remaining_sequence + flankingVariant["after"]

def repair_before_flanks():
    for marker_name in wrong_before_flank:
        marker = data["markers"].get(marker_name)

        STRsize = marker["STRlength"]
            
        for variant in marker["lengthVariants"]:
            allele_str = str(variant["numberOfRepeats"])
            if '.' in allele_str:
                integer_part, decimal_part = allele_str.split('.')
                integer_part = int(integer_part)
                decimal_part = int(decimal_part[0])
            else:
                integer_part = int(allele_str)
                decimal_part = 0

            if integer_part > 0:
                seq_size = (integer_part * STRsize) + decimal_part
                
                for sequenceVariant in variant["sequenceVariants"]:
                    sequence = sequenceVariant["sequence"]
                    if len(sequence) > seq_size:
                        to_move = len(sequence) - seq_size
                        to_move_sequence = sequence[:to_move] # vezmem prvych theSize pismen, chcem ich presunut
                        sequenceVariant["sequence"] = sequence[to_move:]
                            
                        for flankingVariant in sequenceVariant["flankingRegionsVariants"]:
                            flankingVariant["before"] += to_move_sequence

def repair_ref_allele():
    marker = data["markers"].get(wrong_ref)
    ref_allele = marker.get('referenceAllele', {})
    whole_seq = ref_allele.get('before', '')

    marker['chromosome'] = 15
    marker['STRlength'] = 5
    marker['startCoordinate'] = 96831015

    sequence_to_find = "TCTTT" * 5

    index = whole_seq.find(sequence_to_find)

    before_part = whole_seq[:index]  
    after_part = whole_seq[index + len(sequence_to_find):]

    ref_allele['sequence'] = sequence_to_find
    ref_allele['before'] = before_part
    ref_allele['after'] = after_part

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
        str_size = details.get('STRlength', 4)  # defaultne 4
        ref_seq = reference_allele.get('sequence', '')

        repeats = find_repeats(ref_seq, str_size)
        details['repeats'] = repeats

        ref_before = reference_allele.get('before', '')
        ref_after = reference_allele.get('after', '')

        allele_variants = details.get('lengthVariants', [])
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

def shorten_before_flank():
    for marker_name in long_before_flank:
        marker = data["markers"].get(marker_name)
        for variant in marker["lengthVariants"]:                
            for sequenceVariant in variant["sequenceVariants"]:  
                for flankingVariant in sequenceVariant["flankingRegionsVariants"]:
                    if marker_name == 'D22S1045':
                        flankingVariant["before"] = flankingVariant["before"][6:] # skipnem prvych 6 pismen
                    elif marker_name == 'PentaD':
                        flankingVariant["before"] = flankingVariant["before"][5:] # skipnem prvych 5 pismen

json_file_path = 'data/transformed_data.json'

# DATA Z JSNU
with open(json_file_path, 'r') as file:
    data = json.load(file)

# zjednotenie strandov vsetkych alel. variant, ktore su z opacneho
unify_strand()

# oprava 3' flanking oblasti, ktore boli v primarnych datach nespravne zaradene k rep. oblasti
repair_after_flanks()

# oprava 5' flanking oblasti, ktore boli v primarnych datach nespravne zaradene k rep. oblasti
repair_before_flanks() 

# oprava ref. alely lokusu Penta E
repair_ref_allele()

# osekanie flanking oblasti vsetkych ref. alel, ponechanie iba casti, ktore boli citane aj v CZ databaze
adjust_ref_allele()

# skratenie before flanking regions u lokusov D22 a PentaD, pretoze v referencii take dlhe nie su (neda sa porovnavat)
shorten_before_flank()

# oprava nadbytocnych referencnych repeats lokusu D21
d21 = data['markers'].get("D21S11")
d21['repeats'] = ["TCTA", "TCTG"]

with open(json_file_path, 'w') as file:
    json.dump(data, file, indent=4)