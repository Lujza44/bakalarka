import pandas as pd
import json
import math

json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)

rows = [[]]

# iterovanie cez vsetky markery v JSNE ZORADENE podla 'chromosome'
for marker, marker_info in sorted(data['markers'].items(), key=lambda x: x[1].get('chromosome', 0)):
    if marker == 'D19S433': continue

    str_length = marker_info['STRlength']
    repeats = marker_info['repeats']
    if not str_length or not repeats:
        continue

    reference_allele = marker_info.get('referenceAllele', {})
    chromosome = marker_info.get('chromosome', 0)

    ref_num_of_repeats = reference_allele.get('numberOfRepeats', 0)
    ref_sequence = list(reference_allele.get('sequence', ''))
    ref_before = list(reference_allele.get('before', ''))
    ref_after = list(reference_allele.get('after', ''))

    before_length = len(ref_before)
    after_length = len(ref_after)
    
    max_number_of_repeats = 0 # na zistenie o kolko treba rozsirit ref. sekvenciu
    for allele_var in marker_info['lengthVariants']:
        if allele_var["numberOfRepeats"] > max_number_of_repeats:
            max_number_of_repeats = allele_var["numberOfRepeats"]
    
    if max_number_of_repeats > ref_num_of_repeats: # ak je ref. sekv. prikratka, vlozime prazdne miesta
        ref_sequence = ref_sequence + [''] * ((math.ceil(max_number_of_repeats) - ref_num_of_repeats) * str_length)

    numbering = ['' for _ in range(before_length)]
    for num in range(1, ref_num_of_repeats + 1):
        numbering.append(num)
        for i in range(str_length - 1):
            numbering.append('')
    rows.append([marker, 'Allele', 'Chr' + str(chromosome)] + numbering)
    rows.append(['', '', ''] + ref_before + ref_sequence + ref_after)



    # iterovanie cez vsetky varianty alel ZORADENE podla 'allele'
    for allele_var in sorted(marker_info['lengthVariants'], key=lambda x: x['numberOfRepeats']): 
        allele = allele_var['numberOfRepeats']
        for seq_var in allele_var['sequenceVariants']:
            sequence = seq_var['sequence']
            if seq_var['flankingRegionsVariants']: # ak existuju flanking region varianty
                for flank_var in seq_var['flankingRegionsVariants']:
                    before = flank_var['before']
                    after = flank_var['after']
                    before_indexes = flank_var["beforeSNPIndices"]
                    before_rs = flank_var["beforeRsNumbers"]
                    after_indexes = flank_var["afterSNPIndices"]
                    after_rs = flank_var["afterRsNumbers"]

                    before = list(before)
                    sequence = list(sequence)
                    after = list(after)

                    if len(ref_sequence) > len(sequence): # posun after flanking oblasti, aby bol zarovnany s ref.
                        shift = len(ref_sequence) - len(sequence)
                        row = before + sequence + [''] * shift + after
                    else:
                        row = before + sequence + after
                    
                    rows.append(['', '', ''] + row)

            else: # ak neexistuju
                row = ['' for _ in range(before_length)] + list(sequence)
                rows.append(['', '', ''] + row)
    rows.append([]) # prazdne riadky na oddelenie
    rows.append([])
    rows.append([])

df = pd.DataFrame(rows)

csv_file_path = 'data/raw_vis6.csv'

df.to_csv(csv_file_path, index=False)