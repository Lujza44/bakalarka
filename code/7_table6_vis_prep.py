import pandas as pd
import json
import math

def align(dna_sequence, substring_size, repeats):
    result = []
    i = 0
    counter = 0
    while i < len(dna_sequence):
        if i + substring_size > len(dna_sequence):
            result.append(list(dna_sequence[i:]))
            break
        current_substring = dna_sequence[i:i+substring_size]
        if current_substring in repeats and counter == 0: 
            result.extend(list(dna_sequence[i:i+substring_size]))
            i += substring_size
        else:
            if counter > 0 and current_substring in repeats:
                fill_size = substring_size - counter
                result.extend([''] * fill_size)
                counter = 0
                continue
            result.append(dna_sequence[i])
            counter = (counter + 1) % substring_size
            i += 1
    if counter != 0 and counter < substring_size:
        result.extend([''] * (substring_size - counter))

    return result


def get_rs_numbers(length_variants):
    # Initialize sets for storing unique indices and numbers
    unique_before_rs_numbers = []
    unique_after_rs_numbers = []

    # Assuming 'length_variants' is your JSON data loaded as a Python dictionary
    for item in length_variants:
        for sequence_variant in item.get("sequenceVariants", []):
            for flanking_variant in sequence_variant.get("flankingRegionsVariants", []):
                if "beforeRsNumbers" in flanking_variant:
                    for rs_tuple in flanking_variant["beforeRsNumbers"]:
                        if rs_tuple not in unique_before_rs_numbers:
                            unique_before_rs_numbers.append(rs_tuple)
                if "afterRsNumbers" in flanking_variant:
                    for rs_tuple in flanking_variant["afterRsNumbers"]:
                        if rs_tuple not in unique_after_rs_numbers:
                            unique_after_rs_numbers.append(rs_tuple)
                
    print("Unique Before RS Numbers:", unique_before_rs_numbers)
    print("Unique After RS Numbers:", unique_after_rs_numbers)




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
    start_coordinate = marker_info.get('startCoordinate', 0)

    ref_num_of_repeats = reference_allele.get('numberOfRepeats', 0)
    ref_sequence = list(reference_allele.get('sequence', ''))
    ref_before = list(reference_allele.get('before', ''))
    ref_after = list(reference_allele.get('after', ''))

    before_length = len(ref_before)
    seq_length = len(ref_sequence)
    after_length = len(ref_after)
    
    # rozsirenenie ref. sekvencie, ak treba
    max_number_of_repeats = 0 # na zistenie o kolko treba rozsirit ref. sekvenciu
    for allele_var in marker_info['lengthVariants']:
        if allele_var["numberOfRepeats"] > max_number_of_repeats:
            max_number_of_repeats = allele_var["numberOfRepeats"]
    fillers = (math.ceil(max_number_of_repeats) - ref_num_of_repeats) * str_length
    if max_number_of_repeats > ref_num_of_repeats: # ak je ref. sekv. prikratka, rozsirime ju = vlozime prazdne miesta
        ref_sequence = ref_sequence + [''] * fillers
    
    # cislovanie repeats v ref. sekvencii
    numbering = ['' for _ in range(before_length)]
    for num in range(1, ref_num_of_repeats + 1):
        numbering.append(num)
        for i in range(str_length - 1):
            numbering.append('')

    # GRCh38 suradnice jednotlivych nukleotidov v ref. sekvencii
    coordinates = ['', '', 'GRCh38 coordinates']
    for num in range(start_coordinate - before_length, start_coordinate + seq_length):
        coordinates.append(num)
    coordinates += [''] * fillers
    for num in range(start_coordinate + seq_length, start_coordinate + seq_length + after_length):
        coordinates.append(num)

    # cislovanie vzdialenosti od repeat region
    distance = ['', '', 'Distance from repeat region']
    for i in range(0, before_length):
        distance.append(before_length - i)
    distance += [''] * fillers + [''] * seq_length
    for i in range(1, after_length + 1):
        distance.append(i)

    # append uvodnych riadkov
    rows.append([marker, 'Allele', 'Chr' + str(chromosome)] + numbering)
    rows.append(['', '', 'Reference sequence'] + ref_before + ref_sequence + ref_after)
    rows.append(coordinates)
    rows.append(distance)

    print(marker)
    get_rs_numbers(marker_info['lengthVariants'])

    # iterovanie cez vsetky varianty alel ZORADENE podla 'allele'
    for allele_var in sorted(marker_info['lengthVariants'], key=lambda x: x['numberOfRepeats']): 
        allele = allele_var['numberOfRepeats']
        if allele == 0: continue

        for seq_var in allele_var['sequenceVariants']:
            sequence = seq_var['sequence']

            # TODO align sequence - skusit podla ref. sekv.
            if allele % 1 != 0:
                sequence = align(sequence, str_length, repeats)

            if seq_var['flankingRegionsVariants']: # ak existuju flanking region varianty, treba ich pripojit do riadku
                for flank_var in seq_var['flankingRegionsVariants']:
                    before = flank_var['before']
                    after = flank_var['after']

                    before = list(before)
                    sequence = list(sequence)
                    after = list(after)

                    if len(ref_sequence) > len(sequence): # posun after flanking oblasti, aby bol zarovnany s ref.
                        shift = len(ref_sequence) - len(sequence)
                        row = before + sequence + [''] * shift + after
                    else:
                        row = before + sequence + after
                    
                    rows.append(['', allele, ''] + row)

            else: # ak neexistuju, tak namiesto nich prazdne miesta
                row = ['' for _ in range(before_length)] + list(sequence)
                rows.append(['', allele, ''] + row)
    rows.append([]) # prazdny riadkok na oddelenie od nasledujuceho markera

df = pd.DataFrame(rows)

csv_file_path = 'data/raw_vis6.csv'

df.to_csv(csv_file_path, index=False)