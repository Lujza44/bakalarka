import pandas as pd
import json
import re


def to_bracket_notation(dna_sequence, substring_size, repeats):
    result = ""
    i = 0
    count = 1
    last_space = -1
    while i < len(dna_sequence):
        if i + substring_size > len(dna_sequence):
            result += dna_sequence[i:] # pripojenie zvysku sekvencie, ak je < STRsize
            break
        current_substring = dna_sequence[i:i+substring_size]
        if current_substring in repeats: # ak na tejto pozicii je repeat, pripojim do vysledku a presuniem sa na dalsi
            if i + substring_size * 2 <= len(dna_sequence) and dna_sequence[i:i+substring_size] == dna_sequence[i+substring_size:i+substring_size*2]:
                count += 1
            else:
                if count > 1:
                    result += "[" + current_substring + "]" + str(count) + " "
                else:
                    result += current_substring + " "
                count = 1 
                last_space = len(result) - 1
            i += substring_size
        else:
            if len(result) - last_space >= substring_size + 1:
                result += " "
                last_space = len(result) - 1
            result += dna_sequence[i] # pripojim len pismenko ak na tejto pozicii nie je repeat
            next_substring = dna_sequence[i+1:i+1+substring_size]
            if next_substring in repeats:
                result += " "
                last_space = len(result) - 1
            i += 1
    return result.strip()

def chromosome_key(chromosome_str):
    match = re.search(r'\d+', chromosome_str)
    if match:
        return int(match.group())
    return 0  # Default value for strings without numbers

#json_file_path = 'data/transformed_data.json'
json_file_path = 'data/repaired.json'
with open(json_file_path, 'r') as file:
    data = json.load(file)

rows = [[]]

# iterovanie cez vsetky markery ZORADENE podla 'chromosome'
for marker, marker_info in sorted(data['markers'].items(), key=lambda x: chromosome_key(x[1].get('chromosome', '0'))):
    sum_count = 0
    sum_freq = 0

    str_length = marker_info['STRsize']
    repeats = marker_info['repeats']

    if not str_length or not repeats:
        continue

    # iterovanie cez vsetky varianty alel ZORADENE podla 'allele'
    for allele_var in sorted(marker_info['alleleVariants'], key=lambda x: x['allele']): 
        allele = allele_var['allele']
        for seq_var in allele_var['sequenceVariants']:
            sequence = seq_var['sequence']

            sequence = to_bracket_notation(seq_var['sequence'], str_length, repeats)
            
            if seq_var['flankingRegionsVariants']: # ak existuju flanking region varianty, kazdy bude mat vlastny riadok
                for flank_var in seq_var['flankingRegionsVariants']:
                    before = flank_var['before']
                    after = flank_var['after']
                    count = flank_var['count']
                    frequency = flank_var['frequency']
                    sum_count += count
                    sum_freq += frequency
                    rows.append([marker, allele, sequence, "", count, frequency, before, after])
            else: # ak neexistuju, bunky ostanu prazdne
                rows.append([marker, allele, sequence, "", "", "", "", ""])
    rows.append(["", "", "", "", sum_count, sum_freq, "", ""])
    rows.append([])

df = pd.DataFrame(rows, columns=['Locus', 'Allele', 'Bracketed Repeat Region', 'Flanking Region Variants from GRCh38', 'Counts', 'Frequencies', '5\'-Flanking Region', '3\'-Flanking Region'])

csv_file_path = 'data/output_data.csv'
df.to_csv(csv_file_path, index=False)