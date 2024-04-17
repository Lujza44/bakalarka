import pandas as pd
import json

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


# DATA Z JSNU
json_file_path = 'data/transformed_data.json'
with open(json_file_path, 'r') as file:
    data = json.load(file)

rows = [[]]

# iterovanie cez vsetky markery v JSNE ZORADENE podla 'chromosome'
for marker, marker_info in sorted(data['markers'].items(), key=lambda x: x[1].get('chromosome', 0)):
    if marker == 'D19S433': continue
    sum_count = 0

    str_length = marker_info['STRlength']
    repeats = marker_info['repeats']

    reference_allele = marker_info.get('referenceAllele', {})
    chromosome = marker_info.get('chromosome', 0)

    if not str_length or not repeats:
        continue

    # iterovanie cez vsetky varianty alel ZORADENE podla 'allele'
    for allele_var in sorted(marker_info['lengthVariants'], key=lambda x: x['numberOfRepeats']): 
        allele = allele_var['numberOfRepeats']
        for seq_var in allele_var['sequenceVariants']:
            sequence = seq_var['sequence']

            sequence = to_bracket_notation(seq_var['sequence'], str_length, repeats)
            
            if seq_var['flankingRegionsVariants']: # ak existuju flanking region varianty, kazdy bude mat vlastny riadok
                for flank_var in seq_var['flankingRegionsVariants']:
                    before = flank_var['before']
                    after = flank_var['after']
                    
                    before_indices = []
                    before_rs = []
                    after_indices = []
                    after_rs = []

                    if "beforeRsNumbers" in flank_var:
                        before_tuples = flank_var["beforeRsNumbers"]
                        for snp in before_tuples:
                            before_indices.append(snp[0])
                            before_rs.append(snp[1])
                    
                    if "afterRsNumbers" in flank_var:
                        after_tuples = flank_var["afterRsNumbers"]
                        for snp in after_tuples:
                            after_indices.append(snp[0])
                            after_rs.append(snp[1])

                    before_indices = ", ".join(str(num) for num in before_indices)
                    after_indices = ", ".join(str(num) for num in after_indices)

                    before_rs = [x for x in before_rs if x != '']
                    after_rs = [x for x in after_rs if x != '']
                    rs_numbers = ", ".join(before_rs + after_rs)

                    count = flank_var['count']
                    frequency = flank_var['frequency']
                    sum_count += count

                    rows.append([marker, allele, sequence, rs_numbers, count, frequency, before, before_indices, after, after_indices])
            else: # ak neexistuju, bunky ostanu prazdne
                rows.append([marker, allele, sequence])
    rows.append(["", "", "", "", sum_count])
    rows.append([])

# skonvertovanie riadkov do objektu DataFrame
df = pd.DataFrame(rows, columns=['Locus', 'Allele', 'Bracketed Repeat Region', 'Flanking Region Variants from GRCh38', 'Counts', 'Frequencies', '5\'-Flanking Region', '5\' SNP indexes', '3\'-Flanking Region', '3\' SNP indexes'])

# DataFrame zapisany do CSV
csv_file_path = 'data/vis/raw_vis.csv'
df.to_csv(csv_file_path, index=False)