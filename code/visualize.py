import pandas as pd
import json
import re
import requests


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

def find_SNPs(s1, s2, coordinate, before, STRsize, allele):
    diff_positions = []
    for i in range(min(len(s1),len(s2))):
        if s1[i] != s2[i]:
            if before: # before flanking region
                diff_positions.append(coordinate - (len(s1) - i))
            else: # after flanking region
                diff_positions.append(i + coordinate + STRsize * allele) # allele z ref. allele
    return diff_positions
    
def get_rs_number(chromosome, position):
    """Query Ensembl REST API for rs number at a specific genomic position."""
    server = "https://rest.ensembl.org"
    # Adjusted to use the overlap region endpoint, specifying variation as the feature type
    ext = f"/overlap/region/human/{chromosome}:{position}-{position}?feature=variation"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(f"{server}{ext}", headers=headers)
        response.raise_for_status()
    except requests.exceptions.HTTPError as err:
        print(f"HTTP error occurred: {err}")
        return None
    except Exception as err:
        print(f"An error occurred: {err}")
        return None
    
    data = response.json()
    
    # Extract rs numbers if present
    rs_numbers = [item['id'] for item in data if 'id' in item]
    
    return rs_numbers


json_file_path = 'data/transformed_data.json'
#json_file_path = 'data/repaired.json'
with open(json_file_path, 'r') as file:
    data = json.load(file)

rows = [[]]

# iterovanie cez vsetky markery ZORADENE podla 'chromosome'
for marker, marker_info in sorted(data['markers'].items(), key=lambda x: x[1].get('chromosome', 0)):
    if marker == 'D19S433': continue
    sum_count = 0
    sum_freq = 0

    str_length = marker_info['STRsize']
    repeats = marker_info['repeats']

    reference_allele = marker_info.get('referenceAllele', {})
    chromosome = marker_info.get('chromosome', 0)

    # toto sa bude porovnavat
    ref_before = reference_allele.get('before', '')
    ref_after = reference_allele.get('after', '')

    # parametre
    start_coordinate = marker_info.get('startCoordinate', '')
    ref_allele_number = reference_allele.get('allele', '')

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
                    rs_numbers = ""

                    SNPs_before = find_SNPs(ref_before, before, start_coordinate, True, str_length, ref_allele_number) # zoznam
                    SNPs_after = find_SNPs(ref_after, after, start_coordinate, False, str_length, ref_allele_number) # zoznam

                    if SNPs_before:
                        for snp in SNPs_before:
                            SNPs.append((marker, chromosome, snp))
                    if SNPs_after:
                        for snp in SNPs_after:
                            SNPs.append((marker, chromosome, snp))

                    SNPs = SNPs_before + SNPs_after

                    for snp in SNPs:
                        rs_numbers = ', '.join(get_rs_number(chromosome, snp))
                        #print(f"Position {chromosome}:{snp} has rs number(s): {'No SNPs found' if not rs_numbers else rs_numbers}")
                        # TODO co s No SNPs found?

                    count = flank_var['count']
                    frequency = flank_var['frequency']
                    sum_count += count
                    sum_freq += frequency
                    rows.append([marker, allele, sequence, rs_numbers, count, frequency, before, after])
            else: # ak neexistuju, bunky ostanu prazdne
                rows.append([marker, allele, sequence, "", "", "", "", ""])
    rows.append(["", "", "", "", sum_count, sum_freq, "", ""])
    rows.append([])

df = pd.DataFrame(rows, columns=['Locus', 'Allele', 'Bracketed Repeat Region', 'Flanking Region Variants from GRCh38', 'Counts', 'Frequencies', '5\'-Flanking Region', '3\'-Flanking Region'])

csv_file_path = 'data/output_data.csv'
df.to_csv(csv_file_path, index=False)

#excel_file_path = 'data/output_data.xlsx'
#df.to_excel(excel_file_path, index=False, engine='openpyxl')