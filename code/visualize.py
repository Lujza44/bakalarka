import pandas as pd
import json
import re

def to_bracket_notation(sequence, str_length):
    result = []
    current_str = sequence[:str_length]
    count = 1
    
    for i in range(str_length, len(sequence), str_length): # prechadzam sekvenciou po jednotlivych n-ticiach nukleotidov
        next_str = sequence[i:i+str_length]
        
        if next_str == current_str: # ak najblizsia repeticia je rovnaka ako aktualna
            count += 1  # zvysim pocet
        else: # ak je najblizia repeticia nova, tie predtym zapisem
            if count > 1:  # v zatvorkovej notacii ak je repeticii viac
                result.append(f"[{current_str}]{count}")
            else:  # bez zatvoriek ak je iba jedna
                result.append(current_str)
            
            current_str = next_str
            count = 1
    
    if count > 1: # pre poslednu poziciu
        result.append(f"[{current_str}]{count}")
    else:
        result.append(current_str)
    
    return ' '.join(result)

def chromosome_key(chromosome_str):
    """Extracts numerical part from chromosome string for sorting."""
    # This regular expression matches the first sequence of digits in the string
    match = re.search(r'\d+', chromosome_str)
    if match:
        return int(match.group())
    return 0  # Default value for strings without numbers

def create_rows():
    rows = []

    for marker, marker_info in sorted(data['markers'].items(), key=lambda x: chromosome_key(x[1]['referenceAllele'].get('chromosome', '0'))): # iterovanie cez vsetky markery v jsne
        str_length = marker_info['referenceAllele'].get('STRsize', None)
        if not str_length:  
            continue  
        
        for allele_var in marker_info['alleleVariants']: # iterovanie cez vsetky varianty alel
            allele = allele_var['allele']
            for seq_var in allele_var['sequenceVariants']:
                sequence = seq_var['sequence']

                if allele.is_integer(): # TODO aj pre necele cisla
                    sequence = to_bracket_notation(seq_var['sequence'], str_length)
                
                if seq_var['flankingRegionsVariants']: # ak existuju flanking region varianty, kazdy bude mat vlastny riadok
                    for flank_var in seq_var['flankingRegionsVariants']:
                        before = flank_var['before']
                        after = flank_var['after']
                        count = flank_var['count']
                        frequency = flank_var['frequency']
                        rows.append([marker, allele, sequence, "", count, frequency, before, after])
                else: # ak neexistuju, bunky ostanu prazdne
                    rows.append([marker, allele, sequence, "", count, frequency, "", ""])
        rows.append([])
        return rows

json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)

df = pd.DataFrame(create_rows(), columns=['Locus', 'Allele', 'Bracketed Repeat Region', 'Flanking Region Variants from GRCh38', 'Counts', 'Frequencies', '5\'-Flanking Region', '3\'-Flanking Region'])

csv_file_path = 'data/output_data.csv'
df.to_csv(csv_file_path, index=False)