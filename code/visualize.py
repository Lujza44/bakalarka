import pandas as pd
import json
import re

'''
def to_bracket_notation(sequence, str_length, allele):
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

def to_bracket_notation(dna_sequence, substring_size, repeats):
    result = ""
    i = 0
    while i < len(dna_sequence):
        # Check if the remaining sequence is shorter than the substring size
        if i + substring_size > len(dna_sequence):
            # Append the remaining sequence and a space, then break
            result += dna_sequence[i:] + " "
            break
        # Extract the current substring of the specified size
        current_substring = dna_sequence[i:i+substring_size]
        # If the current substring is in the repeats list, append it to result and skip its length
        if current_substring in repeats:
            result += current_substring + " "
            i += substring_size
        else:
            # If not in repeats, append the current letter and a space, then move to the next letter
            result += dna_sequence[i] + " "
            i += 1
    return result.strip()

def to_bracket_notation(dna_sequence, substring_size, repeats):
    result = ""
    i = 0
    current_subsequence = ""
    count = 1  # Initialize count to 1 to correctly count the first occurrence

    while i < len(dna_sequence):
        if i + substring_size > len(dna_sequence):
            if dna_sequence[i:] == current_subsequence:
                count += 1
            else:
                result += (current_subsequence if current_subsequence else dna_sequence[i:]) + str(count) + " "
            break

        substring = dna_sequence[i:i+substring_size]
        if substring in repeats:
            if substring == current_subsequence:
                count += 1
            else:
                if current_subsequence:  # Check if it's not the first subsequence
                    result += current_subsequence + str(count) + " "
                current_subsequence = substring
                count = 1
            i += substring_size
        else:
            if current_subsequence:
                result += current_subsequence + str(count) + " "
                current_subsequence = ""
                count = 1
            else:
                result += dna_sequence[i] + "1 "
            i += 1

    if current_subsequence:  # For the last subsequence if it matches
        result += current_subsequence + str(count) + " "

    return result.strip()
'''

def to_bracket_notation(dna_sequence, substring_size, repeats):
    result = ""
    i = 0
    current_subsequence = ""
    count = 0
    single_letters = ""  # To accumulate single letters

    while i < len(dna_sequence):
        if i + substring_size > len(dna_sequence):
            if current_subsequence and current_subsequence == dna_sequence[i:]:
                count += 1
            if current_subsequence:
                result += (current_subsequence if count == 1 else "[" + current_subsequence + "]" + str(count)) + " "
            elif single_letters:  # Append accumulated single letters
                result += single_letters + dna_sequence[i:] + " "
            else:
                result += dna_sequence[i:] + " "
            break

        substring = dna_sequence[i:i+substring_size]
        if substring in repeats:
            if substring == current_subsequence:
                count += 1
            else:
                if current_subsequence:
                    result += (current_subsequence if count == 1 else "[" + current_subsequence + "]" + str(count)) + " "
                elif single_letters:  # Append accumulated single letters before the subsequence
                    result += single_letters + " "
                    single_letters = ""
                current_subsequence = substring
                count = 1
            i += substring_size
        else:
            if current_subsequence:
                result += (current_subsequence if count == 1 else "[" + current_subsequence + "]" + str(count)) + " "
                current_subsequence = ""
                count = 0
            single_letters += dna_sequence[i]  # Accumulate single letters
            i += 1

    # For the last subsequence or single letters, if any
    if current_subsequence:
        result += (current_subsequence if count == 1 else "[" + current_subsequence + "]" + str(count)) + " "
    elif single_letters:
        result += single_letters + " "

    return result.strip()


def chromosome_key(chromosome_str):
    """Extracts numerical part from chromosome string for sorting."""
    # This regular expression matches the first sequence of digits in the string
    match = re.search(r'\d+', chromosome_str)
    if match:
        return int(match.group())
    return 0  # Default value for strings without numbers

json_file_path = 'data/transformed_data.json'
with open(json_file_path, 'r') as file:
    data = json.load(file)

rows = [[]]

for marker, marker_info in sorted(data['markers'].items(), key=lambda x: chromosome_key(x[1]['referenceAllele'].get('chromosome', '0'))): # iterovanie cez vsetky markery v jsne
    sum_count = 0
    sum_freq = 0

    str_length = marker_info['STRsize']
    repeats = marker_info['repeats']

    if not str_length or not repeats:
        continue

    for allele_var in marker_info['alleleVariants']: # iterovanie cez vsetky varianty alel
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