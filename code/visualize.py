import pandas as pd
import json
import pysam

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

# toto appenduje konkretne SNPs do pola
def find_SNPs(s1, s2, coordinate, before, STRsize, allele, chromosome):
    diff_positions = []
    for i in range(min(len(s1),len(s2))):
        if s1[i] != s2[i]:
            if before: # before flanking region
                diff_positions.append((str(chromosome), coordinate - (len(s1) - i)))
            else: # after flanking region
                diff_positions.append((str(chromosome), i + coordinate + STRsize * allele)) # allele z ref. allele
    return diff_positions
    

def get_rs_number(vcf_path, chromosome, position):
    vcf_in = pysam.VariantFile(vcf_path)
    rs_number = None

    for record in vcf_in.fetch(chromosome, position -1, position):
        rs_number = record.id

    vcf_in.close()
    return rs_number

vcf_path = 'data/00-common_all.vcf.gz'

'''
for chromosome, position in diff_positions:
    rs_number = get_rs_number(vcf_path, chromosome, position)
    print(f"Position {chromosome}:{position} has rs number(s): {'No SNPs found' if not rs_number else rs_number}")
'''

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

    str_length = marker_info['STRlength']
    repeats = marker_info['repeats']

    reference_allele = marker_info.get('referenceAllele', {})
    chromosome = marker_info.get('chromosome', 0)

    # toto sa bude porovnavat
    ref_before = reference_allele.get('before', '')
    ref_after = reference_allele.get('after', '')

    # parametre
    start_coordinate = marker_info.get('startCoordinate', '')
    ref_allele_number = reference_allele.get('numberOfRepeats', '')

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

                    SNPs_before = find_SNPs(ref_before, before, start_coordinate, True, str_length, ref_allele_number, chromosome) # zoznam
                    SNPs_after = find_SNPs(ref_after, after, start_coordinate, False, str_length, ref_allele_number, chromosome) # zoznam

                    if SNPs_before:
                        for snp in SNPs_before:
                            SNPs.append((chromosome, snp))
                    if SNPs_after:
                        for snp in SNPs_after:
                            SNPs.append((chromosome, snp))

                    SNPs = SNPs_before + SNPs_after

                    rs_numbers = []
                    for chromosome, position in SNPs:
                        rs_num = get_rs_number(vcf_path, chromosome, position)
                        if rs_num is not None: 
                            rs_numbers.append(rs_num)
                    rs_numbers = ', '.join(rs_numbers)
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