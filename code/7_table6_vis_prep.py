import pandas as pd
import json


json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)

rows = [[]]

# iterovanie cez vsetky markery v JSNE ZORADENE podla 'chromosome'
for marker, marker_info in sorted(data['markers'].items(), key=lambda x: x[1].get('chromosome', 0)):
    if marker == 'D19S433': continue

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
            
            if seq_var['flankingRegionsVariants']: # ak existuju flanking region varianty, kazdy bude mat vlastny riadok
                for flank_var in seq_var['flankingRegionsVariants']:
                    before = flank_var['before']
                    after = flank_var['after']

                    before_indexes = flank_var["beforeSNPIndices"]
                    before_rs = flank_var["beforeRsNumbers"]
                    after_indexes = flank_var["afterSNPIndices"]
                    after_rs = flank_var["afterRsNumbers"]

                    rs_numbers = ", ".join(before_rs + after_rs)
                    before_indexes = ", ".join(str(num) for num in before_indexes)
                    after_indexes = ", ".join(str(num) for num in after_indexes)

                    rows.append([marker, allele, chromosome, sequence, rs_numbers, before, before_indexes, after, after_indexes])
            else: # ak neexistuju, bunky ostanu prazdne
                rows.append([marker, allele, sequence])
    rows.append([])

df = pd.DataFrame(rows, columns=['Locus', 'Allele', 'Bracketed Repeat Region', 'Flanking Region Variants from GRCh38', 'Counts', 'Frequencies', '5\'-Flanking Region', '5\' SNP indexes', '3\'-Flanking Region', '3\' SNP indexes'])

csv_file_path = 'data/raw_vis6.csv'

df.to_csv(csv_file_path, index=False)