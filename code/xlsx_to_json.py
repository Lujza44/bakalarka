import pandas as pd
import json
import re

xlsx_file_path = 'data/414_2021_2685_MOESM8_ESM.xlsx'
df = pd.read_excel(xlsx_file_path, sheet_name = 'Table S6', engine='openpyxl') # precitanie konnkretneho sheetu = Table S6, odtial budem cerpat referecncne sekvencie
df.iloc[:, 0] = df.iloc[:, 0].ffill() # forward fill na vyplnenie stlpca s nazvami markerov tam kde su merged bunky

ref_seq_dict = {}


def extract_number(chromosome_str):
    # Use regex to find all digits in the string
    numbers = re.findall(r'\d+', chromosome_str)
    # Join the digits and convert to an integer
    if numbers:
        return int(''.join(numbers))
    else:
        return None

# spracovanie dat z ref. databazy
for index, row in df.iterrows():
    if row.iloc[2] == 'Reference sequence': # budem pracovat iba s riadkami, ktore obsahuju ref. sekvenciu
        ref_before = "" # toto bude ref. 5' flanking region
        ref_seq = ""    # toto bude ref. sekvencia repet. oblasti
        ref_after = ""  # toto bude ref. 3' flanking region
        str_size = 0    # dlzka jedneho STR
        chromosome = 0
        start_coordinate = 0 # TODO zatial nevyuzita suradnica zaciatku repetitivnej oblasti v ref. genome

        if index + 3 < len(df) and index - 1 >= 0: 
            coordinates = df.iloc[index + 3] # budem sa pozerat do riadka o 3 nizsie, tam su koordinaty jednotlivych nt z ref. genomu
            str_indexes = df.iloc[index - 1] # budem sa pozerat do riadka o 1 vyssie, tam je cislovanie repeticii
            numeric_values = pd.to_numeric(str_indexes, errors='coerce').dropna()
            number_of_repetitions = int(numeric_values.max())

            is_storing_ref_seq = False
            increment_str_size = False

            for i in range(3, len(row)):
                if pd.notnull(pd.to_numeric(str_indexes.iloc[i], errors='coerce')) and str(str_indexes.iloc[i]).strip() == '1':
                    is_storing_ref_seq = True # ak sme narazili na zaciatok repetitivnej oblasti, ukladame nt do ref_seq
                    increment_str_size = True # a pocitame vzdialenost do dalsieho oznaceneho STR
                    start_coordinate = coordinates.iloc[i] # TODO zatial nevyuzita suradnica zaciatku repetitivnej oblasti v ref. genome
                    chromosome = extract_number(df.iloc[index - 1, 2])

                if pd.notnull(pd.to_numeric(str_indexes.iloc[i], errors='coerce')) and str(str_indexes.iloc[i]).strip() == '2':
                    increment_str_size = False

                if pd.notnull(pd.to_numeric(coordinates.iloc[i], errors='coerce')): # check, ci nt ma koordinatu v ref. genome, teda nie je inzercia
                    if is_storing_ref_seq:
                        ref_seq += str(row.iloc[i])
                    else:
                        ref_before += str(row.iloc[i])
                    
                    if increment_str_size:
                        str_size += 1

                     
        ref_after = ref_seq[str_size*number_of_repetitions:] 
        ref_seq = ref_seq[:str_size*number_of_repetitions] 
        
        ref_seq_dict[row.iloc[0].replace(" ", "")] = (chromosome, str_size, start_coordinate, number_of_repetitions, ref_seq, ref_before, ref_after)

# ulozenie dat z ref. databazy do jsnu
with open('data/transformed_data.json', 'r') as json_file:
    data = json.load(json_file)

for key, value in ref_seq_dict.items():
    chromosome, str_size, start_coordinate, allele, sequence, before, after = value # rozbalenie n-tice    
    data["markers"][key]["chromosome"] = chromosome
    data["markers"][key]["STRsize"] = str_size
    data["markers"][key]["startCoordinate"] = start_coordinate
    data["markers"][key]["referenceAllele"] = {
        "allele": allele,
        "sequence": sequence,
        "before": before,
        "after": after
    }

with open('data/transformed_data.json', 'w') as json_file:
    json.dump(data, json_file, indent=4)