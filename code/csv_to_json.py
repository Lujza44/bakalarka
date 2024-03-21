import pandas as pd
import json

# funkcia, ktora transformuje data z .csv formatu do .json formatu
def transform_to_json(data):
    json_structure = {"markers": {}}
    
    for _, row in data.iterrows(): # iteracia vsetkymi riadkami .csv suboru
        marker = row['marker']
        allele = row['allele']
        sequence = row['sequence'] if row['sequence'] != "Null" else ""
        count = row['count_seq']
        frequency = row['frequency']

        flank_sequence = row['flank_sequence'] if row['flank_sequence'] != "Null" else ""
        flank_count = row['flank_count_seq']
        flank_frequency = row['flank_frequency']


        
        if marker not in json_structure['markers']: # inicializujem marker ak este v strukture nie je
            #json_structure['markers'][marker] = []
            #json_structure['markers'][marker] = {"ref_seq": "", "alleles": []}
            json_structure['markers'][marker] = {
                "referenceAllele": {}, 
                "alleleVariants": []}
        
        # najdenie alebo inicializacia alely
        alleles = json_structure['markers'][marker]["alleleVariants"]
        allele_entry = next((a for a in alleles if a['allele'] == allele), None)
        if allele_entry is None:
            allele_entry = {
                "allele": allele, 
                "sequenceVariants": []}
            alleles.append(allele_entry)
        
        # najdenie alebo inicializacia sekvencnej varianty alely
        variants = allele_entry['sequenceVariants']
        variant_entry = next((v for v in variants if v['sequence'] == sequence), None)
        if variant_entry is None:
            variant_entry = {
                "sequence": sequence,
                "count": count,
                "frequency": frequency,
                "flankingRegionsVariants": []
            }
            variants.append(variant_entry)
        
        # pridanie flanking region variant, ak existuju
        if pd.notnull(flank_sequence) and pd.notnull(flank_count) and pd.notnull(flank_frequency):
            allele_index = flank_sequence.find(sequence)
            before = flank_sequence[:allele_index]  # Substring before the allele sequence
            after = flank_sequence[allele_index + len(sequence):]            
            variant_entry['flankingRegionsVariants'].append({
                "before": before,
                "after": after,
                "count": int(flank_count),
                "frequency": flank_frequency
            })
    
    return json_structure



csv_file_path = 'data/sql_fin_data.csv'
data = pd.read_csv(csv_file_path) # data z .csv

json_data = transform_to_json(data) # .csv -> .json

json_file_path = 'data/transformed_data.json'
with open(json_file_path, 'w') as f:
    json.dump(json_data, f, indent=4) # ulozenie dat do .json suboru