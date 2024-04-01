import json

wrong_after_flank = ['D1S1656', 'D5S818', 'D7S820', 'vWA', 'D13S317', 'D18S51', 'D19S433', 'D21S11'] # TODO pozor D19 aj aj !!!
wrong_before_flank = ['PentaD', 'D19S433']
# TODO special penta E

json_file_path = 'data/transformed_data.json'

def process_markers(data, wrong_after_flank):
    for marker_name in wrong_after_flank:
        marker = data["markers"].get(marker_name)
        if not marker:
            continue

        STRsize = marker["referenceAllele"]["STRsize"]
        
        for variant in marker["alleleVariants"]:
            allele_str = str(variant["allele"])
            if '.' in allele_str:
                integer_part, decimal_part = allele_str.split('.')
                integer_part = int(integer_part)
                decimal_part = int(decimal_part[0])  # Assuming only one digit after decimal
            else:
                integer_part = int(allele_str)
                decimal_part = 0

            if integer_part > 0:
                theSize = (integer_part * STRsize) + decimal_part
            
                for sequenceVariant in variant["sequenceVariants"]:
                    sequence = sequenceVariant["sequence"]
                    if len(sequence) > theSize:
                        remaining_sequence = sequence[theSize:]
                        sequenceVariant["sequence"] = sequence[:int(theSize)]
                        
                        for flankingVariant in sequenceVariant["flankingRegionsVariants"]:
                            flankingVariant["after"] = remaining_sequence + flankingVariant["after"]

with open(json_file_path, 'r') as file:
    data = json.load(file)

process_markers(data, wrong_after_flank)

with open(json_file_path, 'w') as file:
    json.dump(data, file, indent=4)