import json

def to_bracket_notation(sequence, str_length):
    result = []
    current_str = sequence[:str_length]
    count = 1
    
    for i in range(str_length, len(sequence), str_length): # prechadzam sekvenciou po jednotlivych nt
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


with open('data/transformed_data.json', 'r') as json_file:
    data = json.load(json_file)

    