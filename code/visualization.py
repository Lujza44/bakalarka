import csv

file_path1 = 'marker_auto_strview.csv'
file_path2 = 'marker_auto_strview_flankingreg.csv'
file_path3 = 'new.csv'


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


with open(file_path1, mode='r', encoding='utf-8') as str, \
     open(file_path2, mode='r', encoding='utf-8') as str_flankingreg, \
     open(file_path3, mode='w', encoding='utf-8', newline = '') as new_file :
     
    str = csv.reader(str)
    str_flankingreg = csv.reader(str_flankingreg)
    csv_writer = csv.writer(new_file)

    next(str, None) #preskocenie headeru
    next(str_flankingreg, None) #preskocenie headeru

    HEADER = ["locus", "allele", "bracketed repeat region", "flanking region variants", "counts (n)", "frequencies","5'-flanking region", "3'-flanking region"]
    csv_writer.writerow(HEADER)

    for row1, row2 in zip(str, str_flankingreg):

        if row1[2] == 'Null':
            continue

        if int(row1[1]) % 1 != 0:
            continue

        repetition_only = row1[2]
        whole_sequence = row2[2]
            

        # identifikacia flanking regions
        
        start_index = whole_sequence.find(repetition_only)
        end_index = start_index + len(repetition_only)

        string_before = whole_sequence[:start_index] # 5' flanking region
        string_after = whole_sequence[end_index:]    # 3' flanking region


        # TODO 

        str_length = len(repetition_only)//int(row1[1])

        repetition = to_bracket_notation(repetition_only, str_length)

        row_to_write = (row1[0],row1[1], repetition, '', row2[4], row2[5], string_before, string_after)
        csv_writer.writerow(row_to_write)