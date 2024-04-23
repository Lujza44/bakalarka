def ORIGINAL_to_bracket_notation(dna_sequence, substring_size, repeats):
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


def test_align(dna_sequence, substring_size, repeats):
    result = []
    part = []
    i = 0
    while i < len(dna_sequence):
        if i + substring_size > len(dna_sequence):
            alined_rest =  align_partial_rep(list(current_substring), list(dna_sequence[i:]))
            result += alined_rest
            break
        current_substring = dna_sequence[i:i+substring_size]
        if current_substring in repeats:
            result += current_substring
            i += substring_size
        else:
            part += dna_sequence[i]
            next_substring = dna_sequence[i+1:i+1+substring_size]
            if next_substring in repeats:
                aligned = align_partial_rep(list(next_substring), part)
                result += aligned
            i += 1
    return result

    
def align_partial_rep(next_substring, part):
    result = [""] * len(next_substring)  # Initialize result list with empty strings
    align_iter = iter(part)  # Create an iterator over align_list
    next_elem = next(align_iter, None)  # Get the first element from the iterator

    for i, item in enumerate(next_substring):
        if item == next_elem:
            result[i] = item  # Place the item at the correct position
            next_elem = next(align_iter, None)  # Move to the next element in align_list
    return result


sequence = "CCTATCTATCATCTA"
repeats = ["CCTA", "TCTA"]
STRsize = 4

sequence = "TCTATCTATCTATCTATC"
repeats = ["CCTA", "TCTA"]
STRsize = 4

print(test_align(sequence,STRsize,repeats))