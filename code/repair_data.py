
####### TOTO CELE BY SA asi MOHLO PRESUNUT DO xlsx_to_json.py #################

import json

wrong_strand = ['D1S1656', 'D2S1338', 'FGA', 'D5S818', 'CSF1PO', 'D7S820', 'vWA', 'PentaE', 'D19S433'] # DONE

wrong_after_flank = ['D1S1656', 'D5S818', 'D7S820', 'vWA', 'D13S317', 'D18S51', 'D19S433', 'D21S11'] # TODO pozor D19 aj aj !!!
wrong_before_flank = ['PentaD', 'D19S433']

def complement_dna(s):
    trans_table = str.maketrans('ATCG', 'TAGC')
    return s.translate(trans_table)

def reverse_string(s):
    return s[::-1]

json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)

for marker, info in data['markers'].items():
    ref_allele = info['referenceAllele']
    
    if marker in wrong_strand:
        
        seq = ref_allele['sequence']
        bef = ref_allele['before']
        aft = ref_allele['after']

        ref_allele['sequence'] = reverse_string(complement_dna(seq))
        ref_allele['before'] = reverse_string(complement_dna(aft))
        ref_allele['after'] = reverse_string(complement_dna(bef))

'''
    for allele_variant in info['alleleVariants']:
        if allele_variant['sequenceVariants']:
            first_variant = allele_variant['sequenceVariants'][0]
            if first_variant['flankingRegionsVariants']:
                first_flank = first_variant['flankingRegionsVariants'][0]
                if first_flank['before']:
                    length_before = len(first_flank['before'])
                    ref_allele['before'] = ref_allele['before'][-length_before:]
                if first_flank['after']:
                    length_after = len(first_flank['after'])
                    ref_allele['after'] = ref_allele['after'][:length_after]
                break 
'''
    # here will be the reduction of flanking region sequences, because it needs to be done for every marker
    # iterate through allele variants of the marker, take the first one that has flanking region variants that is NOT an empty
    # list and the flanking variants are NOT empty strings. take the length of the flanking region before (let's call it LENGHT now)
    # and from the reference allele cut the before sequence so that only the end of it, that is LENGHT letters long, remains.
    # so something like before = before[:LENGHT]. then do a similar thing with after string, only now the beginning remains.
    # so take the lenght of the flanking region after (LENGHT1) of the first allele variant that has one and do
    # after = after[LENGHT1:]
        
    # now that i have unified the strand and cut the flanking regions of each marker i can repair the allele variants where
    # there is a bit of flanking region in the repetitive region of the sequence (and it should not be here)
    # check if the str_size * allele (the number of allele variant) isn't < than actual length of the sequence of the allele variant
    # if it is, then there is something extra in the sequence, don't do anything with this comment yet


with open(json_file_path, 'w') as file:
    json.dump(data, file, indent=4)