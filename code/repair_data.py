import json

wrong_strand = ['D1S1656', 'D2S1338', 'FGA', 'D5S818', 'CSF1PO', 'D7S820', 'vWA', 'PentaE', 'D19S433']

wrong_flanks = ['D1S1656', 'D5S818', 'D7S820', 'vWA', 'D13S317', 'D18S51', 'D19S433', 'D21S11', 'PentaD']


def complement_dna(s):
    trans_table = str.maketrans('ATCG', 'TAGC')
    return s.translate(trans_table)

def reverse_string(s):
    return s[::-1]

json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)

for marker in data['markers']:
    if marker in wrong_strand:
        ref_allele = data['markers'][marker]['referenceAllele']
        
        seq = ref_allele['sequence']
        bef = ref_allele['before']
        aft = ref_allele['after']

        ref_allele['sequence'] = reverse_string(complement_dna(seq))
        ref_allele['before'] = reverse_string(complement_dna(aft))
        ref_allele['after'] = reverse_string(complement_dna(bef))

with open(json_file_path, 'w') as file:
    json.dump(data, file, indent=4)