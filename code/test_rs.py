import json


def find_SNPs(s1, s2, coordinate, before, STRsize, allele):
    diff_positions = []
    for i in range(min(len(s1),len(s2))):
        if s1[i] != s2[i]:
            if before: # before flanking region
                #diff_positions.append(coordinate - (len(s1) - i))
                diff_positions.append(i)
            else: # after flanking region
                diff_positions.append(i + coordinate + STRsize * allele) # allele z ref. allele
                #diff_positions.append(i)
            #print(s1[i],s2[i])
    if diff_positions != []: print(diff_positions)



json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)



for marker, details in data['markers'].items():
    print(marker)
    if marker == 'D19S433': continue
    reference_allele = details.get('referenceAllele', {})

    # toto sa bude porovnavat
    ref_before = reference_allele.get('before', '')
    ref_after = reference_allele.get('after', '')

    # parametre
    str_size = details.get('STRsize', 4)
    start_coordinate = details.get('startCoordinate', '')
    ref_allele = reference_allele.get('allele', '')


    allele_variants = details.get('alleleVariants', [])
    for variant in allele_variants:
        for sequence_variant in variant.get('sequenceVariants', []):
            flanking_variants = sequence_variant.get('flankingRegionsVariants', [])
            for flank_var in flanking_variants:
                before = flank_var.get('before', '')
                after = flank_var.get('after', '')
                find_SNPs(ref_before, before, start_coordinate, True, str_size, ref_allele)
                find_SNPs(ref_after, after, start_coordinate, False, str_size, ref_allele)
    print()



#with open(json_file_path, 'w') as file:
#    json.dump(data, file, indent=4)