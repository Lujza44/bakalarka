import pandas as pd
import json
import pysam


def find_SNPs_before(s1, s2, coordinate):
    SNPs_before_coordinates = []
    SNPs_before_indices = []
    for i in range(min(len(s1),len(s2))):
        if s1[i] != s2[i]:
            SNPs_before_coordinates.append(coordinate - (len(s1) - i)) 
            SNPs_before_indices.append(i)
    return SNPs_before_coordinates, SNPs_before_indices # koordinaty, indexy

def find_SNPs_after(s1, s2, coordinate, STRsize, allele):
    SNPs_after_coordinates = []
    SNPs_after_indices = []
    for i in range(min(len(s1),len(s2))):
        if s1[i] != s2[i]:
            SNPs_after_coordinates.append(i + coordinate + STRsize * allele) 
            SNPs_after_indices.append(i)
    return SNPs_after_coordinates, SNPs_after_indices # koordinaty, indexy
    
def get_rs_number(vcf_path, chromosome, position):
    vcf_in = pysam.VariantFile(vcf_path)
    rs_number = None

    for record in vcf_in.fetch(str(chromosome), int(position) -1, int(position)):
        rs_number = record.id

    vcf_in.close()
    return rs_number

vcf_path = 'data/00-common_all.vcf.gz'

json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)


for marker, details in data['markers'].items():
    if marker == "D19S433": continue # TODO

    reference_allele = details.get('referenceAllele', {})
    
    ref_before = reference_allele.get('before', '')
    ref_after = reference_allele.get('after', '')
    start_coordinate = details.get('startCoordinate', '')
    ref_allele_number = reference_allele.get('numberOfRepeats', '')
    str_length = details['STRlength']
    chromosome = details['chromosome']

    allele_variants = details.get('lengthVariants', [])
    for variant in allele_variants:
        for sequence_variant in variant.get('sequenceVariants', []):
            flanking_variants = sequence_variant.get('flankingRegionsVariants', [])
            if flanking_variants:
                for flank_variant in flanking_variants:
                    before = flank_variant['before']
                    after = flank_variant['after']

                    SNPs_before_coordinates, SNPs_before_indices = find_SNPs_before(ref_before, before, start_coordinate)
                    SNPs_after_coordinates, SNPs_after_indices = find_SNPs_after(ref_after, after, start_coordinate, str_length, ref_allele_number)

                    rs_numbers_before = []
                    rs_numbers_after = []
                    for coordinate in SNPs_before_coordinates:
                        rs_num = get_rs_number(vcf_path, chromosome, coordinate)
                        if rs_num is not None: 
                            rs_numbers_before.append(rs_num)

                    for coordinate in SNPs_after_coordinates:
                        rs_num = get_rs_number(vcf_path, chromosome, coordinate)
                        if rs_num is not None: 
                            rs_numbers_after.append(rs_num)
                    
                    flank_variant['beforeSNPIndices'] = SNPs_before_indices
                    flank_variant['beforeRsNumbers'] = rs_numbers_before
                    flank_variant['afterSNPIndices'] = SNPs_after_indices
                    flank_variant['afterRsNumbers'] = rs_numbers_after

with open(json_file_path, 'w') as file:
    json.dump(data, file, indent=4)