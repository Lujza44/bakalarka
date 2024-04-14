import json

'''
import requests

def get_rs_number(chromosome, position):
    """Query Ensembl REST API for rs number at a specific genomic position."""
    server = "https://rest.ensembl.org"
    # Adjusted to use the overlap region endpoint, specifying variation as the feature type
    ext = f"/overlap/region/human/{chromosome}:{position}-{position}?feature=variation"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(f"{server}{ext}", headers=headers)
        response.raise_for_status()
    except requests.exceptions.HTTPError as err:
        print(f"HTTP error occurred: {err}")
        return None
    except Exception as err:
        print(f"An error occurred: {err}")
        return None
    
    data = response.json()
    
    # Extract rs numbers if present
    rs_numbers = [item['id'] for item in data if 'id' in item]
    
    return rs_numbers
'''

json_file_path = 'data/transformed_data.json'

with open(json_file_path, 'r') as file:
    data = json.load(file)


diff_positions = []

# toto appenduje konkretne SNPs do pola
def find_SNPs(s1, s2, coordinate, before, STRsize, allele, chromosome):
    for i in range(min(len(s1),len(s2))):
        if s1[i] != s2[i]:
            if before: # before flanking region
                diff_positions.append((str(chromosome), coordinate - (len(s1) - i)))
            else: # after flanking region
                diff_positions.append((str(chromosome), i + coordinate + STRsize * allele)) # allele z ref. allele

# toto mi hlada SNPs
for marker, details in data['markers'].items():
    print(marker)
    if marker == 'D19S433': continue
    reference_allele = details.get('referenceAllele', {})
    chromosome = details.get('chromosome', 0)

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
                find_SNPs(ref_before, before, start_coordinate, True, str_size, ref_allele, chromosome)
                find_SNPs(ref_after, after, start_coordinate, False, str_size, ref_allele, chromosome)


print(diff_positions)


import pysam

def get_rs_number(vcf_path, chromosome, position):
    vcf_in = pysam.VariantFile(vcf_path)
    rs_number = None

    for record in vcf_in.fetch(chromosome, position -1, position):
        rs_number = record.id

    vcf_in.close()
    return rs_number

vcf_path = 'data/00-common_all.vcf.gz'

for chromosome, position in diff_positions:
    rs_number = get_rs_number(vcf_path, chromosome, position)
    print(f"Position {chromosome}:{position} has rs number(s): {'No SNPs found' if not rs_number else rs_number}")



'''
positions = [('10', '129294243'), ('12', '12297179')]

for chromosome, position in diff_positions:
    rs_numbers = get_rs_number(chromosome, position)
    print(f"Position {chromosome}:{position} has rs number(s): {'No SNPs found' if not rs_numbers else ', '.join(rs_numbers)}")
'''