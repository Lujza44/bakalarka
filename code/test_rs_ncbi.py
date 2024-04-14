import requests

'''
def get_rs_numbers(chromosome, position):
    """Retrieve rs numbers for SNPs at a specific genomic position from NCBI."""
    url = "https://api.ncbi.nlm.nih.gov/variation/v0/coords/human"
    
    # Prepare the query. Note: You might need to adjust parameters based on the API documentation.
    query = {
        "chr": chromosome,
        "pos": position,
    }
    
    try:
        response = requests.get(url, params=query)
        response.raise_for_status()  # Raises an exception for HTTP errors
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
        return None
    
    data = response.json()
    
    # Extract and return the rs numbers from the response
    rs_numbers = [variant['refsnp_id'] for variant in data['data']]
    return rs_numbers

# Example usage:
chromosome = "17"
position = 41223094  # Example position, adjust as needed
rs_numbers = get_rs_numbers(chromosome, position)
print(rs_numbers)
'''


import pysam

def get_rs_number(vcf_path, chromosome, position):
    vcf_in = pysam.VariantFile(vcf_path)
    #rs_numbers = []

    for record in vcf_in.fetch(chromosome, position -1, position):
        rs_number = record.id
        #if record.pos + 1 == position:  # Check if the position matches exactly
        #    rs_numbers.append(record.id)

    vcf_in.close()
    return rs_number

vcf_path = 'data/00-common_all.vcf.gz'
chromosome = '2'  # Example chromosome
position = 68011922  # Example position

#rs_number = get_rs_number(vcf_path, chromosome, position)
#print(rs_number)



positions = [('10', 129294243), ('12', 12297179)]

for chromosome, position in positions:
    rs_number = get_rs_number(vcf_path, chromosome, position)
    print(f"Position {chromosome}:{position} has rs number(s): {'No SNPs found' if not rs_number else rs_number}")