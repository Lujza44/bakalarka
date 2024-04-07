import requests

def check_ensembl_api_status():
    url = "https://rest.ensembl.org/info/ping?content-type=application/json"
    response = requests.get(url)
    if response.ok:
        data = response.json()
        if data and data.get('ping') == 1:
            print("Ensembl API is reachable and responding.")
        else:
            print("Ensembl API response is unexpected:", data)
    else:
        print("Failed to reach Ensembl API.")

check_ensembl_api_status()
