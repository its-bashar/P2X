import requests
import os
import json
import pandas as pd

# Directory to save the downloaded files
output_dir = './h5ad_files/'
os.makedirs(output_dir, exist_ok=True)

# Load dataset IDs from CSV file
dataset_ids = pd.read_csv('unique_datasets.csv', header=None)[0].tolist()

domain_name = "cellxgene.cziscience.com"
api_url_base = f"https://api.{domain_name}"

def fetch_and_download_h5ad(dataset_id):
    try:
        # Formulate request and fetch a list of published Dataset Versions metadata
        dataset_versions_path = f"/curation/v1/datasets/{dataset_id}/versions"
        url = f"{api_url_base}{dataset_versions_path}"
        res = requests.get(url=url)
        res.raise_for_status()
        res_content = res.json()

        # Extract the assets from the response
        assets = res_content[0]['assets']

        # Print the URLs and related information for verification
        print(f"Dataset ID: {dataset_id}")
        print("Dataset Assets URLs:")
        for asset in assets:
            asset_url = asset['url']
            asset_filetype = asset.get('filetype', 'unknown')  # Example: type of asset, if available
            asset_filesize = asset.get('filesize', 'unknown')  # Example: size of asset, if available
            asset_filename = os.path.basename(asset_url)
            print(f"Filetype: {asset_filetype}, Filesize: {asset_filesize}, URL: {asset_url}, Filename: {asset_filename}")

        # Function to download a file from a URL
        def download_file(url, output_path):
            response = requests.get(url, stream=True)
            response.raise_for_status()
            with open(output_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
            print(f"Downloaded {output_path}")

        # Download only H5AD assets
        for asset in assets:
            if asset.get('filetype') == 'H5AD':
                asset_url = asset['url']
                asset_filename = os.path.basename(asset_url)
                output_path = os.path.join(output_dir, asset_filename)
                download_file(asset_url, output_path)
    except Exception as e:
        print(f"An error occurred while processing dataset ID {dataset_id}: {e}")

# Loop through each dataset ID and process
for dataset_id in dataset_ids:
    fetch_and_download_h5ad(dataset_id)

