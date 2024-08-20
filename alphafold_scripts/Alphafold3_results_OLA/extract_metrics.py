import os
import json
import pandas as pd

# Base directory containing the folders
base_dir = "Alphafold3_results"

# List to store extracted data
data = []

# Function to extract metrics from JSON file
def extract_metrics(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    metrics = {
        'iptm': data.get('iptm', None),
        'ptm': data.get('ptm', None),
        'ranking_score': data.get('ranking_score', None),
        'fraction_disordered': data.get('fraction_disordered', None),
        'has_clash': data.get('has_clash', None),
        'num_recycles': data.get('num_recycles', None),
        'chain_iptm': data.get('chain_iptm', []),
        'chain_pair_iptm': data.get('chain_pair_iptm', []),
        'chain_pair_pae_min': data.get('chain_pair_pae_min', []),
        'chain_ptm': data.get('chain_ptm', []),
    }
    return metrics

# Traverse through directories
for folder in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder)
    if os.path.isdir(folder_path):
        # Extract subunit information and type from folder name
        parts = folder.split('_')
        if len(parts) == 6:
            first_subunit = parts[1]
            second_subunit = parts[2]
            third_subunit = parts[3]
            structure_type = parts[4]

            # JSON file path
            json_file = f"{folder}_summary_confidences_0.json"
            json_path = os.path.join(folder_path, json_file)

            # Extract metrics
            if os.path.exists(json_path):
                metrics = extract_metrics(json_path)
                # Add additional information to metrics
                metrics.update({
                    'first_subunit': first_subunit,
                    'second_subunit': second_subunit,
                    'third_subunit': third_subunit,
                    'structure_type': structure_type
                })
                data.append(metrics)

# Convert collected data to DataFrame
df = pd.DataFrame(data)

# Save DataFrame to CSV
output_csv = 'P2X_Models_Metrics.csv'
df.to_csv(output_csv, index=False)

print(f"Metrics extracted and saved to {output_csv}")