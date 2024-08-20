import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def process_data(data, trim_residues=None):
    plddts = np.array(data['atom_plddts'])
    pae = np.array(data['pae'])
    
    # Count unique residues
    residues = list(zip(data['token_chain_ids'], data['token_res_ids']))
    unique_residues = len(set(residues))
    
    if trim_residues:
        plddts = plddts[:trim_residues]
        pae = pae[:trim_residues, :trim_residues]
        unique_residues = min(unique_residues, trim_residues)
    
    avg_plddt = np.mean(plddts)
    avg_pae = np.mean(pae)

    return avg_plddt, avg_pae, plddts, pae, unique_residues

# Load data
ola_data = load_json('fold_p2x3_ola_full_data_0.json')
trimmed_data = load_json('fold_p2x3_p2x3_p2x3_trimmed_full_data_0.json')

# Process data
ola_avg_plddt, ola_avg_pae, ola_plddts, ola_pae, ola_num_residues = process_data(ola_data, trim_residues=1065)
trimmed_avg_plddt, trimmed_avg_pae, trimmed_plddts, trimmed_pae, trimmed_num_residues = process_data(trimmed_data)

# Print results
print(f"OLA structure - Average pLDDT: {ola_avg_plddt:.2f}, Average PAE: {ola_avg_pae:.2f}, Number of Residues: {ola_num_residues}")
print(f"Trimmed structure - Average pLDDT: {trimmed_avg_plddt:.2f}, Average PAE: {trimmed_avg_pae:.2f}, Number of Residues: {trimmed_num_residues}")
# Visualizations
# 1. Bar plot for average pLDDT scores
plt.figure(figsize=(10, 6))
plt.bar(['OLA structure', 'Trimmed structure'], [ola_avg_plddt, trimmed_avg_plddt])
plt.title('Comparison of Average pLDDT Scores')
plt.ylabel('Average pLDDT')
plt.savefig('plddt_comparison.png')
plt.close()

# 2. PAE heatmaps
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

sns.heatmap(ola_pae, ax=ax1, cmap='viridis')
ax1.set_title('PAE Heatmap - OLA structure')
ax1.set_xlabel('Residue')
ax1.set_ylabel('Residue')

sns.heatmap(trimmed_pae, ax=ax2, cmap='viridis')
ax2.set_title('PAE Heatmap - Trimmed structure')
ax2.set_xlabel('Residue')
ax2.set_ylabel('Residue')

plt.tight_layout()
plt.savefig('pae_heatmaps.png')
plt.close()

# 3. pLDDT distribution plot
plt.figure(figsize=(12, 6))
plt.hist(ola_plddts, bins=50, alpha=0.5, label='OLA structure')
plt.hist(trimmed_plddts, bins=50, alpha=0.5, label='Trimmed structure')
plt.title('Distribution of pLDDT Scores')
plt.xlabel('pLDDT Score')
plt.ylabel('Frequency')
plt.legend()
plt.savefig('plddt_distribution.png')
plt.close()
