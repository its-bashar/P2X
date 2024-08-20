import os
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Function to read JSON files and extract scores
def extract_scores(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return {
        'iptm': data['iptm'],
        'ptm': data['ptm'],
        'ranking_score': data['ranking_score'],
        'fraction_disordered': data['fraction_disordered']
    }

# Collect data
data = {}
counter = 0  # Initialize the counter variable
for seed in range(1, 11):
    data[seed] = {}
    for ranking in range(5):
        counter += 1
        file_path = f'fold_p2x1_seed_{seed}/fold_p2x1_seed_{seed}_summary_confidences_{ranking}.json'
        if os.path.exists(file_path):
            data[seed][ranking] = extract_scores(file_path)
print(f"Total files being compared: {counter}")

# Prepare data for analysis
scores = {
    'iptm': [],
    'ptm': [],
    'ranking_score': [],
    'fraction_disordered': []
}

for seed in data:
    for ranking in data[seed]:
        for score_type in scores:
            scores[score_type].append(data[seed][ranking][score_type])

# Calculate statistics
stats = {}
for score_type in scores:
    stats[score_type] = {
        'mean': np.mean(scores[score_type]),
        'std': np.std(scores[score_type])
    }

# Print statistics
for score_type, stat in stats.items():
    print(f"{score_type}:")
    print(f"  Mean: {stat['mean']:.4f}")
    print(f"  Std Dev: {stat['std']:.4f}")
    print()

# Visualization
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Distribution of Scores Across Seeds')

for i, (score_type, values) in enumerate(scores.items()):
    row = i // 2
    col = i % 2
    sns.boxplot(data=values, ax=axs[row, col])
    axs[row, col].set_title(score_type)
    axs[row, col].set_ylabel('Score')

plt.tight_layout()
plt.show()
