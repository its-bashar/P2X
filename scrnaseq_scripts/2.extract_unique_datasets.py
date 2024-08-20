import pandas as pd

# Load the CSV file
df = pd.read_csv('filtered_p2x_genes.csv')

# Check if all datasets are Homo sapiens
all_homo_sapiens = (df['organism_name'] == 'Homo sapiens').all()

# Group dataset IDs by publication citation
grouped = df.groupby('publication_citation')['dataset_id'].unique()

# Print out the unique publications and their unique datasets
for publication, datasets in grouped.items():
    print(f"Publication: {publication}")
    print(f"Unique Dataset IDs: {', '.join(datasets)}")
    print("-" * 40)

# Save the unique Homo sapiens dataset IDs to a file
homo_sapiens_datasets = df[df['organism_name'] == 'Homo sapiens']['dataset_id'].unique()
pd.Series(homo_sapiens_datasets).to_csv('unique_datasets.csv', index=False, header=False)

if not all_homo_sapiens:
    print("Note: Not all datasets are Homo sapiens, but a saved file unique_datasets.csv contains only the Homo sapiens datasets.")
else:
	print("Note: all datasets are Homo sapiens, file saved as unique_datasets.csv")

