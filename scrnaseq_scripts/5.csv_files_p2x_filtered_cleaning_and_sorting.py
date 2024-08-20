import os
import pandas as pd

# List of P2X gene IDs in the desired order
p2x_gene_ids = [
    'ENSG00000108405', 'ENSG00000187848', 'ENSG00000109991',
    'ENSG00000135124', 'ENSG00000083454', 'ENSG00000099957',
    'ENSG00000089041'
]

# Directory containing the CSV files
input_dir = 'CELLxGENE_csv_files_p2x_filtered'

# Process each CSV file in the input directory
for file_name in os.listdir(input_dir):
    if file_name.endswith('_p2x_filtered.csv'):
        file_path = os.path.join(input_dir, file_name)
        try:
            # Load the CSV file into a DataFrame
            print(f"Processing file: {file_path}")
            df = pd.read_csv(file_path, index_col=0)

            # Track missing P2X genes
            missing_genes = []

            # Check for missing P2X gene names and append rows with NaN values if needed
            for gene_id in p2x_gene_ids:
                if gene_id not in df.index:
                    df.loc[gene_id] = [np.nan] * df.shape[1]
                    missing_genes.append(gene_id)

            # Print missing P2X genes if any were found and added
            if missing_genes:
                print(f"Missing P2X genes in {file_name}: {missing_genes}")
                print(f"Added missing P2X genes with NaN values.")

            # Sort the DataFrame rows by the desired order of P2X gene IDs
            df = df.loc[p2x_gene_ids]

            # Save the updated DataFrame back to the CSV file
            df.to_csv(file_path)

            # Print success message
            print(f"File processed and saved: {file_path}")

        except Exception as e:
            # Print error message
            print(f"Error processing file {file_path}: {e}")

