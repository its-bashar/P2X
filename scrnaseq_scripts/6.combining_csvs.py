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
output_file = 'combined_p2x_filtered.csv'

# Initialize an empty DataFrame to store the combined data
combined_df = pd.DataFrame(index=p2x_gene_ids)

# Process each CSV file in the input directory
for file_name in os.listdir(input_dir):
    if file_name.endswith('_p2x_filtered.csv'):
        file_path = os.path.join(input_dir, file_name)
        try:
            # Load the CSV file into a DataFrame
            print(f"Processing file: {file_path}")
            df = pd.read_csv(file_path, index_col=0)

            # Ensure the DataFrame has the correct gene names in the correct order
            df = df.loc[p2x_gene_ids]

            # Append dataset identifier to each sample name (column name)
            dataset_id = os.path.splitext(file_name)[0]
            df.columns = [f"{col}_{dataset_id}" for col in df.columns]

            # Concatenate the columns of the current DataFrame to the combined DataFrame
            combined_df = pd.concat([combined_df, df], axis=1)

            # Print success message
            print(f"File processed: {file_path}")

        except Exception as e:
            # Print error message
            print(f"Error processing file {file_path}: {e}")

# Save the combined DataFrame to a CSV file
combined_df.to_csv(output_file)

# Print success message
print(f"Combined data saved to {output_file}")

