import scanpy as sc
import os
import numpy as np
import pandas as pd
import gc

# List of P2X gene IDs
p2x_gene_ids = [
    'ENSG00000108405', 'ENSG00000187848', 'ENSG00000109991',
    'ENSG00000135124', 'ENSG00000083454', 'ENSG00000099957',
    'ENSG00000089041'
]

# Directories
input_dir = 'CELLxGENE_h5ad_files'
output_dir = 'CELLxGENE_csv_files_p2x_filtered'
metadata_dir = 'CELLxGENE_csv_files_metadata'
log_file = 'processed_files.log'

# Create the output and metadata directories if they do not exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(metadata_dir):
    os.makedirs(metadata_dir)

# Load the log file to check which files have already been processed
processed_files = set()
if os.path.exists(log_file):
    with open(log_file, 'r') as f:
        processed_files = set(line.strip() for line in f)

# Processing each .h5ad file in the input directory
for file_name in os.listdir(input_dir):
    if file_name.endswith('.h5ad') and file_name not in processed_files:
        file_path = os.path.join(input_dir, file_name)
        try:
            # Loading the .h5ad file in backed mode
            print(f"Processing file: {file_path}")
            adata = sc.read_h5ad(file_path, backed='r')

            # Checking the shape of the data
            print(f"Data shape: {adata.shape}")

            # Extract the indices of the P2X genes
            p2x_gene_indices = [adata.var_names.get_loc(gene_id) for gene_id in p2x_gene_ids if gene_id in adata.var_names]
            if len(p2x_gene_indices) == 0:
                print(f"No P2X genes found in {file_path}")
                continue

            # Display the indices of P2X genes
            print("P2X gene indices:", p2x_gene_indices)

            # Calculate the number of chunks
            chunk_size = 50000
            total_chunks = (adata.n_obs // chunk_size) + 1
            print(f"Total number of chunks: {total_chunks}")

            chunk_files = []
            metadata_files = []

            # Process the data in chunks
            for i in range(total_chunks):
                start = i * chunk_size
                end = min(start + chunk_size, adata.n_obs)
                print(f"Processing chunk: rows {start} to {end}")

                # Load the current chunk into memory
                chunk = adata[start:end, :].to_memory()

                # Filter the chunk for P2X genes
                filtered_chunk = chunk[:, p2x_gene_indices]

                # Remove samples with zero counts for all P2X genes
                non_zero_indices = np.any(filtered_chunk.X.toarray() != 0, axis=1)
                filtered_chunk = filtered_chunk[non_zero_indices]

                # Convert the filtered chunk to a DataFrame and transpose
                df = filtered_chunk.to_df().T

                # Save the DataFrame as a CSV file
                chunk_file = os.path.join(output_dir, f'filtered_chunk_{i}.csv')
                df.to_csv(chunk_file)
                chunk_files.append(chunk_file)

                # Save the metadata (obs) of the chunk to a separate CSV file
                metadata_df = filtered_chunk.obs
                metadata_file = os.path.join(metadata_dir, f'metadata_chunk_{i}.csv')
                metadata_df.to_csv(metadata_file)
                metadata_files.append(metadata_file)

                # Clear memory
                del chunk
                del filtered_chunk
                del df
                del metadata_df
                gc.collect()

                # Print statements to monitor progress
                print(f"Chunk processed: {start} to {end}")
                print(f"Filtered chunk saved to {chunk_file}")
                print(f"Metadata saved to {metadata_file}")

            # Combine all chunk CSV files into a single CSV file
            combined_csv_file = os.path.join(output_dir, f"{os.path.splitext(file_name)[0]}_p2x_filtered.csv")
            combined_df = pd.concat([pd.read_csv(chunk_file, index_col=0) for chunk_file in chunk_files], axis=1)
            combined_df.to_csv(combined_csv_file)

            # Delete chunk files after combining
            for chunk_file in chunk_files:
                os.remove(chunk_file)

            # Combine all metadata CSV files into a single metadata CSV file
            combined_metadata_file = os.path.join(metadata_dir, f"{os.path.splitext(file_name)[0]}_metadata.csv")
            combined_metadata_df = pd.concat([pd.read_csv(metadata_file, index_col=0) for metadata_file in metadata_files])
            combined_metadata_df.to_csv(combined_metadata_file)

            # Delete metadata files after combining
            for metadata_file in metadata_files:
                os.remove(metadata_file)

            # Log the successfully processed file
            with open(log_file, 'a') as f:
                f.write(file_name + '\n')

            # Print success message
            print(f"Filtered data combined and saved to {combined_csv_file}")
            print(f"Metadata combined and saved to {combined_metadata_file}")

        except Exception as e:
            # Print error message
            print(f"Error processing file {file_path}: {e}")
