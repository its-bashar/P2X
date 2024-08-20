import pandas as pd
import re
import os
import gc
from multiprocessing import Pool, cpu_count
from functools import partial

# Load the filtered data
filtered_df = pd.read_csv('filtered_final_count_matrix.csv')

# Define the function to extract dataset ID
def extract_dataset_id(cell_name):
    match = re.search(r'_([^_]+)_p2x_filtered$', cell_name)
    return match.group(1) if match else None

# Define the function to extract sample ID
def extract_sample_id(cell_name):
    return re.sub(r'_[^_]+_p2x_filtered$', '', cell_name)

# Path to the CELLxGENE_csv_files_metadata folder
metadata_folder = 'CELLxGENE_csv_files_metadata'

# Path to the output CSV file
output_csv = 'combined_metadata.csv'

# Path to the processed samples log file
log_file = 'processed_samples.log'

# Initialize the CSV file with headers if it does not exist
if not os.path.exists(output_csv):
    headers = [
        'sample', 'tissue_ontology_term_id', 'cell_type_ontology_term_id', 'cell_type',
        'disease_ontology_term_id', 'sex_ontology_term_id', 'self_reported_ethnicity_ontology_term_id',
        'development_stage_ontology_term_id', 'assay'
    ]
    with open(output_csv, 'w') as f:
        f.write(','.join(headers) + '\n')

# Load the log file if it exists
if os.path.exists(log_file):
    with open(log_file, 'r') as f:
        processed_samples = set(f.read().splitlines())
else:
    processed_samples = set()

# Extract cell IDs from the filtered DataFrame
cell_ids = filtered_df.iloc[:, 0]

# Create a Series of dataset IDs
dataset_ids = cell_ids.apply(extract_dataset_id)

# Group cell IDs by dataset IDs
grouped_cell_ids = cell_ids.groupby(dataset_ids)

def process_dataset(dataset_group, metadata_folder, output_csv, log_file, processed_samples):
    dataset_id, cells = dataset_group
    print(f'Processing dataset_id: {dataset_id}')
    
    # Path to the metadata file for the current dataset_id
    metadata_file = os.path.join(metadata_folder, f'{dataset_id}_metadata.csv')
    
    if not os.path.exists(metadata_file):
        print(f'Metadata file not found for dataset_id: {dataset_id}')
        return
    
    # Load the metadata file
    metadata_df = pd.read_csv(metadata_file)
    
    # Rename the first column to 'sample'
    metadata_df.rename(columns={metadata_df.columns[0]: 'sample'}, inplace=True)
    
    # Drop duplicate columns if they exist
    metadata_df = metadata_df.loc[:, ~metadata_df.columns.duplicated()]
    
    # Process each sample in the current dataset
    for cell_id in cells:
        sample_id = extract_sample_id(cell_id)
        
        # Skip processing if the sample has already been processed
        if sample_id in processed_samples:
            print(f'Skipping already processed sample_id: {sample_id}')
            continue
        
        print(f'Processing sample_id: {sample_id}')
        
        # Extract the relevant subset of metadata
        subset = metadata_df[
            metadata_df['sample'] == sample_id
        ][[
            'sample', 'tissue_ontology_term_id', 'cell_type_ontology_term_id', 'cell_type',
            'disease_ontology_term_id', 'sex_ontology_term_id', 'self_reported_ethnicity_ontology_term_id',
            'development_stage_ontology_term_id', 'assay'
        ]]
        
        if subset.empty:
            print(f'No metadata found for sample_id: {sample_id}')
        else:
            print(f'Metadata found for sample_id: {sample_id}, appending to CSV')
            # Append the subset to the CSV file
            subset.to_csv(output_csv, mode='a', header=False, index=False)
            
            # Log the processed sample
            with open(log_file, 'a') as log:
                log.write(f'{sample_id}\n')
    
    # Delete the metadata DataFrame from memory
    del metadata_df
    gc.collect()

if __name__ == '__main__':
    # Determine the number of processes to use (up to 60)
    num_processes = min(60, cpu_count())
    
    # Create a partial function with fixed arguments
    process_dataset_partial = partial(
        process_dataset,
        metadata_folder=metadata_folder,
        output_csv=output_csv,
        log_file=log_file,
        processed_samples=processed_samples
    )
    
    # Create a pool of worker processes
    with Pool(num_processes) as pool:
        # Map the process_dataset function to the grouped cell IDs
        pool.map(process_dataset_partial, grouped_cell_ids)

    print(f'Finished processing. Combined metadata saved to {output_csv}')
