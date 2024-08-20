# P2X Heterotrimer Analysis

This repository contains scripts and analysis for exploring hetero-trimeric P2X receptor complexes using AlphaFold3 and single-cell RNA sequencing (scRNA-seq) data. The project aims to investigate the structural feasibility and co-expression patterns of P2X hetero-trimers in human cells.

## Repository Structure

### alphafold_scripts/
- **Alphafold3_results_OLA**: Scripts used for AlphaFold3 prediction results including oleic acids
- **Alphafold3_results_trimmed**: Scripts used for AlphaFold3 prediction results for trimmed sequences
- **Alphafold3_seed_exploring**: Scripts for seed variation analysis of AlphaFold3 Models
- **MSA**: Multiple Sequence Alignment files and scripts for UniProt p2x receptors
- **P2X3_AlphaFold_PDB_Comparison**: Scripts for comparing AlphaFold predictions with known P2X3 crystal structures
- **P2X3_OLA_exploring**: Scripts for investigating oleic acid influence on P2X3 structure
- **P2X5_MSA**: Multiple Sequence Alignment specific to P2X5

### scrnaseq_scripts/
1. `expression_filter_p2x_genes.sh`: Filters expression data for P2X genes
2. `extract_unique_datasets.py`: Extracts unique datasets IDs from the database
3. `download_all_files.py`: Downloads necessary files from the database
4. `h5ad_files_p2x_filtering_converting_to_csvs.py`: Filters h5ad files for P2X genes and converts to CSVs
5. `csv_files_p2x_filtered_cleaning_and_sorting.py`: Cleans and sorts filtered CSV files
6. `combining_csvs.py`: Combines multiple CSV files into one CSV file
7. `check_duplicates.py`: Checks for duplicate samples (this script only checks the duplicates)
8. `checking_and_removing_duplicates.R`: R script for checking and removing duplicates
9. `flipping_rows_with_columns.py`: Transposes the data matrix
10. `Filtering.ipynb`: Jupyter notebook for data filtering
11. `metadata_combiner.py`: Extracts and Combines metadata of the filtered samples from the .h5ead files
12. `Analysis.ipynb`: Jupyter notebook for main data analysis
13. `test_submit_job.sbatch`: SLURM job submission script for HPC test file

## Usage

### AlphaFold3 Analysis
1. Download and Prepare sequences using MSA scripts in the `MSA` and `P2X5_MSA` folders.
2. Run AlphaFold3 predictions and use these scripts to analyise the results using `Alphafold3_results_trimmed` for trimmed sequences and `Alphafold3_results_OLA` for including oleic acids.
3. (optional) To Analyze seed variation using scripts in `Alphafold3_seed_exploring`.
4. Compare results with known structures using MM-align scripts in `P2X3_AlphaFold_PDB_Comparison`.

### scRNA-seq Analysis
1. Run `expression_filter_p2x_genes.sh` to filter for P2X genes.
2. Use `extract_unique_datasets.py` and `download_all_files.py` to prepare the dataset.
3. Process h5ad files using `h5ad_files_p2x_filtering_converting_to_csvs.py`.
4. Clean and combine CSV files using scripts 5-7.
5. Remove duplicates with `checking_and_removing_duplicates.R`.
6. Prepare final dataset using `flipping_rows_with_columns.py` and `Filtering.ipynb`.
7. Combine metadata using `metadata_combiner.py`.
8. Perform main analysis using `Analysis.ipynb`.

For HPC usage, modify and use `test_submit_job.sbatch` as needed.

## Dependencies

- Python 3.9.18
- R (version used in HPC)
- scanpy 1.10.1
- AlphaFold3 (web server used)
- SikuliX 2.0.5 IDE
- EMBL-EBI Clustal Omega Job Dispatcher
- Jalview
- MM-Align 2021/08/16 version

## Data Sources

- CZ CELLxGENE Discover database
- UniProt Knowledgebase
- UCSC Genome Browser
- Protein Data Bank (PDB)

## Citing This Work

If you use the code or methodology from this project, please cite:

Bashar Majed Al-Smadi. (2024). Exploring Hetero-trimeric P2X Receptor Complexes with AlphaFold3 and scRNA-seq. <https://github.com/its-bashar/P2X>

## Contact

Bashar Al-Smadu - Bsmadi20@gmail.com

For any questions or feedback, please open an issue in this GitHub repository.
