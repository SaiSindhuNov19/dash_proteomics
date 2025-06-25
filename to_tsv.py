import pandas as pd
import glob
import os

# Folder where your .parquet files are
input_dir = './'  # or set your folder

# Find all .parquet files
parquet_files = glob.glob(os.path.join(input_dir, '*.parquet'))

# Loop through and convert each file
for parquet_file in parquet_files:
    # Load the parquet file into a DataFrame
    df = pd.read_parquet(parquet_file)

    # Generate output .tsv file name
    base_name = os.path.splitext(parquet_file)[0]
    output_tsv = base_name + '.tsv'

    # Save as TSV
    df.to_csv(output_tsv, sep='\t', index=False)

    print(f"Converted {parquet_file} -> {output_tsv}")

print("All files converted!")

