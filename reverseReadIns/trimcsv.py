"""
T-Sequence Trimmer

This script processes large CSV files containing sequence data by trimming leading 'T' 
nucleotides from each sequence. It works efficiently by reading and processing the file 
in chunks to minimize memory usage.

Inputs:
    - input_file: Path to the input CSV file containing a 'Sequence' column
                  Default: 'fastq_analysis_output.csv'
    
    - output_file: Path where the processed CSV will be saved
                   Default: 'JK85_read1s_trimmed.csv'
    
    - chunk_size: Number of rows to process at once (optional)
                  Default: 100,000 rows

Outputs:
    - A CSV file with identical structure to the input file, but with all leading 'T' 
      nucleotides removed from the 'Sequence' column
    
    - Progress information printed to the console showing the processing status

Usage:
    python trim_t_sequences.py

Notes:
    - The script preserves all columns and rows from the original file
    - Processing is done in memory-efficient chunks
    - Progress updates display the number of rows processed and percentage complete
"""

import pandas as pd

def trim_t_sequence(sequence):
    # Remove leading T's
    trimmed = sequence.lstrip('T')
    return trimmed

def process_csv_in_chunks(input_file, output_file, chunk_size=100000):
    # Get total number of rows for progress tracking
    total_rows = sum(1 for _ in open(input_file)) - 1  # subtract 1 for header
    processed_rows = 0
    
    # Create a CSV writer for the first chunk (will create the file)
    first_chunk = True
    
    # Process the file in chunks
    for chunk_number, chunk in enumerate(pd.read_csv(input_file, chunksize=chunk_size), 1):
        # Trim T's from the Sequence column
        chunk['Sequence'] = chunk['Sequence'].apply(trim_t_sequence)
        
        # Write to output file
        if first_chunk:
            chunk.to_csv(output_file, index=False, mode='w')
            first_chunk = False
        else:
            # Append subsequent chunks without header
            chunk.to_csv(output_file, index=False, mode='a', header=False)
        
        # Update progress
        processed_rows += len(chunk)
        progress = (processed_rows / total_rows) * 100
        print(f"Processed chunk {chunk_number}: {processed_rows:,}/{total_rows:,} rows ({progress:.2f}%)")

# Usage
if __name__ == "__main__":
    input_file = 'fastq_analysis_output.csv'
    output_file = 'JK85_read1s_trimmed.csv'
    print("Starting CSV processing...")
    process_csv_in_chunks(input_file, output_file)
    print("Processing complete! Output saved to:", output_file)