"""
FASTQ Processor for Cell Barcodes and UMIs

This script processes gzipped FASTQ files from single-cell RNA-seq experiments,
extracting cell barcodes, UMIs, and sequence information. It works efficiently
with large files by processing data in memory-managed chunks.

Inputs:
    - input_file: Path to the input FASTQ.gz file
                  Default: "combined_R1.fastq.gz" in the current directory
    
    - output_file: Path where the processed data will be saved as a CSV
                   Default: "fastq_analysis_output.csv"
    
    - chunk_size: Number of reads to process before writing to disk (optional)
                  Default: 100,000 reads

Outputs:
    - A CSV file with four columns:
      * CellID: 16bp cell barcode sequence from the start of each read
      * UMI: 12bp Unique Molecular Identifier following the cell barcode
      * Sequence: The remaining sequence after cell barcode and UMI
      * Average_Quality: Mean Phred quality score for the read
    
    - Progress information printed to the console during processing

Usage:
    python process_fastq.py

Notes:
    - The script assumes reads follow the structure: [16bp CellID][12bp UMI][Sequence]
    - Processing happens in chunks to minimize memory usage
    - Garbage collection is triggered after each chunk to manage memory
    - Progress updates include timestamp information
"""

import os
import gzip
import pandas as pd
from Bio import SeqIO
import gc
from datetime import datetime

def process_fastq(input_file, output_file, chunk_size=100000):
    """
    Process FASTQ file and extract cellID, UMI, and sequence information
    
    Args:
        input_file (str): Path to input FASTQ.gz file
        output_file (str): Path to output CSV file
        chunk_size (int): Number of reads to process in each chunk
    """
    print(f"Starting processing at {datetime.now().strftime('%H:%M:%S')}")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    
    # Initialize counters
    total_reads = 0
    current_chunk = 0
    
    # Parse the FASTQ file
    fastq_parser = SeqIO.parse(gzip.open(input_file, "rt"), format="fastq")
    
    # Initialize data dictionary
    data = []
    
    try:
        for i, record in enumerate(fastq_parser):
            # Extract information (assuming same structure as your original script)
            cell_id = str(record.seq[:16])
            umi = str(record.seq[16:28])
            sequence = str(record.seq[28:])  # Rest of the sequence
            quality = record.letter_annotations["phred_quality"]
            avg_quality = sum(quality) / len(quality)
            
            # Store in data list
            data.append({
                'CellID': cell_id,
                'UMI': umi,
                'Sequence': sequence,
                'Average_Quality': round(avg_quality, 2)
            })
            
            total_reads += 1
            
            # Process chunk if chunk_size is reached
            if (i + 1) % chunk_size == 0:
                current_chunk += 1
                
                # Convert chunk to DataFrame
                df_chunk = pd.DataFrame(data)
                
                # Write to CSV (with header only for first chunk)
                if current_chunk == 1:
                    df_chunk.to_csv(output_file, mode='w', index=False)
                else:
                    df_chunk.to_csv(output_file, mode='a', header=False, index=False)
                
                # Print progress
                print(f"Processed {total_reads:,} reads ({current_chunk} chunks) at {datetime.now().strftime('%H:%M:%S')}")
                
                # Clear data and collect garbage
                data = []
                gc.collect()
    
        # Process any remaining reads
        if data:
            df_chunk = pd.DataFrame(data)
            df_chunk.to_csv(output_file, mode='a', header=(current_chunk == 0), index=False)
            print(f"Processed remaining {len(data):,} reads")
    
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        raise
    
    finally:
        # Print summary statistics
        print("\nProcessing Complete!")
        print(f"Total reads processed: {total_reads:,}")
        print(f"Total chunks processed: {current_chunk + 1}")
        print(f"Final output saved to: {output_file}")
        print(f"Finished at {datetime.now().strftime('%H:%M:%S')}")

if __name__ == "__main__":
    # Set up input and output paths
    current_dir = os.getcwd()
    input_file = "combined_R1.fastq.gz"
    output_file = "fastq_analysis_output.csv"
    
    # Make sure input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found!")
    
    # Process the file
    process_fastq(input_file, output_file)