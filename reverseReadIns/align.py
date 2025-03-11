"""
HIV Sequence Alignment Pipeline

This script processes single-cell RNA-seq data to align trimmed sequences to an HIV reference genome.
It processes data in manageable chunks to handle large datasets efficiently, using Bowtie2 for alignment.

Inputs:
    - JK85_read1s_trimmed.csv: A CSV file containing previously trimmed sequences
                               Required columns: 'CellID', 'UMI', 'Sequence'
    
    - A Bowtie2 index named 'hiv_index' in the current directory

Outputs:
    - alignment_results/ directory containing:
      * final_alignment_results.csv: All high-quality alignments
      * final_alignment_stats.txt: Summary statistics of the alignment process
      * intermediate_results_chunk_*.csv: Intermediate results saved during processing

Workflow:
    1. Reads input CSV in chunks of 100,000 sequences
    2. For each chunk:
       - Validates sequences and creates reverse complements
       - Writes sequences to temporary FASTA file
       - Runs Bowtie2 alignment with --local and --very-sensitive-local settings
       - Parses SAM output and extracts high-quality alignments (MAPQ ≥ 20)
       - Cleans up temporary files
    3. Combines results and saves statistics

Requirements:
    - Python packages: numpy, pandas, biopython
    - External tools: Bowtie2 (must be in PATH)
    - Bowtie2 index files for HIV reference (named 'hiv_index')

Usage:
    python align_hiv_sequences.py
"""

import numpy as np
import pandas as pd
import subprocess
from Bio.Seq import Seq
import os
from pathlib import Path
import datetime
import math

def reverse_complement_sequence(seq):
    """
    Return the reverse complement of a DNA sequence using Biopython
    """
    try:
        if pd.isna(seq):
            return None
        seq_str = str(seq).strip()
        if not seq_str:
            return None
        return str(Seq(seq_str).reverse_complement())
    except Exception as e:
        print(f"Error processing sequence: {seq}")
        print(f"Error details: {str(e)}")
        return None

def parse_sam_file(sam_file):
    """
    Parse SAM file without using pysam
    """
    results = []
    total_reads = 0
    unmapped = 0
    low_mapq = 0
    
    print(f"Parsing SAM file: {sam_file}")
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            
            total_reads += 1
            if total_reads % 50000 == 0:  # Increased progress reporting interval
                print(f"Processed {total_reads} alignments...")
                
            fields = line.strip().split('\t')
            
            if len(fields) < 11:
                continue
                
            flag = int(fields[1])
            mapq = int(fields[4])
            
            if flag & 0x4:
                unmapped += 1
                continue
                
            if mapq < 20:
                low_mapq += 1
                continue
                
            read_id = fields[0]
            ref_pos = int(fields[3])
            cigar = fields[5]
            seq = fields[9]
            
            results.append({
                'read_id': read_id,
                'ref_pos': ref_pos,
                'mapq': mapq,
                'cigar': cigar,
                'sequence': seq
            })
    
    stats = {
        'total_reads': total_reads,
        'unmapped': unmapped,
        'low_mapq': low_mapq,
        'mapped_high_qual': len(results)
    }
    
    return results, stats

def process_chunk(chunk, chunk_num, output_dir):
    """
    Process a single chunk of the dataframe
    """
    print(f"\nProcessing chunk {chunk_num}")
    print(f"Chunk size: {len(chunk)} sequences")
    
    # Data validation
    print("Validating sequences...")
    chunk['Sequence'] = chunk['Sequence'].astype(str).str.strip()
    invalid_seqs = chunk['Sequence'].isna() | (chunk['Sequence'] == '')
    if invalid_seqs.any():
        print(f"Found {invalid_seqs.sum()} invalid sequences in chunk {chunk_num}")
        chunk = chunk[~invalid_seqs].copy()
    
    # Get reverse complement of sequences
    print("Getting reverse complements...")
    chunk['rev_comp_seq'] = chunk['Sequence'].apply(reverse_complement_sequence)
    
    # Remove rows where reverse complement failed
    valid_rev_comp = chunk['rev_comp_seq'].notna()
    if not valid_rev_comp.all():
        print(f"Removed {(~valid_rev_comp).sum()} sequences with failed reverse complement")
        chunk = chunk[valid_rev_comp]
    
    if len(chunk) == 0:
        print(f"No valid sequences remaining in chunk {chunk_num}")
        return [], {'total_reads': 0, 'unmapped': 0, 'low_mapq': 0, 'mapped_high_qual': 0}
    
    # Write sequences to temporary FASTA file
    temp_fasta = output_dir / f'temp_seqs_chunk_{chunk_num}.fa'
    print(f"Writing {len(chunk)} sequences to temporary FASTA: {temp_fasta}")
    with open(temp_fasta, 'w') as f:
        for idx, row in chunk.iterrows():
            f.write(f'>{row.CellID}_{row.UMI}\n{row.rev_comp_seq}\n')
    
    # Run bowtie2 alignment
    sam_output = output_dir / f'alignment_chunk_{chunk_num}.sam'
    print(f"Running bowtie2 alignment for chunk {chunk_num}...")
    alignment_cmd = [
        'bowtie2',
        '-x', 'hiv_index',  # Changed to use index in current directory
        '-f',
        '-U', str(temp_fasta),
        '-S', str(sam_output),
        '--local',
        '--very-sensitive-local'
    ]
    
    try:
        subprocess.run(alignment_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running bowtie2 for chunk {chunk_num}")
        print(f"Command output: {e.output}")
        print(f"Error message: {e.stderr}")
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)
        return [], {'total_reads': 0, 'unmapped': 0, 'low_mapq': 0, 'mapped_high_qual': 0}
    
    # Parse alignment results
    print(f"Parsing alignment results for chunk {chunk_num}...")
    alignment_results, stats = parse_sam_file(sam_output)
    
    # Clean up temporary files
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)
    if os.path.exists(sam_output):
        os.remove(sam_output)
    
    return alignment_results, stats

def main():
    print("Starting process...")
    start_time = datetime.datetime.now()
    print(f"Start time: {start_time}")
    
    # Create output directory
    output_dir = Path('alignment_results')
    output_dir.mkdir(exist_ok=True)
    
    # Read CSV in chunks
    chunk_size = 100000  # Increased to 100,000
    total_results = []
    total_stats = {
        'total_reads': 0,
        'unmapped': 0,
        'low_mapq': 0,
        'mapped_high_qual': 0
    }
    
    print(f"Reading CSV file in chunks of {chunk_size}")
    try:
        chunks = pd.read_csv('JK85_read1s_trimmed.csv', chunksize=chunk_size)
        
        # Get total number of chunks for progress tracking
        with open('JK85_read1s_trimmed.csv', 'r') as f:
            total_lines = sum(1 for _ in f) - 1  # Subtract 1 for header
        total_chunks = math.ceil(total_lines / chunk_size)
        print(f"Total sequences to process: {total_lines}")
        print(f"Total number of chunks: {total_chunks}")
        
        # Process each chunk
        for chunk_num, chunk in enumerate(chunks, 1):
            print(f"\nProcessing chunk {chunk_num} of {total_chunks} ({(chunk_num/total_chunks)*100:.1f}% complete)")
            
            results, stats = process_chunk(chunk, chunk_num, output_dir)
            
            # Update totals
            total_results.extend(results)
            for key in total_stats:
                total_stats[key] += stats[key]
            
            # Save intermediate results
            if chunk_num % 2 == 0:  # Save every 2 chunks since chunks are bigger now
                print("Saving intermediate results...")
                pd.DataFrame(total_results).to_csv(
                    output_dir / f'intermediate_results_chunk_{chunk_num}.csv',
                    index=False
                )
        
        # Save final results
        print("\nSaving final results...")
        results_df = pd.DataFrame(total_results)
        results_df.to_csv(output_dir / 'final_alignment_results.csv', index=False)
        
        # Save final statistics
        end_time = datetime.datetime.now()
        with open(output_dir / 'final_alignment_stats.txt', 'w') as f:
            f.write(f"Alignment started: {start_time}\n")
            f.write(f"Alignment completed: {end_time}\n")
            f.write(f"Total processing time: {end_time - start_time}\n\n")
            f.write(f"Total reads processed: {total_stats['total_reads']}\n")
            f.write(f"Unmapped reads: {total_stats['unmapped']}\n")
            f.write(f"Low mapping quality reads: {total_stats['low_mapq']}\n")
            f.write(f"High quality mapped reads: {total_stats['mapped_high_qual']}\n")
        
        print("\nProcessing complete!")
        print(f"Total processing time: {end_time - start_time}")
        print(f"Results saved in: {output_dir}")
        
    except Exception as e:
        print(f"Error during processing: {str(e)}")
        raise

if __name__ == "__main__":
    main()