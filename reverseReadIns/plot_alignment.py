"""
Alignment Coverage Plotter

This script generates a coverage depth plot from an alignment results CSV file.
It calculates the read coverage at each position in the reference genome by
parsing CIGAR strings and mapping reads to their corresponding positions.

Inputs:
    - csv_file: Path to CSV file containing alignment results
                Required columns: 'ref_pos' (reference position) and 'cigar' (CIGAR string)
                Default: "final_alignment_results.csv"
    
    - start_pos: Starting position on the reference to visualize (1-based)
                 Default: 1
    
    - end_pos: Ending position on the reference to visualize
               Default: 10300

Outputs:
    - A PNG image file 'alignment_coverage.png' showing the coverage depth
      across the specified region of the reference genome
    
    - Returns the maximum coverage depth as a numeric value

Usage:
    python plot_alignment_coverage.py

Notes:
    - Coverage is calculated by interpreting CIGAR strings to determine
      the span of each read on the reference
    - Only operations that consume the reference (M, =, X) are counted
    - The output plot includes a filled area under the coverage line for better visualization
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_alignment_coverage(csv_file, start_pos=1, end_pos=10300):
    # Read the CSV file with the header row
    df = pd.read_csv(csv_file)
    
    # Convert ref_pos to integer
    df['ref_pos'] = pd.to_numeric(df['ref_pos'], errors='coerce')
    
    # Function to parse CIGAR string and get read length
    def get_read_length(cigar):
        length = 0
        current_num = ''
        for char in cigar:
            if char.isdigit():
                current_num += char
            elif char in 'MIDNSHP=X':
                if char in 'M=X':  # These operations consume both query and reference
                    length += int(current_num)
                current_num = ''
        return length
    
    # Calculate coverage
    coverage = np.zeros(end_pos - start_pos + 1)
    
    # For each read, increment the coverage array
    for _, row in df.iterrows():
        read_start = int(row['ref_pos'])
        read_length = get_read_length(row['cigar'])
        read_end = read_start + read_length
        
        # Only count positions within our window
        if read_end >= start_pos and read_start <= end_pos:
            start_idx = max(0, read_start - start_pos)
            end_idx = min(end_pos - start_pos, read_end - start_pos)
            coverage[start_idx:end_idx] += 1
    
    # Create the plot
    plt.figure(figsize=(15, 5))
    plt.plot(range(start_pos, end_pos + 1), coverage, color='blue')
    plt.fill_between(range(start_pos, end_pos + 1), coverage, alpha=0.3)
    
    # Customize the plot
    plt.xlabel('Reference Position')
    plt.ylabel('Coverage Depth')
    plt.title('Read Coverage Distribution')
    plt.grid(True, alpha=0.3)
    
    # Add some padding to the y-axis
    plt.margins(y=0.1)
    
    # Save the plot
    plt.savefig('alignment_coverage.png', dpi=300, bbox_inches='tight')
    plt.close()
    return max(coverage)  # Return maximum coverage for reference

# Usage
if __name__ == "__main__":
    csv_file = "final_alignment_results.csv"
    max_coverage = plot_alignment_coverage(csv_file)
    print(f"Maximum coverage depth: {max_coverage}")