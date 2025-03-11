#Find polyC regions in GTFfile
import sys
from Bio import SeqIO

def find_poly_c(sequence, length=5):
    poly_c_ranges = []
    seq_str = str(sequence)
    i = 0
    while i < len(seq_str):
        if seq_str[i:i+length] == 'C' * length:
            start = i + 1  # 1-based indexing
            end = i + length
            poly_c_ranges.append((start, end))
            i += length  # Move past this poly-C sequence
        else:
            i += 1
    return poly_c_ranges

def main(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        poly_c_ranges = find_poly_c(record.seq)
        print(f"Sequence: {record.id}")
        print("Non-overlapping poly-C (CCCCC) ranges:")
        for start, end in poly_c_ranges:
            print(f"  {start}-{end}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py HIV_VIRUS.fa")
    else:
        main(sys.argv[1])