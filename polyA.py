#takes in .fa file and outputs gtf file corresponding to POLYA sequences
from Bio import SeqIO
import argparse

#takes in input genome file, window_size (how many base pairs), minimum number of A's needed (default 10)
def find_adenine_rich_sequences(fasta_file, window_size=20, adenine_threshold=10, jump_size=None):
    if jump_size is None:
        jump_size = window_size // 2  # Default to half the window size

    adenine_rich_sequences = [] #store sequences

    #iterate through each sequence in FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq) #convert to string
        i = 0 #start at position 0

        #continue, until not enough sequence left for full window.
        while i < len(sequence) - window_size + 1:

            window = sequence[i:i + window_size] #current window
            if window.count('A') >= adenine_threshold:
                adenine_rich_sequences.append({
                    'gene': record.id,
                    'start': i + 1,  # 1-based coordinate
                    'end': i + window_size,
                    'sequence': window,
                    'adenine_count': window.count('A')
                })
                i += jump_size  # Jump by the specified amount (avoid recording similar problems)
            else:
                i += 1  # Move to the next base pair if not adenine-rich

    return adenine_rich_sequences

#help with command-line interface
def main():
    parser = argparse.ArgumentParser(description="Find Adenine-rich sequences in a genome")
    # add default arguments
    parser.add_argument("fasta_file", help="Path to the input FASTA file")
    parser.add_argument("--window-size", type=int, default=20, help="Size of the sliding window (default: 20)")
    parser.add_argument("--adenine-threshold", type=int, default=10,
                        help="Minimum number of Adenines required (default: 10)")
    parser.add_argument("--jump-size", type=int,
                        help="Number of base pairs to jump after finding a match (default: half of window size)")
    parser.add_argument("--output", default="adenine_rich_sequences.txt",
                        help="Output file name (default: adenine_rich_sequences.txt)")

    args = parser.parse_args()

    #call function
    adenine_rich_sequences = find_adenine_rich_sequences(args.fasta_file, args.window_size, args.adenine_threshold,
                                                         args.jump_size)

    with open(args.output, 'w') as f:
        #write header
        f.write("Chromosome\tStart\tEnd\tSequence\tAdenine Count\n")
        #write sequences to file
        for seq in adenine_rich_sequences:
            f.write(f"{seq['gene']}\t{seq['start']}\t{seq['end']}\t{seq['sequence']}\t{seq['adenine_count']}\n")
            print(f"Found Adenine-rich sequence at {seq['gene']}:{seq['start']}-{seq['end']}")
    print(f"\nFound {len(adenine_rich_sequences)} Adenine-rich sequences. Full results written to {args.output}")


if __name__ == "__main__":
    main()