#This script processes an existing GTF file by shifting all genomic coordinates backward by 300 base pairs
# (or a user-specified amount),
# keeping the file in GTF format. It maintains all original GTF fields and attributes
# without modification, filters out any entries that would have negative start positions after adjustment,
# and outputs a simple completion message without detailed statistics.

import argparse
import csv

def adjust_gtf_ranges(input_file, output_file, adjustment=300):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t', lineterminator='\n')

        for row in reader:
            if len(row) < 9:  # Skip malformed lines
                continue

            seqname, source, feature, start, end, score, strand, frame, attributes = row

            # Convert start and end to integers and adjust
            new_start = int(start) - adjustment
            new_end = int(end) - adjustment

            # Only include if new_start is positive
            if new_start > 0:
                # Update the row with new start and end
                row[3] = str(new_start)
                row[4] = str(new_end)

                # Write the adjusted row
                writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Adjust GTF file ranges by subtracting a fixed value")
    parser.add_argument("input_file", help="Path to the input GTF file")
    parser.add_argument("output_file", help="Path to the output adjusted GTF file")
    parser.add_argument("--adjustment", type=int, default=300, help="Number of base pairs to subtract (default: 300)")

    args = parser.parse_args()

    adjust_gtf_ranges(args.input_file, args.output_file, args.adjustment)

    print(f"Adjustment complete. New GTF file written to {args.output_file}")


if __name__ == "__main__":
    main()