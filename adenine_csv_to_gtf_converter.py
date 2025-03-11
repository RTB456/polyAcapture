#"This script takes a tab-delimited CSV file containing
#adenine sequence data, shifts the genomic coordinates backward by
# 300 base pairs (or a user-specified amount), and converts the data into GTF format
# for visualization in IGV. It preserves the sequence and adenine count information
# in the attributes field, filters out positions that would become negative after adjustment,
# and provides detailed statistics about the processing
#results including warnings if no output is generated."

import argparse
import csv

def adjust_and_convert_to_gtf(input_file, output_file, adjustment=300, source_name="AdenineRichFinder"):
    total_lines = 0
    adjusted_lines = 0
    skipped_lines = 0

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')

        # Skip the header
        next(reader)

        for row in reader:
            total_lines += 1
            if len(row) < 5:  # Skip malformed lines
                skipped_lines += 1
                continue

            chromosome, start, end, sequence, adenine_count = row

            # Convert start and end to integers and adjust
            new_start = int(start) - adjustment
            new_end = int(end) - adjustment

            # Only include if new_start is positive
            if new_start > 0:
                # Prepare GTF fields
                seqname = chromosome
                source = source_name
                feature = "adenine_rich_region"
                score = "."
                strand = "."
                frame = "."
                attributes = f'sequence "{sequence}"; adenine_count "{adenine_count}"'

                # Write the GTF line
                gtf_line = f"{seqname}\t{source}\t{feature}\t{new_start}\t{new_end}\t{score}\t{strand}\t{frame}\t{attributes}\n"
                outfile.write(gtf_line)

                adjusted_lines += 1
            else:
                skipped_lines += 1

    return total_lines, adjusted_lines, skipped_lines


def main():
    parser = argparse.ArgumentParser(description="Adjust adenine-rich sequence ranges and convert to GTF format")
    parser.add_argument("input_file", help="Path to the input adenine-rich sequences file")
    parser.add_argument("output_file", help="Path to the output GTF file")
    parser.add_argument("--adjustment", type=int, default=300, help="Number of base pairs to subtract (default: 300)")
    parser.add_argument("--source", default="AdenineRichFinder",
                        help="Source name for the GTF file (default: AdenineRichFinder)")

    args = parser.parse_args()

    total, adjusted, skipped = adjust_and_convert_to_gtf(args.input_file, args.output_file, args.adjustment,
                                                         args.source)

    print(f"Conversion and adjustment complete. GTF file written to {args.output_file}")
    print(f"Total lines processed: {total}")
    print(f"Lines adjusted and written: {adjusted}")
    print(f"Lines skipped: {skipped}")

    if adjusted == 0:
        print("\nWARNING: No lines were written to the output file.")
        print("This could be because all start positions were <= 300.")
        print(f"Consider using a smaller adjustment value (current: {args.adjustment}).")


if __name__ == "__main__":
    main()