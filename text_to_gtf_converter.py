#converts text to gtf file
import argparse
import csv

def convert_to_gtf(input_file, output_file, source_name="AdenineRichFinder"):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Skip the header line
        next(infile)

        # Create a CSV reader
        reader = csv.reader(infile, delimiter='\t')

        for row in reader:
            chromosome, start, end, sequence, adenine_count = row

            # GTF fields
            seqname = chromosome
            source = source_name
            feature = "adenine_rich_region"
            score = "."
            strand = "."
            frame = "."

            # Create attributes
            attributes = f'sequence "{sequence}"; adenine_count "{adenine_count}"'

            # Write the GTF line
            gtf_line = f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n"
            outfile.write(gtf_line)


def main():
    parser = argparse.ArgumentParser(description="Convert adenine-rich sequences text file to GTF format")
    parser.add_argument("input_file", help="Path to the input text file with adenine-rich sequences")
    parser.add_argument("output_file", help="Path to the output GTF file")
    parser.add_argument("--source", default="AdenineRichFinder",
                        help="Source name for the GTF file (default: AdenineRichFinder)")

    args = parser.parse_args()

    convert_to_gtf(args.input_file, args.output_file, args.source)

    print(f"Conversion complete. GTF file written to {args.output_file}")


if __name__ == "__main__":
    main()