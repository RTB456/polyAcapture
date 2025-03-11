# polyAcapture
### Mapping PolyA Sequences (./mac239analysis/polyAcapture)

### Overview
This section describes the workflow for identifying and visualizing polyadenine (PolyA) sequences in viral genomes, accounting for the 300bp offset inherent in single-cell sequencing.
### Tools
* **polyA.py** (formerly HeavyAregions.py): Identifies adenine-rich regions in genome FASTA files
* **gtf_coordinate_shifter.py** (formerly adjustGTFregion.py): Adjusts coordinates in existing GTF files
* **adenine_tsv_to_gtf_converter.py** (formerly adjustAdenineRichRegions.py): Converts adenine-rich region data to GTF format with coordinate adjustment
* **text_to_gtf_converter.py**: Converts tab-delimited adenine sequence data to GTF format
* **polyC.py**: Additional tool for identifying cytosine-rich regions (similar to polyA.py)

### Workflow Steps
1. **Prepare Reference Data**
   * Load Mac239 FASTA file in IGV as reference genome
   * Upload Mac239 GTF file for gene annotations

2. **Gather Alignment Files**
   * Obtain mac239.bam and mac239.bam.bai files from scPathoQuant (viral genome aligner)
   * Note: CellRanger BAM files can be used but require additional processing (see note below)

3. **Identify Adenine-Rich Regions**
   * Run **polyA.py** on the Mac239 genome FASTA file
   * The script uses a sliding window approach (default: 20bp) to find regions with high adenine content (default threshold: 10 adenines)
   * Output is a tab-delimited file with coordinates and sequence information

4. **Adjust Coordinates for Single-Cell Sequencing Offset**
   * Option 1: Use **gtf_coordinate_shifter.py** to adjust an existing GTF file by shifting coordinates 300bp upstream
   * Option 2: Use **adenine_csv_to_gtf_converter.py** to process the csv adenine-rich regions file, adjust coordinates, and output in GTF format
   * These adjustments account for the 300bp offset caused by polyT primer binding in single-cell sequencing

5. **Visualize Results**
   * Load the adjusted GTF file in IGV
   * These annotations mark sequences 300bp upstream of PolyA sequences in the reference genome
   * Visually inspect read distribution to identify reads aligned 300bp downstream of polyT primer binding sites

### Note on Using CellRanger BAM Files
To use BAM files from CellRanger instead of scPathoQuant:
1. Determine the exact coordinates where your viral genome starts in the combined reference
2. Use a BAM manipulation tool (samtools or PySam in Python) to extract reads that align to those specific coordinates and beyond

## Processing & Analyzing Reads

This section describes tools for processing raw sequencing data, preparing it for analysis, and aligning it to reference genomes.

### Tools Overview

* **process_fastq.py**: Processes FASTQ files from single-cell RNA-seq, extracting cell barcodes, UMIs, and sequences
* **trim_t_sequences.py**: Removes leading T nucleotides from sequences in large CSV files
* **align_hiv_sequences.py**: Aligns processed sequences to an HIV reference genome using Bowtie2
* **plot_alignment_coverage.py**: Generates coverage depth plots from alignment results

### Read Processing Workflow

1. **Extract Sequence Components** (process_fastq.py)
   - Parses gzipped FASTQ files and extracts 16bp cell barcodes, 12bp UMIs, and sequence content
   - Outputs a CSV file with sequence information and quality scores
   - Processes data in memory-efficient chunks

2. **Trim Leading T's** (trim_t_sequences.py)
   - Removes leading T nucleotides that may result from polyA capture
   - Preserves all other sequence information
   - Handles large files by processing in chunks

3. **Align to Reference** (align_hiv_sequences.py)
   - Creates reverse complements of sequences
   - Aligns to HIV reference genome using Bowtie2
   - Filters for high-quality alignments (MAPQ ≥ 20)
   - Generates comprehensive alignment statistics

4. **Visualize Coverage** (plot_alignment_coverage.py)
   - Creates coverage depth plots from alignment results
   - Parses CIGAR strings to accurately map reads to reference positions
   - Outputs publication-quality visualization

### ReverseReadIns

The `align_hiv_sequences.py` script implements a reverse-read insertion workflow, which is critical for properly aligning sequences from polyA-captured RNA:

1. **Reverse Complement Creation**
   - Single-cell RNA-seq with polyA capture produces reads in the opposite orientation to the reference
   - The script creates reverse complements of all sequences before alignment
   - This ensures proper alignment to the reference genome

2. **Sequence Validation**
   - Checks for and filters out invalid sequences
   - Reports on sequences where reverse complement operations fail

3. **FASTA Generation**
   - Creates temporary FASTA files with properly formatted headers containing cell and UMI information
   - Uses the reverse complemented sequences for alignment

4. **Bowtie2 Alignment**
   - Uses local alignment mode with high sensitivity
   - Parameters optimized for viral genome mapping

This reverse-read approach is necessary because polyA capture protocols bind to the 3' end of mRNAs, resulting in reads that run in the opposite direction of the reference genome.