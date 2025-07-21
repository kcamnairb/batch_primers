
---

# Batch Primers 

This repository provides Python scripts for designing PCR and qPCR primers for gene knockout and RT-qPCR experiments. The scripts leverage `primer3` and `Biopython` libraries to automate primer design for multiple genes based on GFF and FASTA input files.

---

## Features

1. **Primer Design for Gene Knockout (KO):**
   - Automates the design of fusion PCR primers for creating gene knockouts.
   - Generates primers optimized for `pyrG` and `ptrA` selection markers.
   - Ensures proper primer placement based on gene coordinates and strand orientation.

2. **RT-qPCR Primer Design:**
   - Designs primers for RT-qPCR experiments.
   - Includes functionality to create primers overlapping exon-exon junctions for genes with multiple exons.
   - Provides detailed primer quality checks, including specificity and placement.

3. **Error Handling:**
   - Scripts handle errors gracefully, logging problematic gene entries for review.
   - Outputs diagnostic information to help refine input data or parameters.

4. **Output Formats:**
   - Primer designs are saved in multiple formats:
     - CSV files for easy review and downstream processing.
     - GFF files for visualizing primer locations in genome browsers.

---

## Installation

### Prerequisites

1. **Python 3.x**  
   Ensure you have Python 3.x installed.

2. **Required Python Libraries**  
   Install the required libraries using `pip`:
   ```bash
   pip install primer3-py biopython
   ```

3. **Files Required for Input:**
   - GFF file: Gene annotation file containing coordinates and feature information.
   - FASTA file: Genome sequence file corresponding to the GFF file.

---

## Usage

### 1. Gene Knockout Primer Design (`cds_KO_primers.py`)

This script generates primers for fusion PCR to knock out all genes in a genome.

#### **Command:**
```bash
python cds_KO_primers.py <gff_file> <fasta_file>
```

#### **Input:**
- `<gff_file>`: GFF file containing gene annotations.
- `<fasta_file>`: FASTA file containing the genome sequence.

#### **Output:**
- `*_KO_primers_pyrG_selection.csv`: Primer designs for `pyrG` selection.
- `*_KO_primers_ptrA_selection.csv`: Primer designs for `ptrA` selection.
- `genes_with_errors`: List of genes where primer design failed.
- `primers.gff`: GFF file showing primer locations.

---

### 2. RT-qPCR Primer Design (`qpcr_primers_overlapping_exons.py`)

This script generates primers for RT-qPCR, including primers overlapping exon-exon junctions for genes with multiple exons.

#### **Command:**
```bash
python qpcr_primers_overlapping_exons.py <gff_file> <fasta_file>
```

#### **Input:**
- `<gff_file>`: GFF file containing gene annotations.
- `<fasta_file>`: FASTA file containing the genome sequence.

#### **Output:**
- `*_qPCR_primers.csv`: Primer designs with details on product size, penalties, and specificity.
- `genes_with_errors`: List of genes where primer design failed.
- `primers.gff`: GFF file showing primer locations.

---

## Output Description

1. **CSV Files:**
   - Contain primer sequences, penalties, product sizes, and additional metadata for each primer.
   - Example columns:
     - `primer name`: Name of the primer.
     - `primer sequence`: Primer sequence.
     - `penalty`: Primer3 penalty score indicating primer quality.
     - `product size`: Expected PCR product size.

2. **GFF Files:**
   - Annotates primer locations in the genome for visualization in genome browsers.
   - Example fields:
     - `Name`: Primer identifier (e.g., `gene-5F`).
     - `strand`: Indicates forward (`+`) or reverse (`-`) strand.

3. **Error Logs:**
   - `genes_with_errors`: Text file listing genes where primer design encountered issues.

---
