# Neutral MSA Pipeline

## Overview

This pipeline generates clean, phylogenetically-suitable neutral multiple sequence alignments (MSAs) for branch-specific positive selection testing using the branch-specific positive selection test from:

**Berrio, A., Haygood, R. & Wray, G.A.** (2020). Identifying branch-specific positive selection throughout the regulatory genome using an appropriate proxy neutral. *BMC Genomics*, 21, 359. https://doi.org/10.1186/s12864-020-6752-4

## Data Source & Preparation

Non-Functional Regions (NFRs) were obtained from the adaptiPhy repository:
- **Source:** [adaptiPhy/deprecated/NFR.sorted.bed](https://github.com/wodanaz/adaptiPhy/blob/master/deprecated/NFR.sorted.bed)

The original NFR coordinates were lifted over to human hg38 assembly (GRCh38) using UCSC LiftOver:
- **Tool:** https://genome.ucsc.edu/cgi-bin/hgLiftOver
- **Output:** hg38_NFR.bed file containing coordinates of non-functional neutral regions in hg38

## Automated Workflow

All steps are automated via `Neutral_MSAs.sh`:

### Step 1: Extract_Neutral_Bases.py
Extracts MSAs from hg38 NFR coordinates. Produces raw FASTA files containing sequences for each NFR region.

### Step 2: cleaning_fa.py
Removes duplicate sequences, low-quality entries, and sequences with excessive gaps/unknown bases. Produces cleaned FASTA files.

### Step 3: clipkit.sh
Trims alignment positions using ClipKIT (kpi-gappy mode). Removes gappy and non-informative sites while retaining phylogenetic signal. Produces `.clipkit` files.

### Step 4: Files_Stats.py
Generates summary statistics (sequence count, alignment length) for all processed alignments. Outputs results to Excel format for downstream analysis.

### Step 5: Generate_References.py

Creates concatenated neutral reference alignments by combining multiple NFR sequences horizontally. For each query alignment from the CREs pipeline, this script:

Finds suitable NFR reference files that contain all the same species
Concatenates reference files to reach ~3000 bp target alignment length
Produces neutral reference MSAs matched to the same species composition as the query

These neutral references serve as proxy alignments representing evolution under no selection, enabling statistical comparison with CRE alignments in positive selection tests.
Output: Reference alignments saved with summary statistics in reference_summary.xlsx

## Usage

```bash
./Neutral_MSAs.sh
```

## Output

- **Neutral Folder** (trimmed and ready neutral MSAs for analysis)
- **Neutral_Stats.xlsx** (sequence counts and lengths needed for next step)
- **References Folder** (concatenated neutral reference alignments paired with CREs)
- **reference_summary.xlsx** (reference construction details and statistics)

## References

- Berrio et al. (2020). Branch-specific positive selection test. https://doi.org/10.1186/s12864-020-6752-4
- adaptiPhy NFR coordinates. https://github.com/wodanaz/adaptiPhy
- UCSC LiftOver tool. https://genome.ucsc.edu/cgi-bin/hgLiftOver
- ClipKIT alignment trimming. https://github.com/JLSteenwyk/ClipKIT
