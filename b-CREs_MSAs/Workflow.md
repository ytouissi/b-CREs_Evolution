# Primate Sequence Alignment Pipeline

Iterative pipeline to extract and optimize multiple sequence alignments (MSAs) for cis-regulatory elements (CREs) across primate species.

## Installation

### Install Python Dependencies

```bash
pip install biopython requests openpyxl
```

### Install clipkit

Install clipkit via conda:

```bash
conda install -c bioconda clipkit
```

### Install PHAST (for msa_view)

Install PHAST via conda:

```bash
conda install -c bioconda phast
```

Verify the installation:

```bash
which msa_view
```

This command will output the full path to msa_view on your system. Copy this path as you will need it for the next step.

### Configuring msa_view in extraction scripts

Open the following files in a text editor:
- `fa_extractor_cactus447.py`
- `fa_extractor_multiz470.py`
- `fa_extractor_ucsc30.py`

In each file, locate this line:

```python
msa_view_cmd = ["/usr/bin/msa_view", temp_maf_path, "--out-format", "FASTA", "--seqs", SEQ_LIST]
```

Replace `/usr/bin/msa_view` with the path you obtained from the `which msa_view` command. For example:

```python
msa_view_cmd = ["/home/user/miniconda3/bin/msa_view", temp_maf_path, "--out-format", "FASTA", "--seqs", SEQ_LIST]
```

Update this line in all three extraction scripts.

## Required Input Files

- `CREs.bed` - BED file with region coordinates
- `primates.txt` - Primate species names for cactus447 track
- `primates_multiz.txt` - Primate species names for multiz470 track

## Pipeline Workflow

### Step 1: Extract from cactus447 Track

```bash
python fa_extractor_cactus447.py
```

Extracts alignments from the cactus447 database for all regions in CREs.bed.

### Step 2: Clean Sequences

```bash
python cleaning_fa.py
```

Removes gap characters and filters sequences with less than 70% coverage.

### Step 3: Trim Alignments

```bash
bash clipkit.sh
```

Removes uninformative columns using kpic-gappy mode.

### Step 4: Count Sequences (cactus447)

```bash
python CountSeqs.py
```

Counts species in each alignment and generates a BED file with regions that have fewer than 61 species.

Before running, set in `CountSeqs.py`:

```python
THRESHOLD = 61
```

### Step 5: Extract from multiz470 Track

```bash
python fa_extractor_multiz470.py
```

Extracts alignments for regions with <61 species using the multiz470 database. Repeat Steps 2-3.

### Step 6: Count Sequences (multiz470)

```bash
python CountSeqs.py
```

Counts species and generates a BED file with regions that have fewer than 27 species.

Before running, update in `CountSeqs.py`:

```python
THRESHOLD = 27
```

### Step 7: Extract from ucsc30 Track

```bash
python fa_extractor_ucsc30.py
```

Extracts alignments for regions with <27 species using the ucsc30 database. Repeat Steps 2-3.

### Step 8: Select Best Alignments

```bash
python SelectBestMSA.py
```

Compares all three track versions for each region and selects the one with the most species. Copies best alignments to `Adaptify_MSAs` and generates an Excel report with coordinates, sequence counts, and track information.

## Project Structure

```
├── CREs.bed
├── primates.txt
├── primates_multiz.txt
├── fa_extractor_cactus447.py
├── fa_extractor_multiz470.py
├── fa_extractor_ucsc30.py
├── cleaning_fa.py
├── clipkit.sh
├── CountSeqs.py
├── SelectBestMSA.py
└── Adaptify_MSAs/
```

## Configuration

- **msa_view path:** Update in all three extraction scripts with your path from `which msa_view`
- **Coverage threshold:** In `cleaning_fa.py`, default is 0.7 (70%)
- **Sequence thresholds:** In `CountSeqs.py`, set `THRESHOLD = 61` for Step 4 and `THRESHOLD = 27` for Step 6

