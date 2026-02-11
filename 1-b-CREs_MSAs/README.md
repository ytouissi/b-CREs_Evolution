# Multiple Sequence Alignment

Automated workflow to extract and optimize multiple sequence alignments (MSAs) for cis-regulatory elements (CREs) across primate species.

## Before You Start

Download all files from this repository and save them to a folder named `Multiple_Sequence_Alignment` on your desktop.

## Installation

### Requirements

- Python 3 or higher
- conda package manager

### Install PHAST

```bash
conda install -c bioconda phast
```

Verify installation:

```bash
which msa_view
```

Copy the output path - you will need this for configuration.

### Install Python Libraries

```bash
pip install biopython openpyxl
```

### Install clipkit

```bash
conda install -c bioconda clipkit
```

## Get UCSC Browser API Key

1. [Sign up for UCSC account](https://genome.ucsc.edu/cgi-bin/hgLogin?hgLogin.do.signupPage=1) or [log in](https://genome.ucsc.edu/cgi-bin/hgLogin?hgLogin.do.displayLoginPage=1)
2. Go to [Track Data Hubs](https://genome.ucsc.edu/cgi-bin/hgHubConnect#hubDeveloper)
3. Click on "Hub Development" section at the bottom
4. Your API key is displayed there
5. Copy this key for the next step

## Configuration

### Edit config.json

Open `config.json` in a text editor and update it:

```json
{
  "msa_view_path": "/path/from/which/msa_view",
  "api_key": "your_ucsc_api_key_here"
}
```

Replace:
- `/path/from/which/msa_view` with the result from `which msa_view`
- `your_ucsc_api_key_here` with your UCSC API key from Hub Development

## Run the Pipeline

```bash
bash Get_MSAs.sh
```

The script runs all extraction, cleaning, trimming, and analysis steps automatically. When complete:

- Final results are in `Adaptify_MSAs/` with a report saved as `MSA_Selection_Report.xlsx`
- Intermediate files (raw alignments and individual track final alignments) are in intermediate output folders

## What the Pipeline Does

The pipeline automates the following workflow:

1. Extracts alignments from cactus447 database
2. Cleans sequences (removes gaps, filters low coverage)
3. Trims uninformative columns
4. Counts sequences in each alignment
5. For regions with <61 species, extracts from multiz470 database
6. Cleans and trims multiz470 data
7. Counts sequences again
8. For regions with <27 species, extracts from ucsc30 database
9. Cleans and trims ucsc30 data
10. Selects the best version of each region (most species)
11. Generates a report

For detailed explanations of each step, see [WORKFLOW.md](Workflow.md).
