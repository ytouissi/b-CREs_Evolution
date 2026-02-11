# Youness-Toussi Research Projects

Repository containing automated bioinformatics pipelines for genomic analysis.

## Projects Overview

### 1. CREs MSAs (Multiple Sequence Alignments)

Automated workflow to extract and optimize sequence alignments for cis-regulatory elements across primate species.

**What it does:**
- Extracts alignments from three different genomic databases (cactus447, multiz470, ucsc30)
- Intelligently selects the best alignment for each region based on sequence coverage
- Cleans, trims, and optimizes alignments for phylogenetic analysis

**Quick start:**
```bash
cd 1-b-CREs_MSAs
bash run_pipeline.sh
```

See [1-b-CREs_MSAs/README.md](1-b-CREs_MSAs/README.md) for complete instructions.

---

### 2. Neutral MSAs

Analysis of neutral sequences for evolutionary comparison.

---

### 3. Positive Selection Test

Identifies genomic regions under positive selection across species.

---

### 4. Differential Expression

Analysis of gene expression differences between conditions.

---

## Getting Started

Each project folder contains its own README with specific instructions. Navigate to the project folder you want to work with and follow the setup steps.

### General Requirements

- Python 3 or higher
- conda package manager
- Unix-like system (macOS, Linux)

### Common Dependencies

Most projects require:
```bash
pip install biopython openpyxl
conda install -c bioconda phast clipkit
```

## Repository Structure

```
youness-toussi/
├── 1-b-CREs_MSAs/
│   ├── README.md
│   ├── WORKFLOW.md
│   ├── config.json
│   ├── run_pipeline.sh
│   └── [scripts and data files]
├── 2-Neutral_MSAs/
├── 3-Positive_Selection_Test/
├── 4-Differential_Expression/
└── README.md (this file)
```

## License

[Add your license here]

## Contact

For questions or issues, contact Youness Toussi

