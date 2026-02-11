# Adaptive Selection Analysis (HyPhy)

Detects branch-specific positive selection in b-CRE regions using HyPhy likelihood ratio tests.

## Overview

This script performs branch-specific positive selection testing on query b-CREs against neutral reference alignments. It tests for positive selection on human lineage and great ape clade branches using phylogenetic models.

## Input Data

- **Queries/** folder - Quality b-CRE MSA files from Step 1
- **References/** folder - Concatenated neutral MSAs from Step 3
- **243Primates_Tree.nwk** - Master phylogenetic tree (same directory)
- **null4-fgrnd_spec.bf** - HyPhy null model template (same directory)
- **alt4-fgrnd_spec.bf** - HyPhy alternative model template (same directory)

## Workflow

1. Prunes master tree to match query-reference species
2. Runs HyPhy null and alternative models for human and great ape branches
3. Calculates likelihood ratio tests (LRT) and p-values
4. Outputs results with significance calls

## Output

- **evolutionary_analysis_results.xlsx** - Results with LRT statistics, p-values, and significance for each branch
- **processed_bCREs.json** - Tracks completed analyses (allows script to resume if interrupted)
- **script_log.txt** - Analysis log

## Usage

```bash
python3 hyphy.py
```

## Requirements

- Python 3.x with BioPython, pandas, scipy, ete3
- HyPhy installed and accessible from PATH

## Next Step

Results ready for identifying positively selected regulatory regions in human and great ape lineages.
