# Term-Specific Permutation Test for ORA Enrichment

## Overview

This script performs a term-specific permutation sanity check for Over-Representation Analysis (ORA) using g:Profiler. It evaluates whether enrichment terms found in genes regulated by Human Selected Regions (b-CREs) are statistically robust or likely to appear by chance.

**Table S7** contains the enrichment terms enriched in genes regulated by Human Selected Regions (b-CREs) across three stages: Adult, Fetal, and Shared. The script uses these terms as targets and tests how often they appear in random permutations of the same size.

---

## Requirements

### R Packages

- `gprofiler2`
- `readxl`
- `ggplot2`
- `parallel`

Install missing packages with:

```r
install.packages(c("gprofiler2", "readxl", "ggplot2"))
```

---

## Input Files

All files must be placed in the **same directory as the script**:

| File | Description |
|---|---|
| `Enrichment_terms_permutation.R` | The main R script |
| `Table_S7.xlsx` | Enrichment results for HPS b-CREs (Adult, Fetal, Shared) |
| `adult_HPS-bCREs.txt` | Foreground: genes regulated by Adult HPS b-CREs |
| `adult_non_HPS-bCREs.txt` | Background genes regulated by Adult non-HPS b-CREs |
| `fetal_HPS-bCREs.txt` | Foreground genes regulated by Fetal HPS b-CREs |
| `fetal_non_HPS-bCREs.txt` | Background genes regulated by Fetal non-HPS b-CREs |
| `shared_HPS-bCREs.txt` | Foreground genes regulated by Shared HPS b-CREs |
| `shared_non_HPS_b-CREs.txt` | Background genes regulated by Shared non-HPS b-CREs |


---

## How to Run

Place the script and all input files listed above in the same directory.

---

## What the Script Does

For each stage (Adult, Fetal, Shared):

1. Reads the significant enrichment terms from **Table S7** (adjusted p-value < 0.05).
2. Runs **1,000 permutations**  each time randomly sampling foreground genes from the background and running a full g:Profiler ORA.
3. Counts how many times each real enrichment term appears across the 1,000 random draws.
4. Computes an **empirical p-value** per term: how often a random draw reproduces that term.


---

## Output Files

| File | Description |
|---|---|
| `term_specific_permutation_results.csv` | Full results table with real adjusted p-value, number of times seen in permutations, and empirical p-value |
---
