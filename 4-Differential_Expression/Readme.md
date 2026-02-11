# Functional Analysis of Positively Selected b-CREs

Comprehensive analysis of genes regulated by human-positively selected brain cis-regulatory elements (HPS b-CREs) across gene expression, developmental timing, and co-expression networks.

## Overview

Following branch-specific positive selection detection (HyPhy), these analyses characterize functional consequences of HPS b-CREs through comparative transcriptomic analysis across human, chimpanzee, and macaque spanning developmental and adult stages.

## Input Files Required

Excel/Data files (same directory as scripts):
- `aat8077_tabless1s3.xlsx` - Developmental metadata from Zhu et al. (2018) supplementary materials
- `adult_data_filtered.csv` - Adult stage expression data
- `dev_data_filtered.csv` - Developmental stage expression data
- `Genes breakpoint.xlsx` - Developmental breakpoint scores with metadata (Bakken et al. 2016)
- `Gene Modules.xlsx` - WGCNA co-expression module assignments with gene IDs and metadata (Sousa et al. 2017)

Gene lists (one gene ID per line):
- `shared_HPS-bCREs.txt`, `shared_non_HPS-bCREs.txt` - Genes regulated by Shared b-CREs
- `fetal_HPS-bCREs.txt`, `fetal_non_HPS-bCREs.txt` - Genes regulated by Fetal b-CREs
- `adult_all_HPS-bCREs.txt`, `adult_all_non_HPS-bCREs.txt` -Genes regulated by Adult b-CREs

## Analyses

### 1. Delta Expression (Developmental)
**Script:** `Delta_Expression.py`

Calculates Δ Z-score (human minus macaque/chimpanzee expression Z-scores) across developmental time points. Tests if HPS-regulated genes show different expression patterns than controls (Mann-Whitney U test, FDR 5%).

**Input:** `aat8077_tabless1s3.xlsx`, gene lists, `dev_data_filtered.csv`

**Output:** Density plots, statistical results

### 2. Delta Expression (Adult)
**Script:** `Delta_Differential_Adult.py`

Calculates log2FC for adult brains across species pairs. Compares HPS vs non-HPS gene expression distributions (Wilcoxon rank-sum test, FDR 5%).

**Input:** `aat8077_tabless1s3.xlsx`, gene lists, `adult_data_filtered.csv`

**Output:** Density plots, statistical results

### 3. Breakpoint Analysis
**Script:** `breakpoint.py`

Identifies developmental stages with abrupt expression changes. Compares Δ breakpoint scores (human minus macaque) between HPS and control genes (Wilcoxon rank-sum test, FDR 5%).

**Input:** `Genes breakpoint.xlsx`, gene lists

**Output:** Distribution plots, timing comparisons

### 4. WGCNA Enrichment
**Script:** `WGCNA.py`

Tests if HPS-regulated genes are enriched in specific co-expression modules using logistic linear regression. Identifies modules containing concentrated HPS regulatory signal (Bonferroni and FDR 5% correction).

**Input:** `Gene Modules.xlsx`, gene lists

**Output:** Module enrichment results (odds ratios, p-values) for shared/fetal/adult gene lists

## Data Sources

**PsychENCODE Human Brain Evolution Platform**
- http://evolution.psychencode.org/ (Processed data → mRNA-seq)

Expression data and metadata from supplementary materials of:
1. Sousa et al. (2017). Science 358(6366):1027-1032
2. Zhu et al. (2018). Science 362(6420)
3. Bakken et al. (2016). Nature 535(7612):367-375

## Running Analyses

```bash
python3 Delta_Expression.py
python3 Delta_Differential_Adult.py
python3 breakpoint.py
python3 WGCNA.py
```

## Statistical Methods

- Mann-Whitney U test (developmental expression)
- Wilcoxon rank-sum test (adult expression, breakpoint timing)
- Logistic linear regression (co-expression enrichment, WGCNA)
- Multiple testing: FDR 5% (Benjamini-Hochberg) and Bonferroni correction

## References

1. Sousa, A. M. M., Zhu, Y., Raghanti, M. A., Kitchen, R. R., Onorati, M., Tebbenkamp, A. T. N., Stutz, B., Meyer, K. A., Li, M., Kawasawa, Y. I., Liu, F., Perez, R. G., Mele, M., Carvalho, T., Skarica, M., Gulden, F. O., Pletikos, M., Shibata, A., Stephenson, A. R., & Sestan, N. (2017). Molecular and cellular reorganization of neural circuits in the human lineage. *Science*, 358(6366), 1027–1032. https://doi.org/10.1126/science.aan3456

2. Zhu, Y., Sousa, A. M. M., Gao, T., Skarica, M., Li, M., Santpere, G., Esteller-Cucala, P., Juan, D., Ferrández-Peral, L., Gulden, F. O., Yang, M., Miller, D. J., Marques-Bonet, T., Kawasawa, Y. I., Zhao, H., & Sestan, N. (2018). Spatiotemporal transcriptomic divergence across human and macaque brain development. *Science*, 362(6420). https://doi.org/10.1126/science.aat8077

3. Bakken, T. E., Miller, J. A., Ding, S., Sunkin, S. M., Smith, K. A., Ng, L., Szafer, A., Dalley, R. A., Royall, J. J., Lemon, T., Shapouri, S., Aiona, K., Arnold, J., Bennett, J. L., Bertagnolli, D., Bickley, K., Boe, A., Brouner, K., Butler, S., & Lein, E. S. (2016). A comprehensive transcriptional map of primate brain development. *Nature*, 535(7612), 367–375. https://doi.org/10.1038/nature18637
