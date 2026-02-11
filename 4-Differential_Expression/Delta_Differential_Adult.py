#!/usr/bin/env python3
"""
shared HAR Delta Expression Analysis - Individual-level version
Keeps individuals as replicates, computes gene-level log2FCs,
and tests HAR vs control distributions.
Generates both PNG and PDF outputs.
"""

import pandas as pd
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')


def load_data():
    """Load expression data and gene lists"""
    print("Loading expression data...")
    expr_data = pd.read_csv('adult_data_filtered.csv', index_col=0)
    print(f"Expression data shape: {expr_data.shape}")

    # Extract metadata from sample names
    samples = expr_data.columns.tolist()
    metadata_rows = []
    for sample in samples:
        parts = sample.split('.')
        individual = parts[0] if len(parts) > 0 else sample
        brain_region = parts[1] if len(parts) > 1 else 'Unknown'
        sample_str = str(sample).upper()
        if sample_str.startswith('HSB'):
            species = 'Human'
        elif sample_str.startswith('PTB'):
            species = 'Chimpanzee'
        elif sample_str.startswith('RMB'):
            species = 'Macaque'
        else:
            species = 'Unknown'
        metadata_rows.append({
            'Sample': sample,
            'individual': individual,
            'brain_region': brain_region,
            'Species': species
        })
    metadata = pd.DataFrame(metadata_rows)

    # Load gene lists
    def load_genes(filename):
        genes = set()
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    genes.add(line)
        return genes

    har_genes = load_genes('shared_HPS-bCREs.txt')
    control_genes = load_genes('shared_non_HPS-bCREs.txt')

    return metadata, expr_data, har_genes, control_genes


def prepare_individual_means(metadata, expr_data):
    """Average across brain regions for each individual"""
    print("Averaging expression across brain regions per individual...")
    individual_profiles = {}
    for (individual, species), group in metadata.groupby(['individual', 'Species']):
        samples = group['Sample'].tolist()
        if len(samples) > 0:
            avg_expr = expr_data[samples].mean(axis=1)
            individual_profiles[(individual, species)] = avg_expr
    return individual_profiles


def compute_log2fc(individual_profiles, species1, species2):
    """Compute per-gene log2FC between species1 and species2"""
    species1_expr = pd.DataFrame([v for (ind, sp), v in individual_profiles.items() if sp == species1])
    species2_expr = pd.DataFrame([v for (ind, sp), v in individual_profiles.items() if sp == species2])

    common_genes = set(species1_expr.columns) & set(species2_expr.columns)
    species1_expr = species1_expr[list(common_genes)]
    species2_expr = species2_expr[list(common_genes)]

    log2fc = np.log2(species1_expr.mean(axis=0) + 1) - np.log2(species2_expr.mean(axis=0) + 1)
    return log2fc


def analyze_comparison(individual_profiles, har_genes, control_genes, species1, species2):
    """Compare HAR vs control genes in log2FC distributions"""
    log2fc = compute_log2fc(individual_profiles, species1, species2)

    available_genes = set(log2fc.index)
    har_available = har_genes & available_genes
    control_available = control_genes & available_genes

    har_vals = log2fc[list(har_available)]
    control_vals = log2fc[list(control_available)]

    # Wilcoxon rank-sum test
    stat, p_val = ranksums(har_vals, control_vals)

    result = {
        'comparison': f"{species1}_vs_{species2}",
        'har_mean': har_vals.mean(),
        'control_mean': control_vals.mean(),
        'difference': har_vals.mean() - control_vals.mean(),
        'p_value': p_val,
        'har_count': len(har_vals),
        'control_count': len(control_vals),
        'har_vals': har_vals,
        'control_vals': control_vals
    }
    return result


def plot_all_comparisons(results_list, pvals_adj=None):
    """Plot all comparisons side by side with adjusted p-values"""
    if not results_list:
        print("No results to plot")
        return

    colors = {'HAR': '#4CCBA1', 'Control': '#ee6c4d'}
    n = len(results_list)
    cols = min(2, n)
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 5*rows), squeeze=False)
    axes = axes.flatten()

    for i, res in enumerate(results_list):
        ax = axes[i]
        har_vals = res['har_vals']
        control_vals = res['control_vals']

        # Use increased bandwidth to make curves wider
        sns.kdeplot(har_vals, fill=True, color=colors['HAR'], alpha=0.4, ax=ax, bw_adjust=1.5)
        sns.kdeplot(control_vals, fill=True, color=colors['Control'], alpha=0.4, ax=ax, bw_adjust=1.5)

        # Get current x-axis limits and expand them by 40% to make curves flatter
        xlim = ax.get_xlim()
        x_range = xlim[1] - xlim[0]
        new_xlim = (xlim[0] - 0.001 * x_range, xlim[1] + 0.001 * x_range)
        ax.set_xlim(new_xlim)

        ax.set_xlabel('log2FC Gene Expression ({} - {})'.format(res['comparison'].split('_vs_')[0], res['comparison'].split('_vs_')[1]), fontsize=12, fontweight='bold', fontfamily='Trebuchet MS')
        ax.set_ylabel('Density', fontsize=12, fontweight='bold', fontfamily='Trebuchet MS')
        
        # Add adjusted p-value at the bottom, closer to x-axis
        if pvals_adj is not None:
            # Get y-axis limit to position text near the bottom, closer to x-axis
            ylim = ax.get_ylim()
            y_pos = ylim[1] * 0.02  # Position at 8% of max y (closer to x-axis)
            x_pos = sum(new_xlim) / 2  # Middle of x-axis
            p_val = pvals_adj[i]
            exponent = int(np.floor(np.log10(abs(p_val))))
            mantissa = p_val / (10 ** exponent)
            p_text = f"Adj \n P-value =\n{mantissa:.2f}×10$^{{{exponent}}}$"  # Added \n for line break

            ax.text(0, y_pos, p_text,  # Changed x_pos to 0 to center above x=0
                fontsize=11, fontweight='bold', fontfamily='Trebuchet MS',
                ha='center', va='bottom')    
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Hide unused subplots
    for j in range(i+1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle('Differential gene expression associated with positively selected shared (adult) b-CREs in human brain.', fontsize=16, fontweight='bold', fontfamily='Trebuchet MS')
        # Create proxy artists for a single, clean figure legend (with subtle outlines)
    har_patch = Rectangle((0, 0), 1, 1, fc=colors['HAR'], alpha=0.5, ec=colors['HAR'], lw=1.0)
    control_patch = Rectangle((0, 0), 1, 1, fc=colors['Control'], alpha=0.5, ec=colors['Control'], lw=1.0)
    
    fig.legend(handles=[har_patch, control_patch], 
               labels=['Genes regulated by HS b-CREs.', 'Genes regulated by Non-HS b-CREs.'],
               loc='upper right', 
               bbox_to_anchor=(0.95, 0.95),
               fontsize=14,
               prop={'family': 'Trebuchet MS'},
               frameon=False)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("shared_HAR_comparisons.png", dpi=300)
    plt.savefig("shared_HAR_comparisons.pdf", format='pdf', bbox_inches='tight')
    plt.show()


def main():
    metadata, expr_data, har_genes, control_genes = load_data()
    individual_profiles = prepare_individual_means(metadata, expr_data)

    comparisons = []
    if any(sp == 'Human' for (_, sp) in individual_profiles) and any(sp == 'Macaque' for (_, sp) in individual_profiles):
        comparisons.append(('Human', 'Macaque'))
    if any(sp == 'Human' for (_, sp) in individual_profiles) and any(sp == 'Chimpanzee' for (_, sp) in individual_profiles):
        comparisons.append(('Human', 'Chimpanzee'))
    if any(sp == 'Chimpanzee' for (_, sp) in individual_profiles) and any(sp == 'Macaque' for (_, sp) in individual_profiles):
        comparisons.append(('Chimpanzee', 'Macaque'))

    results_list = []
    for sp1, sp2 in comparisons:
        res = analyze_comparison(individual_profiles, har_genes, control_genes, sp1, sp2)
        results_list.append(res)
        print(f"{res['comparison']}: HAR mean={res['har_mean']:.3f}, "
              f"Control mean={res['control_mean']:.3f}, "
              f"Δ={res['difference']:.3f}, p={res['p_value']:.2e}")

    pvals_adj = None
    if len(results_list) > 1:
        pvals = [r['p_value'] for r in results_list]
        reject, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')
        for i, r in enumerate(results_list):
            print(f"{r['comparison']} adjusted p={pvals_adj[i]:.2e}, HPS-bCREs={reject[i]}")

    # Generate all outputs with adjusted p-values
    plot_all_comparisons(results_list, pvals_adj)


if __name__ == "__main__":
    main()