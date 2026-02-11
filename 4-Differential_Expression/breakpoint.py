#!/usr/bin/env python3
"""
HAR Gene Breakpoint Analysis - Dual Graph Comparison
Left: Fetal b-CREs
Right: Shared (Fetal) b-CREs
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, gaussian_kde, mannwhitneyu, shapiro, normaltest, kstest
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

def load_data():
    """Load breakpoint data and gene lists"""
    print("Loading breakpoint data from Excel file...")
    
    try:
        # Load the Bakken breakpoint data
        breakpoint_data = pd.read_excel('Genes breakpoint.xlsx', sheet_name='Sheet1')
        print(f"Loaded {len(breakpoint_data)} rows of breakpoint data")
        
        # Check data structure
        print(f"Columns: {list(breakpoint_data.columns)}")
        print(f"Species: {breakpoint_data['species'].unique()}")
        
        regions = breakpoint_data['region'].unique()
        print(f"Brain regions ({len(regions)}): {sorted(regions)}")
        
        # Show data distribution by region and species
        print("\nData distribution by region:")
        for region in sorted(regions):
            region_data = breakpoint_data[breakpoint_data['region'] == region]
            species_counts = region_data['species'].value_counts()
            human_count = species_counts.get('human', 0)
            macaque_count = species_counts.get('macaque', 0)
            rat_count = species_counts.get('rat', 0)
            print(f"  {region}: Human={human_count}, Macaque={macaque_count}, Rat={rat_count}")
        
    except FileNotFoundError:
        print("\nERROR: Excel file 'Genes breakpoint.xlsx' not found.")
        print("Please make sure all data files are in the same folder as the script.")
        return None
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        return None
    
    # Load gene lists for FETAL
    print(f"\nLoading FETAL gene lists...")
    def load_genes(filename):
        try:
            with open(filename, 'r') as f:
                genes = set()
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        genes.add(line)
                return genes
        except FileNotFoundError:
            print(f"\nERROR: Gene list file '{filename}' not found.")
            print("Please make sure all data files are in the same folder as the script.")
            return None

    fetal_har_genes = load_genes('fetal_HPS-bCREs.txt')
    fetal_control_genes = load_genes('fetal_non_HPS-bCREs.txt')

    if fetal_har_genes is None or fetal_control_genes is None:
        return None, None, None, None, None
    
    print(f"Fetal HAR genes: {len(fetal_har_genes)}")
    print(f"Fetal Control genes: {len(fetal_control_genes)}")
    
    # Load gene lists for SHARED FETAL
    print(f"\nLoading SHARED FETAL gene lists...")
    shared_har_genes = load_genes('shared_HPS-bCREs.txt')
    shared_control_genes = load_genes('shared_non_HPS-bCREs.txt')

    if shared_har_genes is None or shared_control_genes is None:
        return None, None, None, None, None
    
    print(f"Shared Fetal HAR genes: {len(shared_har_genes)}")
    print(f"Shared Fetal Control genes: {len(shared_control_genes)}")
    
    return breakpoint_data, fetal_har_genes, fetal_control_genes, shared_har_genes, shared_control_genes

def calculate_delta_breakpoints(breakpoint_data):
    """Calculate delta breakpoints (human - macaque) for each gene-region pair"""
    print("\nCalculating h_m_escore (human - macaque) exactly like paper's R code...")
    
    # Filter for human and macaque data only
    human_data = breakpoint_data[breakpoint_data['species'] == 'human'].copy()
    macaque_data = breakpoint_data[breakpoint_data['species'] == 'macaque'].copy()
    
    print(f"Human data points: {len(human_data)}")
    print(f"Macaque data points: {len(macaque_data)}")
    
    # Create gene-region keys for matching
    human_data['gene_region'] = human_data['gene'] + '_' + human_data['region']
    macaque_data['gene_region'] = macaque_data['gene'] + '_' + macaque_data['region']
    
    # Merge human and macaque data
    merged_data = pd.merge(
        human_data[['gene', 'region', 'gene_region', 'bp (event score)']],
        macaque_data[['gene_region', 'bp (event score)']],
        on='gene_region',
        suffixes=('_human', '_macaque')
    )
    
    print(f"Gene-region pairs with both human and macaque data: {len(merged_data)}")
    
    # Calculate h_m_escore (human - macaque) exactly like the paper
    merged_data['h_m_escore'] = merged_data['bp (event score)_human'] - merged_data['bp (event score)_macaque']
    
    print(f"h_m_escore range: {merged_data['h_m_escore'].min():.3f} to {merged_data['h_m_escore'].max():.3f}")
    print(f"Mean h_m_escore: {merged_data['h_m_escore'].mean():.3f}")
    
    # Show some examples
    print("\nExample h_m_escore calculations:")
    for _, row in merged_data.head(5).iterrows():
        print(f"{row['gene']} ({row['region']}): Human={row['bp (event score)_human']:.3f}, "
              f"Macaque={row['bp (event score)_macaque']:.3f}, h_m_escore={row['h_m_escore']:.3f}")
    
    return merged_data

def test_normality(data, name):
    """Test normality of data distribution"""
    print(f"\n--- Testing normality for {name} (n={len(data)}) ---")
    
    # Remove NaN values
    clean_data = data.dropna()
    
    if len(clean_data) < 8:
        print(f"Too few data points ({len(clean_data)}) for normality testing")
        return False
    
    # Shapiro-Wilk test (best for small samples)
    if len(clean_data) <= 5000:
        shapiro_stat, shapiro_p = shapiro(clean_data)
        print(f"Shapiro-Wilk test: W = {shapiro_stat:.4f}, p = {shapiro_p:.2e}")
        shapiro_normal = shapiro_p > 0.05
    else:
        shapiro_normal = None
        print("Shapiro-Wilk test: Skipped (sample too large)")
    
    # D'Agostino's normality test
    if len(clean_data) >= 20:
        dagostino_stat, dagostino_p = normaltest(clean_data)
        print(f"D'Agostino test: stat = {dagostino_stat:.4f}, p = {dagostino_p:.2e}")
        dagostino_normal = dagostino_p > 0.05
    else:
        dagostino_normal = None
        print("D'Agostino test: Skipped (sample too small)")
    
    # Kolmogorov-Smirnov test against normal distribution
    # Standardize the data first
    standardized = (clean_data - clean_data.mean()) / clean_data.std()
    ks_stat, ks_p = kstest(standardized, 'norm')
    print(f"KS test vs normal: D = {ks_stat:.4f}, p = {ks_p:.2e}")
    ks_normal = ks_p > 0.05
    
    # Summary
    tests_passed = []
    if shapiro_normal is not None:
        tests_passed.append(shapiro_normal)
    if dagostino_normal is not None:
        tests_passed.append(dagostino_normal)
    tests_passed.append(ks_normal)
    
    normal_count = sum(tests_passed)
    total_tests = len(tests_passed)
    
    print(f"Normality summary: {normal_count}/{total_tests} tests suggest normal distribution")
    
    # Decision: majority rule
    is_normal = normal_count > total_tests / 2
    print(f"Conclusion: {'NORMAL' if is_normal else 'NON-NORMAL'} distribution")
    
    return is_normal

def analyze_breakpoints_like_paper(merged_data, har_genes, control_genes, dataset_name):
    """
    Analyze breakpoints exactly like the paper's R code:
    - Calculate h_m_escore = human - macaque for each gene-region pair
    - Test data distribution to choose appropriate statistical test
    - Compare HAR vs non-HAR distributions
    """
    print(f"\n--- ANALYZING {dataset_name} ---")
    print("Following their R code: wilcox.test(h_m_escore~har, data=bp.diff)")
    
    # Create the equivalent of their bp.diff dataframe
    bp_diff = merged_data.copy()
    
    print(f"Total gene-region pairs: {len(bp_diff)}")
    print(f"h_m_escore range: {bp_diff['h_m_escore'].min():.3f} to {bp_diff['h_m_escore'].max():.3f}")
    
    # Classify genes as HAR vs non-HAR (like their code)
    bp_diff['har_status'] = bp_diff['gene'].apply(lambda x: 'HAR' if x in har_genes else 'nHAR')
    
    # Get h_m_escore values for each group
    har_scores = bp_diff[bp_diff['har_status'] == 'HAR']['h_m_escore']
    nhar_scores = bp_diff[bp_diff['har_status'] == 'nHAR']['h_m_escore']
    
    print(f"\nHAR gene-region pairs: {len(har_scores)}")
    print(f"Non-HAR gene-region pairs: {len(nhar_scores)}")
    
    if len(har_scores) < 5 or len(nhar_scores) < 5:
        print("Insufficient data for analysis")
        return None
    
    # Remove NaN values
    har_scores = har_scores.dropna()
    nhar_scores = nhar_scores.dropna()
    
    print(f"After removing NaN - HAR: {len(har_scores)}, Non-HAR: {len(nhar_scores)}")
    
    # Test normality of both distributions
    har_normal = test_normality(har_scores, f"{dataset_name} HAR genes")
    nhar_normal = test_normality(nhar_scores, f"{dataset_name} Non-HAR genes")
    
    # Decide on statistical test
    print(f"\n--- CHOOSING STATISTICAL TEST ---")
    both_normal = har_normal and nhar_normal
    
    if both_normal:
        print("Both distributions appear normal → Using t-test")
        print("However, paper used Wilcoxon test → Will run both for comparison")
        primary_test = "parametric"
    else:
        print("At least one distribution is non-normal → Using Wilcoxon test")
        print("This matches the paper's choice of non-parametric test")
        primary_test = "non_parametric"
    
    # Run both tests for completeness
    # Wilcoxon rank-sum test (equivalent to their wilcox.test)
    wilcox_stat, wilcox_p = mannwhitneyu(har_scores, nhar_scores, alternative='two-sided')
    
    # t-test for comparison
    t_stat, t_p = ttest_ind(har_scores, nhar_scores, equal_var=False)
    
    result = {
        'har_count': len(har_scores),
        'nhar_count': len(nhar_scores),
        'har_mean': har_scores.mean(),
        'har_std': har_scores.std(),
        'har_median': har_scores.median(),
        'nhar_mean': nhar_scores.mean(),
        'nhar_std': nhar_scores.std(),
        'nhar_median': nhar_scores.median(),
        'mean_difference': har_scores.mean() - nhar_scores.mean(),
        'median_difference': har_scores.median() - nhar_scores.median(),
        'wilcox_stat': wilcox_stat,
        'wilcox_p': wilcox_p,
        't_stat': t_stat,
        't_p': t_p,
        'har_values': har_scores,
        'nhar_values': nhar_scores,
        'primary_test': primary_test,
        'har_normal': har_normal,
        'nhar_normal': nhar_normal
    }
    
    print(f"\n--- RESULTS ---")
    print(f"HAR gene-regions: n={len(har_scores)}")
    print(f"  Mean h_m_escore = {har_scores.mean():.4f} ± {har_scores.std():.4f}")
    print(f"  Median h_m_escore = {har_scores.median():.4f}")
    
    print(f"Non-HAR gene-regions: n={len(nhar_scores)}")
    print(f"  Mean h_m_escore = {nhar_scores.mean():.4f} ± {nhar_scores.std():.4f}")
    print(f"  Median h_m_escore = {nhar_scores.median():.4f}")
    
    print(f"\nDifferences:")
    print(f"  Mean difference: {result['mean_difference']:.4f}")
    print(f"  Median difference: {result['median_difference']:.4f}")
    
    print(f"\nStatistical tests:")
    print(f"  Wilcoxon test: W = {wilcox_stat:.0f}, p = {wilcox_p:.2e}")
    print(f"  t-test: t = {t_stat:.3f}, p = {t_p:.2e}")
    
    # Primary result based on data distribution
    if primary_test == "non_parametric":
        primary_p = wilcox_p
        print(f"\nPrimary result (Wilcoxon): p = {primary_p:.2e}")
    else:
        primary_p = t_p
        print(f"\nPrimary result (t-test): p = {primary_p:.2e}")
    
    result['primary_p'] = primary_p
    
    # Interpretation matching the paper
    if result['mean_difference'] < 0:
        print(f"\n✓ HAR genes have MORE NEGATIVE h_m_escore (earlier breakpoints in humans)")
        print(f"  This SUPPORTS the paper's conclusion about early HAR gene expression")
    else:
        print(f"\n✗ HAR genes have LESS NEGATIVE h_m_escore (later breakpoints in humans)")
        print(f"  This would CONTRADICT the paper's conclusion")
    
    return result

def create_dual_distribution_curves(result_fetal, result_shared):
    """Create side-by-side distribution curves for Fetal and Shared datasets"""
    if not result_fetal or not result_shared:
        print("Missing results for visualization")
        return
    
    print("Creating dual distribution curves...")
    
    # Set font properties globally
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Trebuchet MS']
    plt.rcParams['font.weight'] = 'normal'
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['figure.titleweight'] = 'bold'
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['figure.titlesize'] = 20
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    
    # Color palette
    colors = {
        'har': '#0077b6',       # Clean, rich blue
        'nhar': '#ee6c4d',      # Clean, vibrant red-orange
        'background': '#F8F9FA'
    }
    
    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 6))
    fig.patch.set_facecolor('white')
    fig.suptitle('Expression breakpoint analysis for genes regulated by positively selected b-CREs.', 
            fontsize=22, fontweight='bold', y=1.02)
    
    # Create proxy artists for a single, shared figure legend
    har_patch = Rectangle((0, 0), 1, 1, fc=colors['har'], alpha=0.5, ec=colors['har'], lw=1.0)
    control_patch = Rectangle((0, 0), 1, 1, fc=colors['nhar'], alpha=0.5, ec=colors['nhar'], lw=1.0)
    
    # Add legend once at the figure level
    fig.legend(handles=[har_patch, control_patch], 
           labels=['Genes regulated by HS b-CREs.', 'Genes regulated by Non-HS b-CREs.'],
           loc='upper right', 
           bbox_to_anchor=(0.98, 0.98),
           prop={'family': 'Trebuchet MS', 'size': 12},
           frameon=False)
    
    # Plot LEFT graph - Fetal b-CREs
    plot_single_distribution(ax1, result_fetal, colors, "Fetal b-CREs")
    
    # Plot RIGHT graph - Shared (Fetal) b-CREs
    plot_single_distribution(ax2, result_shared, colors, "Shared (Fetal) b-CREs")
    
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    
    # Save with high quality
    plt.savefig('har_breakpoint_dual_distribution.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.savefig('har_breakpoint_dual_distribution.pdf', bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.show()
    
    print("Distribution curves saved as:")
    print("- har_breakpoint_dual_distribution.png")
    print("- har_breakpoint_dual_distribution.pdf")

def plot_single_distribution(ax, result, colors, title_text):
    """Plot a single distribution curve on given axis"""
    har_values = result['har_values']
    nhar_values = result['nhar_values']
    
    # Calculate range for smooth curves
    x_min = min(har_values.min(), nhar_values.min()) - 0.3
    x_max = max(har_values.max(), nhar_values.max()) + 0.3
    x_range = np.linspace(x_min, x_max, 500)
    
    # Plot Non-HAR Curve (Red-Orange) - BASE LAYER
    if len(nhar_values) > 1:
        try:
            kde_nhar = gaussian_kde(nhar_values)
            nhar_density = kde_nhar(x_range)
            # Semi-transparent fill
            ax.fill_between(x_range, nhar_density, color=colors['nhar'], alpha=0.4, zorder=2)
            # Subtle line outline
            ax.plot(x_range, nhar_density, color=colors['nhar'], linewidth=1.0, alpha=0.4, zorder=4)
        except Exception:
            ax.hist(nhar_values, bins=30, alpha=0.5, color=colors['nhar'], 
                   density=True, histtype='stepfilled', edgecolor=colors['nhar'], linewidth=1.0)
    
    # Plot HAR Curve (Blue) - OVERLAY LAYER
    if len(har_values) > 1:
        try:
            kde_har = gaussian_kde(har_values)
            har_density = kde_har(x_range)
            # Semi-transparent fill
            ax.fill_between(x_range, har_density, color=colors['har'], alpha=0.4, zorder=3)
            # Subtle line outline
            ax.plot(x_range, har_density, color=colors['har'], linewidth=1.0, alpha=0.4, zorder=5)
        except Exception:
            ax.hist(har_values, bins=30, alpha=0.5, color=colors['har'], 
                   density=True, histtype='stepfilled', edgecolor=colors['har'], linewidth=1.0)
    
    # Styling
    ax.set_xlabel('Δ Developemental Event Score (Human - Macaque)', fontweight='bold', fontsize=16)
    ax.set_ylabel('Density', fontweight='bold', fontsize=16)
    
    # Clean spines - remove top and right
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['bottom'].set_color('black')
    
    # Set normal font weight for tick labels
    ax.tick_params(axis='both', which='major', labelsize=12)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('normal')
    
    # Add title text in top left corner
    ax.text(0.02, 0.98, title_text,
        transform=ax.transAxes,
        fontsize=18, fontweight='bold', fontfamily='Trebuchet MS',
        ha='left', va='top')
    
    # Add p-value annotation at the center bottom
    p_val = result['primary_p']
    exponent = int(np.floor(np.log10(abs(p_val))))
    mantissa = p_val / (10 ** exponent)
    p_text = f"Adj P-value =\n{mantissa:.2f}×10$^{{{exponent}}}$"
    
    # Get y-axis limits to position text
    ylim = ax.get_ylim()
    y_pos = ylim[1] * 0.08  # Position at 8% of max y
    
    ax.text(0, y_pos, p_text,
        fontsize=15, fontweight='bold', fontfamily='Trebuchet MS',
        ha='center', va='bottom')

def main():
    """Main analysis function"""
    print("HAR GENE BREAKPOINT ANALYSIS - DUAL COMPARISON")
    print("=" * 70)
    print("Comparing Fetal b-CREs vs Shared (Fetal) b-CREs")
    print("=" * 70)
    
    # Load data
    data = load_data()
    
    if data[0] is None:
        print("\nScript terminated: Could not load required data files.")
        return
    
    breakpoint_data, fetal_har_genes, fetal_control_genes, shared_har_genes, shared_control_genes = data
    
    # Calculate h_m_escore
    merged_data = calculate_delta_breakpoints(breakpoint_data)
    
    # Analyze FETAL dataset
    print("\n" + "="*70)
    print("ANALYSIS 1: FETAL b-CREs")
    print("="*70)
    result_fetal = analyze_breakpoints_like_paper(merged_data, fetal_har_genes, fetal_control_genes, "FETAL")
    
    # Analyze SHARED FETAL dataset
    print("\n" + "="*70)
    print("ANALYSIS 2: SHARED (FETAL) b-CREs")
    print("="*70)
    result_shared = analyze_breakpoints_like_paper(merged_data, shared_har_genes, shared_control_genes, "SHARED FETAL")
    
    if result_fetal and result_shared:
        # Create dual distribution curves
        create_dual_distribution_curves(result_fetal, result_shared)
        
        # Final summary
        print(f"\n" + "="*70)
        print("FINAL SUMMARY - BOTH ANALYSES")
        print("="*70)
        
        print("\nFETAL b-CREs:")
        print(f"  Data: {result_fetal['har_count']} HAR gene-regions, {result_fetal['nhar_count']} non-HAR gene-regions")
        print(f"  p-value: {result_fetal['primary_p']:.2e}")
        significance = "HPS-bCREs" if result_fetal['primary_p'] < 0.05 else "NOT HPS-bCREs"
        print(f"  Result: {significance}")
        
        print("\nSHARED (FETAL) b-CREs:")
        print(f"  Data: {result_shared['har_count']} HAR gene-regions, {result_shared['nhar_count']} non-HAR gene-regions")
        print(f"  p-value: {result_shared['primary_p']:.2e}")
        significance = "HPS-bCREs" if result_shared['primary_p'] < 0.05 else "NOT HPS-bCREs"
        print(f"  Result: {significance}")
    
    print(f"\nAnalysis complete!")

if __name__ == "__main__":
    main()