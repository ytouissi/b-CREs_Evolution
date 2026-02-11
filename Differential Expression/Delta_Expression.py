#!/usr/bin/env python3
"""
Distribution Curves Only - fetal HAR Delta Expression Analysis
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, gaussian_kde
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

def load_data():
    """Load all required data with correct file handling"""
    print("Loading metadata from Excel file...")
    
    # Load Excel metadata - try different approaches to find the data
    try:
        # First try: load and examine the structure
        excel_raw = pd.read_excel('aat8077_tabless1s3.xlsx', sheet_name='Table S1', header=None)
        
        # Find the row containing "Sample" 
        header_row = None
        for i, row in excel_raw.iterrows():
            row_values = [str(cell) for cell in row.values if pd.notna(cell)]
            if any('Sample' in val for val in row_values):
                header_row = i
                break
        
        if header_row is None:
            # Try reading with different skiprows
            for skip in [0, 1, 2, 3]:
                try:
                    test_df = pd.read_excel('aat8077_tabless1s3.xlsx', sheet_name='Table S1', skiprows=skip)
                    if 'Sample' in test_df.columns or any('sample' in str(col).lower() for col in test_df.columns):
                        metadata = test_df
                        print(f"Found data by skipping {skip} rows")
                        break
                except:
                    continue
            else:
                print("Could not find Sample column in any configuration")
                return None, None, None, None, None
        else:
            # Read with the found header row
            metadata = pd.read_excel('aat8077_tabless1s3.xlsx', sheet_name='Table S1', header=header_row)
            print(f"Found header at row {header_row}")
        
    except FileNotFoundError:
        print("\nERROR: Excel file 'aat8077_tabless1s3.xlsx' not found.")
        print("Please make sure all data files are in the same folder as the script.")
        return None, None, None, None, None
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        return None, None, None, None, None
    metadata.columns = metadata.columns.str.strip()
    
    # Find predicted period column
    period_col = None
    for col in metadata.columns:
        if 'predicted period' in col.lower():
            period_col = col
            break
    
    if period_col is None:
        print("Could not find 'Predicted period' column")
        return None, None, None, None, None
    
    metadata = metadata.rename(columns={period_col: 'predicted_period'})
    
    # Remove rows with missing samples or periods
    metadata = metadata.dropna(subset=['Sample', 'predicted_period'])
    
    # Fix species labels based on sample names
    metadata['Species'] = metadata['Sample'].apply(
        lambda x: 'Human' if str(x).startswith('HSB') else 'Macaque' if str(x).startswith('RMB') else 'Unknown'
    )
    
    # Extract individual and region
    metadata['individual'] = metadata['Sample'].str.split('.').str[0]
    metadata['brain_region'] = metadata['Sample'].str.split('.').str[1]
    
    print(f"Loaded {len(metadata)} samples from Excel")
    print(f"Species distribution: {metadata['Species'].value_counts().to_dict()}")
    
    # Load expression data
    print(f"\nLoading expression data...")
    try:
        expr_data = pd.read_csv('dev_data_filtered.csv', index_col=0)
    except FileNotFoundError:
        print("\nERROR: Expression file 'dev_data_filtered.csv' not found.")
        print("Please make sure all data files are in the same folder as the script.")
        return None, None, None, None, None
    print(f"Expression data shape: {expr_data.shape}")
    
    # Check sample overlap
    metadata_samples = set(metadata['Sample'])
    expr_samples = set(expr_data.columns)
    common_samples = metadata_samples & expr_samples
    
    print(f"Samples in metadata: {len(metadata_samples)}")
    print(f"Samples in expression: {len(expr_samples)}")
    print(f"Common samples: {len(common_samples)}")
    
    if len(common_samples) == 0:
        print("ERROR: No common samples between metadata and expression data!")
        return None, None, None, None, None
    
    # Load gene lists
    print(f"\nLoading gene lists...")
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

    fetal_har_genes = load_genes('shared_fetal_significant.txt')
    control_genes = load_genes('shared_fetal_non_significant.txt')

    if fetal_har_genes is None or control_genes is None:
        return None, None, None, None, None
    
    print(f"fetal HAR genes: {len(fetal_har_genes)}")
    print(f"Control genes: {len(control_genes)}")
    
    # Create period to age mapping (using Human species only)
    period_age_map = {}
    if 'Age' in metadata.columns:
        human_data = metadata[metadata['Species'] == 'Human']
        for period in human_data['predicted_period'].unique():
            period_data = human_data[human_data['predicted_period'] == period]
            if len(period_data) > 0:
                # Get the most common age for this period
                age_value = period_data['Age'].mode()
                if len(age_value) > 0:
                    period_age_map[period] = str(age_value.iloc[0])
    
    print(f"Period to Age mapping (Human): {period_age_map}")
    
    return metadata, expr_data, fetal_har_genes, control_genes, period_age_map

def analyze_fetal_periods(metadata, expr_data, fetal_har_genes, control_genes):
    """Analyze fetal periods only"""
    
    # Define fetal periods (prenatal: 2-7)
    fetal_periods = [2, 3, 4, 5, 6, 7]
    
    print(f"\nAnalyzing fetal periods: {fetal_periods}")
    
    # Filter for common samples
    common_samples = set(metadata['Sample']) & set(expr_data.columns)
    metadata_filtered = metadata[metadata['Sample'].isin(common_samples)]
    expr_filtered = expr_data[list(common_samples)]
    
    # Check which fetal periods have both species
    analyzable_periods = []
    for period in fetal_periods:
        period_data = metadata_filtered[metadata_filtered['predicted_period'] == period]
        species_counts = period_data['Species'].value_counts()
        
        human_count = species_counts.get('Human', 0)
        macaque_count = species_counts.get('Macaque', 0)
        
        print(f"Period {period}: Human={human_count}, Macaque={macaque_count}")
        
        if human_count > 0 and macaque_count > 0:
            analyzable_periods.append(period)
            print(f"  ✓ Can analyze")
        else:
            print(f"  ✗ Cannot analyze (need both species)")
    
    if not analyzable_periods:
        print("No analyzable fetal periods found!")
        return {}
    
    print(f"\nAnalyzable fetal periods: {analyzable_periods}")
    
    # Average expression across brain regions for each individual
    print(f"\nAveraging expression across brain regions...")
    individual_expr = {}
    
    for (individual, species, period), group in metadata_filtered.groupby(['individual', 'Species', 'predicted_period']):
        if period not in analyzable_periods:
            continue
            
        samples = group['Sample'].tolist()
        if len(samples) > 0:
            avg_expr = expr_filtered[samples].mean(axis=1)
            individual_expr[(individual, species, period)] = avg_expr
    
    print(f"Created {len(individual_expr)} individual profiles")
    
    # Group by species and period
    species_period_data = {'Human': {}, 'Macaque': {}}
    
    for (individual, species, period), expr in individual_expr.items():
        if period not in species_period_data[species]:
            species_period_data[species][period] = []
        species_period_data[species][period].append(expr)
    
    # Calculate mean expression for each period
    species_means = {'Human': {}, 'Macaque': {}}
    
    for species in ['Human', 'Macaque']:
        for period, expr_list in species_period_data[species].items():
            if len(expr_list) > 0:
                mean_expr = pd.DataFrame(expr_list).mean(axis=0)
                species_means[species][period] = mean_expr
    
    print(f"Human periods with data: {sorted(species_means['Human'].keys())}")
    print(f"Macaque periods with data: {sorted(species_means['Macaque'].keys())}")
    
    # Calculate developmental Z-scores
    print(f"\nCalculating developmental Z-scores...")
    species_zscores = {}
    
    for species in ['Human', 'Macaque']:
        if len(species_means[species]) >= 2:
            period_matrix = pd.DataFrame(species_means[species])
            # Z-score across periods for each gene
            zscores = period_matrix.apply(lambda x: (x - x.mean()) / (x.std() + 1e-8), axis=1)
            species_zscores[species] = zscores
    
    # Calculate delta Z-scores
    print(f"Calculating delta Z-scores...")
    
    if 'Human' not in species_zscores or 'Macaque' not in species_zscores:
        print("Missing species Z-score data!")
        return {}
    
    human_z = species_zscores['Human']
    macaque_z = species_zscores['Macaque']
    
    common_periods = set(human_z.columns) & set(macaque_z.columns)
    final_periods = [p for p in analyzable_periods if p in common_periods]
    
    print(f"Final analyzable periods: {sorted(final_periods)}")
    
    delta_zscores = {}
    for period in final_periods:
        delta_z = human_z[period] - macaque_z[period]
        delta_z = delta_z.dropna()
        delta_zscores[period] = delta_z
        print(f"Period {period}: {len(delta_z)} genes, range {delta_z.min():.3f} to {delta_z.max():.3f}")
    
    # Analyze HAR enrichment
    print(f"\nAnalyzing HAR enrichment...")
    results = {}
    
    for period in sorted(delta_zscores.keys()):
        delta_z = delta_zscores[period]
        available_genes = set(delta_z.index)
        
        har_available = fetal_har_genes & available_genes
        control_available = control_genes & available_genes
        
        if len(har_available) < 10 or len(control_available) < 10:
            print(f"Period {period}: Insufficient genes (HAR={len(har_available)}, Control={len(control_available)})")
            continue
        
        har_delta = delta_z[list(har_available)]
        control_delta = delta_z[list(control_available)]
        
        # T-test
        t_stat, p_val = ttest_ind(har_delta, control_delta, equal_var=False)
        
        results[period] = {
            'har_count': len(har_available),
            'control_count': len(control_available),
            'har_mean': har_delta.mean(),
            'har_std': har_delta.std(),
            'control_mean': control_delta.mean(),
            'control_std': control_delta.std(),
            'difference': har_delta.mean() - control_delta.mean(),
            't_stat': t_stat,
            'p_value': p_val,
            'har_values': har_delta,
            'control_values': control_delta
        }
        
        print(f"\nPeriod {period}:")
        print(f"  HAR genes: {len(har_available)}, mean Δ Z = {har_delta.mean():.4f} ± {har_delta.std():.4f}")
        print(f"  Control: {len(control_available)}, mean Δ Z = {control_delta.mean():.4f} ± {control_delta.std():.4f}")
        print(f"  Difference: {results[period]['difference']:.4f}")
        print(f"  t = {t_stat:.3f}, p = {p_val:.2e}")
    
    # Multiple testing correction
    if len(results) > 1:
        periods_list = list(results.keys())
        p_values = [results[p]['p_value'] for p in periods_list]
        
        rejected_fdr, p_adj_fdr, _, _ = multipletests(p_values, method='fdr_bh')
        
        print(f"\nMultiple testing correction (FDR):")
        for i, period in enumerate(periods_list):
            results[period]['p_fdr'] = p_adj_fdr[i]
            results[period]['significant_fdr'] = rejected_fdr[i]
            print(f"Period {period}: p_raw = {p_values[i]:.2e}, p_fdr = {p_adj_fdr[i]:.2e}, significant = {rejected_fdr[i]}")
    
    return results

def create_distribution_curves(results, period_age_map=None):
    """Create beautiful distribution curves for each period"""
    if not results:
        print("No results to visualize")
        return
    
    print("Creating distribution curves...")
    
    # Set font properties globally
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Trebuchet MS', 'Arial', 'DejaVu Sans']
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
    
    # --- NEW, HARMONIOUS COLOR PALETTE (same as last, which blends well) ---
    colors = {
        'har': '#0077b6',       # A clean, rich blue
        'control': '#ee6c4d',   # A clean, vibrant red-orange
        'background': '#F8F9FA'
    }
    
    periods = sorted(results.keys())
    n_periods = len(periods)
    
    # Calculate grid dimensions
    cols = min(3, n_periods)
    rows = (n_periods + cols - 1) // cols
    
    # Create figure
    fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 5*rows))
    fig.patch.set_facecolor('white')
    
    # Handle single subplot case
    if n_periods == 1:
        axes = [axes]
    elif rows == 1:
        axes = [axes] if n_periods == 1 else axes
    else:
        axes = axes.flatten()
    
    for i, period in enumerate(periods):
        ax = axes[i]
        result = results[period]
        
        har_values = result['har_values']
        control_values = result['control_values']
        
        # Calculate range for smooth curves
        x_min = min(har_values.min(), control_values.min()) - 0.3
        x_max = max(har_values.max(), control_values.max()) + 0.3
        x_range = np.linspace(x_min, x_max, 500)
        
        # --- REVISED PLOTTING LOGIC FOR FAINT OVERLAP WITH SUBTLE OUTLINES ---
        # Both fills are semi-transparent.
        # Outlines are added back as thin lines of the same color, with some transparency.

        # Plot Control Curve (Red-Orange) - BASE LAYER
        if len(control_values) > 1:
            try:
                kde_control = gaussian_kde(control_values)
                control_density = kde_control(x_range)
                # Semi-transparent fill
                ax.fill_between(x_range, control_density, color=colors['control'], alpha=0.4, zorder=2)
                # Subtle line outline
                ax.plot(x_range, control_density, color=colors['control'], linewidth=1.0, alpha=0.4, zorder=4)
            except Exception:
                ax.hist(control_values, bins=30, alpha=0.5, color=colors['control'], 
                       density=True, histtype='stepfilled', edgecolor=colors['control'], linewidth=1.0)
        
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
        
        # Add period title in top-left corner
        if period_age_map and period in period_age_map:
            title_text = period_age_map[period]
        else:
            title_text = f'Period {period}'
        
        ax.text(0.02, 0.98, title_text, transform=ax.transAxes,
                fontsize=16, fontweight='bold', fontfamily='Trebuchet MS',
                ha='left', va='top')
        
        # Styling
        ax.set_xlabel('Δ Z-score (Human - Macaque)', fontweight='bold', fontsize=16)
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
        
        # Add adjusted p-value annotation
        ylim = ax.get_ylim()
        y_pos = ylim[1] * 0.02
        
        # Get p-value (use p_fdr if available, otherwise p_value)
        p_val = result.get('p_fdr', result.get('p_value'))
        
        if p_val is not None:
            # Format p-value in scientific notation
            exponent = int(np.floor(np.log10(abs(p_val))))
            mantissa = p_val / (10 ** exponent)
            p_text = f"Adj P-value =\n{mantissa:.2f}×10$^{{{exponent}}}$"
            
            ax.text(0, y_pos, p_text, fontsize=16, fontweight='bold', 
                    fontfamily='Trebuchet MS', ha='center', va='bottom')
    
    # Hide extra subplots
    for i in range(n_periods, len(axes)):
        axes[i].set_visible(False)
    
    # Main title
    fig.suptitle('Differential gene expression associated with positively selected shared (fetal) b-CREs across human fetal brain development.', 
                fontsize=22, fontweight='bold', y=1.0)
    
    # Create proxy artists for a single, clean figure legend (with subtle outlines)
    har_patch = Rectangle((0, 0), 1, 1, fc=colors['har'], alpha=0.5, ec=colors['har'], lw=1.0)
    control_patch = Rectangle((0, 0), 1, 1, fc=colors['control'], alpha=0.5, ec=colors['control'], lw=1.0)
    
    fig.legend(handles=[har_patch, control_patch], 
               labels=['Genes regulated by HS b-CREs.', 'Genes regulated by Non-HS b-CREs.'],
               loc='upper right', 
               bbox_to_anchor=(0.95, 0.98),
               fontsize=14,
               frameon=False)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save with high quality
    plt.savefig('fetal_har_distribution_curves.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.savefig('fetal_har_distribution_curves.pdf', bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.show()
    
    print("Distribution curves saved as:")
    print("- fetal_har_distribution_curves.png")
    print("- fetal_har_distribution_curves.pdf")

def main():
    """Main analysis function"""
    print("fetal HAR DISTRIBUTION CURVES ANALYSIS")
    print("=" * 60)
    
    # Load data
    metadata, expr_data, fetal_har_genes, control_genes, period_age_map = load_data()
    
    if metadata is None:
        print("\nScript stopped because data could not be loaded.")
        return
    
    # Run analysis
    results = analyze_fetal_periods(metadata, expr_data, fetal_har_genes, control_genes)
    
    # Create distribution curves only
    if results:
        create_distribution_curves(results, period_age_map)
    
    print("\nDistribution curves analysis complete!")

if __name__ == "__main__":
    main()