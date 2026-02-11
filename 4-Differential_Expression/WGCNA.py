import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Install sklearn if not available
try:
    from sklearn.linear_model import LogisticRegression
except ImportError:
    import subprocess
    import sys
    print("Installing scikit-learn...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-learn"])
    from sklearn.linear_model import LogisticRegression

def benjamini_hochberg_fdr(p_values, alpha=0.05):
    """Apply Benjamini-Hochberg FDR correction"""
    p_values = np.array(p_values)
    n = len(p_values)
    
    # Sort p-values and keep track of original indices
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    
    # Calculate adjusted p-values
    adjusted_p = np.zeros(n)
    for i in range(n-1, -1, -1):
        if i == n-1:
            adjusted_p[sorted_indices[i]] = sorted_p[i]
        else:
            adjusted_p[sorted_indices[i]] = min(
                sorted_p[i] * n / (i + 1),
                adjusted_p[sorted_indices[i+1]]
            )
    
    return np.clip(adjusted_p, 0, 1)

def run_wgcna_enrichment(df, gene_list, gene_type_name):
    """Run WGCNA enrichment analysis for a specific gene list"""
    
    print(f"\n{'='*70}")
    print(f"ANALYZING {gene_type_name.upper()} GENES")
    print(f"{'='*70}")
    
    # Create HAR status for all genes in WGCNA background
    df['HAR_status'] = df['ID'].isin(gene_list).astype(int)
    genes_in_wgcna = df['HAR_status'].sum()
    print(f"Total {gene_type_name} genes in list: {len(gene_list)}")
    print(f"{gene_type_name} genes found in WGCNA background: {genes_in_wgcna}")
    
    if genes_in_wgcna == 0:
        print(f"WARNING: No {gene_type_name} genes found in WGCNA background!")
        return None
    
    # Test enrichment for each module using logistic linear regression
    print(f"\nTesting enrichment using logistic linear regression across {df['Modules'].nunique()} modules...")
    
    modules = sorted(df['Modules'].unique())
    results = []
    
    for module in modules:
        # Create binary predictor: gene is in this module (1) or not (0)
        in_module = (df['Modules'] == module).astype(int)
        
        # Binary outcome: gene is HAR-associated (1) or not (0)
        har_status = df['HAR_status']
        
        # Skip modules with no genes or if no HAR genes exist
        if in_module.sum() == 0 or har_status.sum() == 0:
            continue
        
        # Logistic linear regression as described in the paper
        X = in_module.values.reshape(-1, 1)
        y = har_status.values
        
        try:
            # Fit logistic regression
            lr = LogisticRegression(fit_intercept=True, max_iter=1000)
            lr.fit(X, y)
            
            # Get coefficient (log odds ratio)
            coef = lr.coef_[0][0]
            
            # Calculate p-value using likelihood ratio test
            lr_null = LogisticRegression(fit_intercept=True, max_iter=1000)
            lr_null.fit(np.zeros((len(y), 1)), y)
            
            # Calculate log-likelihoods
            y_pred_full = lr.predict_proba(X)[:, 1]
            y_pred_null = lr_null.predict_proba(np.zeros((len(y), 1)))[:, 1]
            
            # Avoid numerical issues
            epsilon = 1e-15
            y_pred_full = np.clip(y_pred_full, epsilon, 1-epsilon)
            y_pred_null = np.clip(y_pred_null, epsilon, 1-epsilon)
            
            ll_full = np.sum(y * np.log(y_pred_full) + (1-y) * np.log(1-y_pred_full))
            ll_null = np.sum(y * np.log(y_pred_null) + (1-y) * np.log(1-y_pred_null))
            
            # Likelihood ratio test statistic
            lr_stat = 2 * (ll_full - ll_null)
            p_value = 1 - stats.chi2.cdf(lr_stat, df=1)
            
        except Exception as e:
            print(f"Logistic regression failed for module {module}: {e}")
            continue
        
        # Calculate descriptive statistics
        har_in_module = ((har_status == 1) & (in_module == 1)).sum()
        total_in_module = in_module.sum()
        total_har = har_status.sum()
        
        # Expected number of HAR genes in module under null hypothesis
        expected_har_in_module = (total_har * total_in_module) / len(df)
        fold_enrichment = har_in_module / max(expected_har_in_module, 1)
        
        # Odds ratio
        odds_ratio = np.exp(coef)
        
        results.append({
            'Module': module,
            'Total_genes_in_module': total_in_module,
            'HAR_genes_in_module': har_in_module,
            'Total_HAR_genes': total_har,
            'Expected_HAR_in_module': expected_har_in_module,
            'Fold_enrichment': fold_enrichment,
            'Log_odds_ratio': coef,
            'Odds_ratio': odds_ratio,
            'P_value': p_value
        })
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction
    results_df['P_value_bonferroni'] = results_df['P_value'] * len(results_df)
    results_df['P_value_bonferroni'] = results_df['P_value_bonferroni'].clip(upper=1.0)
    
    # FDR correction
    results_df['P_value_fdr'] = benjamini_hochberg_fdr(results_df['P_value'], alpha=0.05)
    
    # Sort by p-value
    results_df = results_df.sort_values('P_value')
    
    # Display summary
    print(f"\nLogistic linear regression analysis complete!")
    print(f"Modules tested: {len(results_df)}")
    print(f"Significant modules (p < 0.05): {(results_df['P_value'] < 0.05).sum()}")
    print(f"Significant modules (Bonferroni p < 0.05): {(results_df['P_value_bonferroni'] < 0.05).sum()}")
    print(f"Significant modules (FDR p < 0.05): {(results_df['P_value_fdr'] < 0.05).sum()}")
    
    # Show top enriched modules
    print(f"\nTop 10 most significantly enriched modules:")
    display_cols = ['Module', 'HAR_genes_in_module', 'Total_genes_in_module', 
                    'Fold_enrichment', 'Odds_ratio', 'P_value', 'P_value_bonferroni', 'P_value_fdr']
    print(results_df[display_cols].head(10))
    
    return results_df


# Read the Excel file with WGCNA modules
print("Reading WGCNA modules data...")
df = pd.read_excel('Gene Modules.xlsx', header=2)

# Clean the data
df = df.dropna(subset=['ID'])
df['ID'] = df['ID'].astype(str).str.strip()

# Handle merged cells - forward fill module numbers
df['Modules'] = df['Modules'].fillna(method='ffill')
df['Modules'] = df['Modules'].astype(int)

print(f"Total genes in WGCNA: {len(df)}")
print(f"Total modules: {df['Modules'].nunique()}")

# Define gene list files and their names
gene_lists = {
    'shared': 'shared_HPS-bCREs.txt',
    'fetal': 'fetal_HPS-bCREs.txt',
    'adult': 'adult_all_HPS-bCREs.txt'
}

# Process each gene list
for gene_type, filename in gene_lists.items():
    try:
        # Read gene list
        print(f"\n\nReading {gene_type} genes list from {filename}...")
        with open(filename, 'r') as f:
            gene_list = [line.strip() for line in f if line.strip()]
        
        # Run enrichment analysis
        results_df = run_wgcna_enrichment(df.copy(), gene_list, gene_type)
        
        if results_df is not None:
            # Save results
            output_filename = f'HAR_WGCNA_enrichment_results_{gene_type}.csv'
            results_df.to_csv(output_filename, index=False)
            print(f"\nResults saved to '{output_filename}'")
        
    except FileNotFoundError:
        print(f"\nWARNING: Could not find file '{filename}' - skipping {gene_type} analysis")
    except Exception as e:
        print(f"\nERROR processing {gene_type} genes: {e}")

print(f"\n\n{'='*70}")
print("ALL ANALYSES COMPLETE")
print(f"{'='*70}")
print("\nMethod used (as per paper):")
print("- Logistic linear regression")
print(f"- Background: {len(df)} genes included in WGCNA analysis")
print("- No regression of exome or gene length (promoter-based interactions)")
print("- Multiple testing: Bonferroni and FDR corrections applied")