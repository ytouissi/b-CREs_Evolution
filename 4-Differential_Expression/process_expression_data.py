import pandas as pd
import os
from multiprocessing import Pool, cpu_count

INPUT_DIR = os.getcwd()
OUTPUT_DIR = os.getcwd()

files_config = [
    {
        'input': 'nhp_development_RPKM_rmTechRep.txt',
        'output': 'dev_data_filtered.csv',
        'sep': '\t',
        'index_col': 0
    },
    {
        'input': 'simFiltered_rpkm_combat_nonParametric.txt',
        'output': 'adult_data_filtered.csv',
        'sep': '\t',
        'index_col': 0
    }
]

def convert_file(config):
    """Convert TXT to CSV with CRLF line endings to match original format"""
    input_file = os.path.join(INPUT_DIR, config['input'])
    output_file = os.path.join(OUTPUT_DIR, config['output'])
    
    print(f"Converting: {config['input']}...")
    
    try:
        df = pd.read_csv(input_file, sep=config['sep'], index_col=config['index_col'])
        
        # Extract clean ENSG IDs from index (remove version and gene name)
        clean_index = []
        for idx in df.index:
            # Format: ENSG00000269933.1|RP3-333A15.2 -> ENSG00000269933
            ensg_id = str(idx).split('|')[0].split('.')[0]
            clean_index.append(ensg_id)
        
        df.index = clean_index
        
        # Remove ProbeID column if it exists
        if 'ProbeID' in df.columns:
            df = df.drop(columns=['ProbeID'])
        
        # Round all numeric columns to 9 decimal places to match original
        numeric_cols = df.select_dtypes(include=['float64']).columns
        for col in numeric_cols:
            df[col] = df[col].round(9)
        
        # Save with CRLF line endings to match original format
        df.to_csv(output_file, lineterminator='\r\n')
        
        file_size = os.path.getsize(output_file) / (1024*1024)
        print(f"✓ {config['output']} - {df.shape[0]} rows × {df.shape[1]} columns - {file_size:.2f} MB")
        
        return {'status': 'success', 'file': config['output'], 'rows': df.shape[0], 'cols': df.shape[1]}
    
    except Exception as e:
        print(f"✗ Error: {e}")
        return {'status': 'error', 'file': config['input'], 'error': str(e)}

def main():
    print("="*70)
    print("Converting TXT to CSV (with CRLF line endings)")
    print("="*70)
    
    with Pool(processes=min(2, cpu_count())) as pool:
        results = pool.map(convert_file, files_config)
    
    print("\n" + "="*70)
    for result in results:
        if result['status'] == 'success':
            print(f"✓ {result['file']}: {result['rows']:,} rows × {result['cols']} cols")
    print("="*70)

if __name__ == '__main__':
    main()