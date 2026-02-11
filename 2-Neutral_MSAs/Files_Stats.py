import os
import pandas as pd
from Bio import SeqIO
import glob
from multiprocessing import Pool, cpu_count

# Path to your folder containing fasta files
folder_path = "Neutral"  # Change if needed

# Function to process a single fasta file
def process_file(fasta_file):
    # Get file name
    file_name = os.path.basename(fasta_file)
    
    # Extract coordinates
    parts = file_name.split('_')
    if len(parts) >= 3:
        chrom = parts[0]
        start = parts[1]
        end = parts[2].split('.')[0]  # Remove .fa extension
        coordinate = f"{chrom}:{start}-{end}"
    else:
        coordinate = "Unknown"
    
    # Count sequences and get length of first sequence
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    num_seq = len(sequences)
    
    if sequences:
        # Get length of first sequence
        length = len(sequences[0].seq)
    else:
        length = 0
        
    return {
        "File": file_name,
        "Coordinates": coordinate,
        "Length": length,
        "Number_of_sequences": num_seq
    }

def main():
    # Get all fasta files in the folder
    fasta_files = glob.glob(os.path.join(folder_path, "*.fa"))
    
    # Use multiprocessing Pool
    num_processes = cpu_count()  # Use all available CPU cores
    print(f"Using {num_processes} processes")
    
    # Process files in parallel
    with Pool(processes=num_processes) as pool:
        results = pool.map(process_file, fasta_files)
    
    # Convert results to DataFrame
    df = pd.DataFrame(results)
    
    # Save to Excel
    output_file = "Neutral_Stats.xlsx"
    df.to_excel(output_file, index=False)
    
    print(f"Analysis complete. Results saved to {output_file}")

if __name__ == "__main__":
    main()
