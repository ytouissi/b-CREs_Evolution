from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os
import multiprocessing as mp
from functools import partial

def process_alignment_file(file_name, input_folder, output_folder, coverage_threshold=0.7):
    """Process a single alignment file."""
    try:
        file_path = os.path.join(input_folder, file_name)
        
        # Load alignment
        alignment = AlignIO.read(file_path, "fasta")
        
        # Replace * with -
        for record in alignment:
            record.seq = Seq(str(record.seq).replace("*", "-").replace("N", "-"))
        
        # Assume first sequence is the reference (e.g., hg38)
        reference = alignment[0]
        ref_length = len(str(reference.seq).replace("-", ""))
        
        # Filter sequences with coverage >= threshold
        seqs_to_keep = []
        for record in alignment:
            seq_length = len(str(record.seq).replace("-", ""))
            if ref_length > 0:
                coverage = seq_length / ref_length
            else:
                coverage = 0
            if coverage >= coverage_threshold:
                seqs_to_keep.append(record)
        
        # Create new alignment with filtered sequences
        cleaned_alignment = MultipleSeqAlignment(seqs_to_keep)
        
        # Save cleaned alignment
        output_path = os.path.join(output_folder, f"cleaned_{file_name}")
        AlignIO.write(cleaned_alignment, output_path, "fasta")
        
        return f"Processed {file_name} successfully"
    except Exception as e:
        return f"Error processing {file_name}: {str(e)}"

def main():
    # Directory with .fa files
    folder_path = "Intermediate_Files/fasta"  # Replace with your folder path
    output_folder = "Intermediate_Files/cleaned"  # Output folder
    os.makedirs(output_folder, exist_ok=True)
    
    # Get all .fa files
    fa_files = [f for f in os.listdir(folder_path) if f.endswith(".fa")]
    
    # Set up the process pool
    num_cores = mp.cpu_count()
    print(f"Using {num_cores} CPU cores")
    
    # Create a partial function with the fixed arguments
    process_file = partial(
        process_alignment_file, 
        input_folder=folder_path,
        output_folder=output_folder
    )
    
    # Process files in parallel
    with mp.Pool(processes=num_cores) as pool:
        results = pool.map(process_file, fa_files)
    
    # Print results
    for result in results:
        print(result)

if __name__ == "__main__":
    main()
