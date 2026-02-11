#!/usr/bin/env python3
import os
import multiprocessing as mp
from pathlib import Path

def count_sequences(file_path):
    with open(file_path, 'r') as f:
        count = sum(1 for line in f if line.startswith('>'))
    
    filename = Path(file_path).stem
    parts = filename.split('_')
    
    if len(parts) >= 3:
        return {
            'chrom': parts[0],
            'start': parts[1],
            'end': parts[2],
            'count': count
        }
    return None

def main():
    # Find all folders ending with _MSAs_Ready
    intermediate_dir = os.path.join(os.getcwd(), "Intermediate_Files")
    
    if not os.path.exists(intermediate_dir):
        print(f"Intermediate_Files folder not found at {intermediate_dir}")
        return
    
    ready_folders = [f for f in os.listdir(intermediate_dir) 
                     if os.path.isdir(os.path.join(intermediate_dir, f)) and f.endswith("_MSAs_Ready")]
    
    if not ready_folders:
        print("No folders ending with '_MSAs_Ready' found in Intermediate_Files directory")
        return
    
    # Process each ready folder
    for folder_name in ready_folders:
        input_folder = os.path.join(intermediate_dir, folder_name)
        
        # Determine threshold based on folder name
        if "Cactus" in folder_name or "cactus" in folder_name:
            THRESHOLD = 60
        elif "Multiz" in folder_name or "multiz" in folder_name:
            THRESHOLD = 27
        else:
            # Skip UCSC30 and any other folders
            print(f"Skipping {folder_name} (not Cactus or Multiz)")
            continue
        
        # Generate output bed file name based on input folder
        folder_prefix = folder_name.replace("_MSAs_Ready", "")
        bed_file = os.path.join(intermediate_dir, f'{folder_prefix}_regions_below_{THRESHOLD}.bed')
        
        # Check if bed file already exists
        if os.path.exists(bed_file):
            print(f"File '{bed_file}' already exists - skipping {folder_name}")
            continue
        
        print(f"\nProcessing {folder_name} with threshold {THRESHOLD}...")
        
        # Get all .fa files
        fasta_files = [os.path.join(input_folder, f) 
                       for f in os.listdir(input_folder) 
                       if f.endswith('.fa')]
        
        if not fasta_files:
            print(f"No .fa files found in {input_folder}")
            continue
        
        # Count sequences in parallel
        num_cores = mp.cpu_count()
        
        with mp.Pool(processes=num_cores) as pool:
            results = pool.map(count_sequences, fasta_files)
        
        results = [r for r in results if r is not None]
        
        # Write regions below threshold to bed file
        below_threshold_count = 0
        with open(bed_file, 'w') as f:
            for region in results:
                if region['count'] < THRESHOLD:
                    f.write(f"{region['chrom']}\t{region['start']}\t{region['end']}\t{region['count']}\n")
                    below_threshold_count += 1
        
        print(f"Total regions analyzed: {len(results)}")
        print(f"Regions below threshold ({THRESHOLD}): {below_threshold_count}")
        print(f"Results saved to: {bed_file}")

if __name__ == "__main__":
    main()