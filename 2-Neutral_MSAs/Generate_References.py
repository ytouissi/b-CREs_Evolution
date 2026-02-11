import os
import traceback
from Bio import AlignIO
import pandas as pd
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
import time
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("reference_creation.log"),
        logging.StreamHandler()
    ]
)

# Directory paths
query_dir = "../1-b-CREs_MSAs/Adaptify_MSAs"  # Directory containing all query .fa files
reference_dir = "Neutral"  # Directory containing all potential reference MSA files
output_dir = "References"  # Directory for output files
excel_filter_file = "Neutral_Stats.xlsx"

# Create directories if they don't exist
for directory in [query_dir, reference_dir, output_dir]:
    if not os.path.exists(directory):
        os.makedirs(directory)

# Output Excel file
output_excel = os.path.join(output_dir, "reference_summary.xlsx")

# Target alignment length
TARGET_LENGTH = 3000

# Cache for alignment information to avoid repeated file reads
alignment_cache = {}

# Function to get species IDs from an MSA file with caching
def get_alignment_info(msa_file):
    """Get species IDs and length from an MSA file with caching"""
    if msa_file in alignment_cache:
        return alignment_cache[msa_file]
    
    try:
        alignment = AlignIO.read(msa_file, 'fasta')
        species_ids = set(record.id for record in alignment)
        alignment_length = alignment.get_alignment_length()
        
        # Verify all sequences have the same length (properly aligned)
        seq_lengths = [len(record.seq) for record in alignment]
        is_aligned = all(length == alignment_length for length in seq_lengths)
        
        if not is_aligned:
            logging.warning(f"MSA file {os.path.basename(msa_file)} has sequences of different lengths!")
            return None
            
        result = (species_ids, alignment_length)
        alignment_cache[msa_file] = result
        return result
        
    except Exception as e:
        logging.error(f"Error reading {msa_file}: {e}")
        return None

# Find suitable reference files for a query
def find_suitable_references(query_species, reference_files, species_df=None):
    """Find reference files that contain ALL of the query species"""
    suitable_refs = []
    
    # First check species counts from Excel if available
    if species_df is not None and 'File' in species_df.columns and 'Number_of_sequences' in species_df.columns:
        file_info = dict(zip(species_df['File'], species_df['Number_of_sequences']))
    else:
        file_info = {}
    
    # Check each reference file
    for ref_file in reference_files:
        ref_basename = os.path.basename(ref_file)
        
        # Quick filter using Excel data if available
        if ref_basename in file_info and file_info[ref_basename] < len(query_species):
            continue
            
        # Get species in this reference
        alignment_info = get_alignment_info(ref_file)
        if alignment_info is None:
            continue
            
        ref_species, ref_length = alignment_info
        
        # Check if this reference contains ALL query species
        if query_species.issubset(ref_species):
            suitable_refs.append((ref_file, ref_length))
            logging.debug(f"Reference {ref_basename} contains all {len(query_species)} query species")
        else:
            missing = len(query_species - ref_species)
            logging.debug(f"Reference {ref_basename} missing {missing} query species, skipping")
    
    return suitable_refs

# Function to concatenate MSA files horizontally
def concatenate_msas(query_species, reference_files, query_records, target_length=TARGET_LENGTH):
    """
    Concatenate multiple MSA files horizontally for the species in query,
    maintaining original sequence order
    """
    # Initialize dictionary to store concatenated sequences by species ID
    concatenated_seqs = {record.id: "" for record in query_records}
    
    # Track reference files used and total length
    used_ref_files = []
    total_length = 0
    
    # Process reference files until we reach target length
    for ref_file, ref_length in reference_files:
        # Check if we've reached target length
        if total_length >= target_length:
            break
            
        # Calculate how much to use from this alignment
        if total_length + ref_length > target_length:
            usable_length = target_length - total_length
        else:
            usable_length = ref_length
        
        try:
            # Read the reference alignment
            alignment = AlignIO.read(ref_file, 'fasta')
            
            # Map species IDs to records for faster lookup
            species_records = {record.id: record for record in alignment}
            
            # Add sequences from this reference for ALL query species
            all_found = True
            for record in query_records:
                species_id = record.id
                if species_id in species_records:
                    # Add sequence (trimmed if needed)
                    sequence = str(species_records[species_id].seq[:usable_length])
                    concatenated_seqs[species_id] += sequence
                else:
                    # This should never happen given our filtering above
                    logging.error(f"Species {species_id} missing from reference {ref_file}!")
                    all_found = False
            
            # Only use this reference if ALL query species were found
            if all_found:
                used_ref_files.append((ref_file, usable_length))
                total_length += usable_length
                logging.debug(f"Added {usable_length} bp from {os.path.basename(ref_file)}, total: {total_length}")
            
        except Exception as e:
            logging.error(f"Error processing reference {ref_file}: {e}")
    
    # Create sequence records for output, preserving original query order
    records = []
    for record in query_records:
        species_id = record.id
        seq = concatenated_seqs[species_id]
        # Verify we have sequence data for this species
        if seq:
            records.append(SeqRecord(Seq(seq), id=species_id, description=""))
        else:
            logging.warning(f"No sequence data for species {species_id}!")
    
    return records, used_ref_files, total_length

# Process a single query
def process_query(query_file, reference_files, species_df=None):
    query_name = os.path.basename(query_file)
    logging.info(f"Processing query: {query_name}")
    
    results = {
        'excel_data': None,
        'success': False,
        'species_count': 0,
        'total_length': 0
    }
    
    try:
        # Get species from query
        try:
            query_alignment = AlignIO.read(query_file, 'fasta')
            query_records = list(query_alignment)  # Keep original records to preserve order
            query_species = set(record.id for record in query_records)
        except Exception as e:
            logging.error(f"Error reading query {query_name}: {e}")
            return results
        
        if not query_species:
            logging.error(f"No species found in query {query_name}")
            return results
        
        logging.info(f"Query {query_name} has {len(query_species)} species")
        
        # Find suitable reference files (that contain ALL query species)
        suitable_refs = find_suitable_references(query_species, reference_files, species_df)
        
        if not suitable_refs:
            logging.error(f"No suitable reference files found for {query_name} (need files with all {len(query_species)} query species)")
            return results
        
        logging.info(f"Found {len(suitable_refs)} suitable reference files for {query_name}")
        
        # Select up to 11 reference files (randomly)
        if len(suitable_refs) > 11:
            selected_refs = random.sample(suitable_refs, 11)
            # Add remaining in case we need them to reach target length
            remaining = [r for r in suitable_refs if r not in selected_refs]
            all_refs = selected_refs + remaining
        else:
            all_refs = suitable_refs
        
        # Concatenate MSAs
        records, used_files, total_length = concatenate_msas(query_species, all_refs, query_records, TARGET_LENGTH)
        
        if not records or len(records) != len(query_species):
            logging.error(f"Failed to get sequences for all query species in {query_name}")
            return results
        
        # Write output file
        output_file = os.path.join(output_dir, query_name)
        with open(output_file, 'w') as out_handle:
            for record in records:
                out_handle.write(f">{record.id}\n")
                out_handle.write(f"{str(record.seq)}\n")
        
        # Prepare data for Excel summary - ONE ROW PER QUERY with comma-separated references
        ref_names = [os.path.basename(ref[0]) for ref in used_files]
        ref_string = ", ".join(ref_names)
        
        results['excel_data'] = {
            'Query': query_name,
            'References_Used': ref_string,
            'Number_Of_References': len(used_files),
            'Total_Length': total_length,
            'Species_Count': len(records)
        }
        
        results['success'] = True
        results['species_count'] = len(records)
        results['total_length'] = total_length
        
        logging.info(f"Created concatenated alignment for {query_name} with {len(records)} species and {total_length} positions")
        
    except Exception as e:
        logging.error(f"Error processing {query_name}: {e}")
        logging.error(traceback.format_exc())
    
    return results

def main():
    start_time = time.time()
    logging.info("Starting horizontal concatenation process")
    
    # Read Excel species count data if available
    species_df = None
    if os.path.exists(excel_filter_file):
        try:
            species_df = pd.read_excel(excel_filter_file)
            logging.info(f"Read Excel file with {len(species_df)} entries")
        except Exception as e:
            logging.error(f"Error reading Excel file: {e}")
    else:
        logging.warning(f"Excel file {excel_filter_file} not found")
    
    # Find all query files
    query_files = []
    if os.path.exists(query_dir):
        for root, _, files in os.walk(query_dir):
            for file in files:
                if file.endswith('.fa'):
                    query_files.append(os.path.join(root, file))
    else:
        logging.error(f"Query directory {query_dir} does not exist!")
        return
    
    if not query_files:
        logging.error("No query files found!")
        return
    
    logging.info(f"Found {len(query_files)} query files")
    
    # Find all reference MSA files
    reference_files = []
    if os.path.exists(reference_dir):
        for root, _, files in os.walk(reference_dir):
            for file in files:
                if file.endswith('.fa'):
                    reference_files.append(os.path.join(root, file))
    else:
        logging.error(f"Reference directory {reference_dir} does not exist!")
        return
    
    if not reference_files:
        logging.error("No reference MSA files found!")
        return
    
    logging.info(f"Found {len(reference_files)} reference MSA files")
    
    # Process queries in parallel
    num_cores = max(1, mp.cpu_count())
    logging.info(f"Using {num_cores} CPU cores for parallel processing")
    
    all_results = []
    with mp.Pool(processes=num_cores) as pool:
        all_results = pool.starmap(
            process_query,
            [(q, reference_files, species_df) for q in query_files]
        )
    
    # Collect Excel data - one row per query
    all_excel_data = []
    successful_queries = 0
    
    for result in all_results:
        if result['success'] and result['excel_data']:
            successful_queries += 1
            all_excel_data.append(result['excel_data'])
    
    # Write Excel summary
    if all_excel_data:
        df = pd.DataFrame(all_excel_data)
        df.to_excel(output_excel, index=False)
        logging.info(f"Excel summary written to {output_excel}")
    
    end_time = time.time()
    logging.info(f"Successfully processed {successful_queries}/{len(query_files)} queries")
    logging.info(f"Total processing time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    # Ensure proper multiprocessing behavior on Windows
    mp.freeze_support()
    try:
        main()
    except Exception as e:
        logging.critical(f"Fatal error: {e}")
        logging.critical(traceback.format_exc())
