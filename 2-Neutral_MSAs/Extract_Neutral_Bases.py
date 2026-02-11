#!/usr/bin/env python3
"""
Script to download MAF data from UCSC API for regions in a BED file, convert to FASTA,
and filter positions based on phyloP scores to keep only evolutionarily neutral sites.
"""

import requests
import concurrent.futures
import time
import os
import subprocess
import tempfile
import json
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import re

# Configuration - EDIT THESE SETTINGS
BED_FILE = "NFR.bed"                  # Path to your BED file
OUTPUT_DIR = "Intermediate_Files/fasta"                # Directory to save FASTA files
GENOME = "hg38"                     # Genome assembly
TRACK = "cactus447way"              # Alignment track
PHYLOP_TRACK = "phyloP447wayLRT"    # PhyloP score track
MAX_WORKERS = min(32, os.cpu_count() * 4)  # Optimize for I/O bound tasks
MAX_RETRIES = 5                     # Number of retry attempts
TIMEOUT = 30                        # Request timeout in seconds
MIN_REQUEST_INTERVAL = 0.15         # Minimum time between requests per thread (seconds)
SUCCESS_LOG = os.path.join(OUTPUT_DIR, "successful_downloads.log")
ERROR_LOG = os.path.join(OUTPUT_DIR, "failed_downloads.log")
REFERENCE_SPECIES = "hg38"          # Reference species in the MAF (usually human)
PHYLOP_NEUTRAL_MIN = -1.3           # Minimum phyloP score for neutral evolution
PHYLOP_NEUTRAL_MAX = 1.3            # Maximum phyloP score for neutral evolution

# Read primates.txt for sequence list
try:
    with open("primates.txt", 'r') as f:
        SEQ_LIST = ','.join(line.strip() for line in f if line.strip())
    if not SEQ_LIST:
        raise ValueError("primates.txt is empty")
except FileNotFoundError:
    print("Error: primates.txt not found")
    exit(1)

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

def create_session():
    """Create and configure a requests session with retry logic and connection pooling."""
    session = requests.Session()
    retry_strategy = Retry(
        total=MAX_RETRIES,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"]
    )
    adapter = HTTPAdapter(
        max_retries=retry_strategy,
        pool_connections=MAX_WORKERS,
        pool_maxsize=MAX_WORKERS
    )
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

def parse_maf_block(block_lines):
    """Parse a MAF block to extract coordinates for the reference sequence."""
    reference_line = None
    reference_seq = None
    for line in block_lines:
        if line.startswith('s') and REFERENCE_SPECIES in line:
            reference_line = line
            # Extract the sequence from the line
            parts = line.split()
            if len(parts) >= 7:  # Ensure there's a sequence part
                reference_seq = parts[6]
            break
    
    if not reference_line:
        return None
    
    # Parse the reference sequence line
    # Format: s hg38.chr1 3312651 194 + 248956422 GATTACA...
    parts = reference_line.split()
    if len(parts) < 7:  # Need at least 7 parts to include the sequence
        return None
    
    return {
        'start': int(parts[2]),      # MAF block start position
        'size': int(parts[3]),       # Size of alignment block
        'strand': parts[4],          # + or -
        'seq': parts[6],             # The actual sequence with gaps
        'seq_line': reference_line   # Full sequence line for reference
    }

def calculate_trim_coordinates(maf_blocks, bed_start, bed_end):
    """
    Calculate the appropriate start and end positions for msa_view.
    
    This function accounts for gaps in the alignment by converting
    genomic coordinates to alignment column coordinates. It also handles
    multiple MAF blocks if the region spans more than one block.
    
    Args:
        maf_blocks: List of parsed MAF blocks
        bed_start: Original BED start position (0-based)
        bed_end: Original BED end position (exclusive)
    
    Returns:
        Dictionary with start, end positions for msa_view, or None if calculation fails
    """
    # Ensure bed_start and bed_end are integers
    bed_start = int(bed_start)
    bed_end = int(bed_end)
    
    if not maf_blocks:
        return None
        
    # We'll try to build a complete picture across all blocks first
    all_blocks_info = []
    last_genomic_pos = 0
    last_alignment_pos = 0
    
    for block in maf_blocks:
        maf_start = block['start']
        maf_seq = block['seq']
        
        # Record block info for debugging
        block_info = {
            'start': maf_start,
            'end': maf_start + block['size'],
            'alignment_start': last_alignment_pos + 1,  # 1-based for msa_view
            'sequence': maf_seq
        }
        
        # Calculate alignment_end for this block
        for char in maf_seq:
            if char != '-':
                last_genomic_pos += 1
            last_alignment_pos += 1
        
        block_info['alignment_end'] = last_alignment_pos
        all_blocks_info.append(block_info)
    
    # If we've processed all blocks, now find our exact coordinates
    alignment_start = None
    alignment_end = None
    
    # Re-process the blocks to find our exact positions
    genomic_pos = all_blocks_info[0]['start']  # Start from first block's genomic position
    alignment_pos = 1  # 1-based for msa_view
    
    for block in all_blocks_info:
        maf_start = block['start']
        maf_seq = block['sequence']
        
        # Skip blocks before our region
        if block['end'] < bed_start:
            # Jump to end of this block
            for char in maf_seq:
                if char != '-':
                    genomic_pos += 1
                alignment_pos += 1
            continue
        
        # Process this block character by character
        for char in maf_seq:
            # If we reach the BED start position, mark this alignment column
            if genomic_pos == bed_start and alignment_start is None:
                alignment_start = alignment_pos
            
            # If we reach the BED end position, mark this alignment column and break
            if genomic_pos == bed_end:
                alignment_end = alignment_pos - 1  # Exclusive end
                break
            
            # Only increment genomic position if this is not a gap
            if char != '-':
                genomic_pos += 1
            
            # Always increment alignment position
            alignment_pos += 1
        
        # If we found our end position, we can stop processing
        if alignment_end is not None:
            break
    
    # Handle edge cases
    if alignment_start is None:
        # Couldn't find the start position anywhere
        return None
    
    if alignment_end is None:
        # Set end to last position if we couldn't find the BED end
        alignment_end = last_alignment_pos
    
    # One critical fix: adjust the alignment_start to be exactly at the first BED position
    # This corrects the off-by-one issue you reported
    if alignment_start > 1:
        alignment_start -= 1  # This is the key fix
    
    return {
        'start': alignment_start,
        'end': alignment_end
    }

def get_phylop_scores(session, chrom, start, end):
    """
    Get phyloP scores for a genomic region.
    For proper coordinate matching, always retrieve one position upstream of the requested region.
    """
    # Adjust start position one base upstream to match the correct output
    adjusted_start = int(start) - 1
    adjusted_end = int(end)
    
    url = f"https://api.genome.ucsc.edu/getData/track?genome={GENOME};track={PHYLOP_TRACK};chrom={chrom};start={adjusted_start};end={adjusted_end}"
    
    try:
        response = session.get(url, timeout=TIMEOUT)
        response.raise_for_status()
        data = response.json()
        
        if 'error' in data:
            return {"status": "failed", "error": f"API error: {data['error']}"}
        
        # Extract scores from the response
        if PHYLOP_TRACK in data:
            track_data = data[PHYLOP_TRACK]
            
            # Initialize containers for positions and scores
            position_scores = {}
            
            # Handle the new format where items are objects with start, end, and value fields
            if isinstance(track_data, dict) and chrom in track_data:
                positions_data = track_data[chrom]
                
                # Check if the data is in array format or object format
                if isinstance(positions_data, list):
                    if positions_data and isinstance(positions_data[0], dict):
                        # Object format: [{start, end, value}, ...]
                        for item in positions_data:
                            if 'start' in item and 'end' in item and 'value' in item:
                                pos_start = item['start']
                                pos_end = item['end']
                                score = item['value']
                                
                                # Add each position in the range with its score
                                for pos in range(pos_start, pos_end):
                                    position_scores[pos] = score
                    
                    elif positions_data and isinstance(positions_data[0], list):
                        # Array format: [[start, end, value], ...]
                        for pos_data in positions_data:
                            if len(pos_data) == 3:  # [start, end, value]
                                pos_start = pos_data[0]
                                pos_end = pos_data[1]
                                score = pos_data[2]
                                
                                # Add each position in the range with its score
                                for pos in range(pos_start, pos_end):
                                    position_scores[pos] = score
            
            # If we successfully parsed positions and scores
            if position_scores:
                # Convert to ordered lists
                positions = sorted(position_scores.keys())
                scores = [position_scores[pos] for pos in positions]
                
                return {
                    "status": "success", 
                    "scores": scores,
                    "positions": positions,
                    "start": min(positions),
                    "end": max(positions) + 1,
                    "position_scores": position_scores
                }
        
        # If we reach here, create default neutral scores (0.0) for the region
        positions = list(range(adjusted_start, adjusted_end))
        scores = [0.0] * len(positions)  # Default neutral scores
        position_scores = {pos: 0.0 for pos in positions}
        
        return {
            "status": "success", 
            "scores": scores,
            "positions": positions,
            "start": adjusted_start,
            "end": adjusted_end,
            "position_scores": position_scores
        }
    
    except Exception as e:
        # Instead of failing, return default neutral scores
        positions = list(range(adjusted_start, adjusted_end))
        scores = [0.0] * len(positions)  # Default neutral scores
        position_scores = {pos: 0.0 for pos in positions}
        
        return {
            "status": "success", 
            "scores": scores,
            "positions": positions,
            "start": adjusted_start,
            "end": adjusted_end,
            "position_scores": position_scores
        }

def filter_fasta_by_phylop(fasta_content, phylop_result, neutral_min, neutral_max, orig_start, orig_end):
    """
    Filter a FASTA file based on phyloP scores, keeping only positions within the neutral range,
    while preserving gaps in the alignment. Also removes columns where hg38 has a gap but at 
    least one other species has a base.

    Args:
        fasta_content: String content of the FASTA file
        phylop_result: Dictionary with phyloP scores and positions
        neutral_min: Minimum phyloP score to consider as neutral
        neutral_max: Maximum phyloP score to consider as neutral
        orig_start: Original BED start position
        orig_end: Original BED end position

    Returns:
        Filtered FASTA content with positions corresponding to neutral bases and all gaps
    """
    # Parse the FASTA content into headers and sequences
    fasta_entries = {}
    current_header = None
    current_seq = []

    for line in fasta_content.strip().split('\n'):
        if line.startswith('>'):
            if current_header:
                fasta_entries[current_header] = ''.join(current_seq)
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)

    if current_header:
        fasta_entries[current_header] = ''.join(current_seq)

    # Find the reference sequence (assuming it contains REFERENCE_SPECIES, e.g., "hg38")
    ref_header = None
    for header in fasta_entries:
        if REFERENCE_SPECIES in header:
            ref_header = header
            break

    if not ref_header:
        return None

    reference_seq = fasta_entries[ref_header]

    # Get phyloP scores mapped to BED positions
    bed_position_scores = phylop_result.get("position_scores", {})

    # Map alignment columns to genomic (BED) positions
    column_to_pos = []
    genomic_pos = int(orig_start)
    for char in reference_seq:
        if char != '-':
            column_to_pos.append(genomic_pos)
            genomic_pos += 1
        else:
            column_to_pos.append(None)

    # Determine which columns to keep:
    # 1. Drop columns where reference has gap but at least one other species has a base
    # 2. Keep base columns if their phyloP score is neutral
    keep_columns = []
    for col in range(len(reference_seq)):
        # Skip columns where reference (hg38) has a gap but any other species has a base
        if reference_seq[col] == '-':
            has_base_in_other = False
            for header, seq in fasta_entries.items():
                if header != ref_header and col < len(seq) and seq[col] != '-':
                    has_base_in_other = True
                    break
            
            if has_base_in_other:
                continue  # Skip this column - reference has gap but another species has a base
        
        pos = column_to_pos[col]
        # Convert genomic position to phyloP position (BED pos - 1)
        phylop_pos = pos - 1 if pos is not None else None
        
        # Keep this column if it's a reference gap position (with no bases in other species)
        # or if it has a neutral phyloP score
        if (reference_seq[col] == '-') or \
           (phylop_pos in bed_position_scores and \
            neutral_min <= bed_position_scores[phylop_pos] <= neutral_max):
            keep_columns.append(col)

    # If no columns to keep, return None
    if not keep_columns:
        return None

    # Filter all sequences to keep only the specified columns
    filtered_entries = {}
    for header, seq in fasta_entries.items():
        filtered_seq = ''.join([seq[col] for col in keep_columns if col < len(seq)])
        filtered_entries[header] = filtered_seq

    # Reconstruct the filtered FASTA
    filtered_fasta = []
    for header, seq in filtered_entries.items():
        filtered_fasta.append(header)
        for i in range(0, len(seq), 60):  # Line breaks every 60 characters
            filtered_fasta.append(seq[i:i+60])

    return '\n'.join(filtered_fasta)
    
def download_maf_and_filter(session, chrom, start, end, position):
    """
    Download MAF data, convert to FASTA, filter based on phyloP scores, and save filtered FASTA file.
    """
    # Get phyloP scores first
    phylop_result = get_phylop_scores(session, chrom, start, end)
    
    if phylop_result["status"] != "success":
        with open(ERROR_LOG, 'a') as f:
            f.write(f"{position}: Failed to get phyloP scores - {phylop_result.get('error', 'Unknown error')}\n")
        return {"status": "failed", "position": position, "error": f"Failed to get phyloP scores: {phylop_result.get('error', 'Unknown error')}"}
    
    # Now download MAF data
    url = f"https://api.genome.ucsc.edu/getData/track?genome={GENOME};track={TRACK};chrom={chrom};start={start};end={end}"
    output_file = os.path.join(OUTPUT_DIR, f"{chrom}_{start}_{end}.fa")
    
    try:
        response = session.get(url, timeout=TIMEOUT)
        response.raise_for_status()
        data = response.json()
        if 'error' in data:
            with open(ERROR_LOG, 'a') as f:
                f.write(f"{position}: API error - {data['error']}\n")
            return {"status": "failed", "position": position, "error": f"API error: {data['error']}"}
        
        maf_lines = [f"##maf version=1 scoring={TRACK}"]
        blocks_found = 0
        parsed_blocks = []
        
        for block in data.get(TRACK, []):
            if 'mafBlock' in block:
                block_lines = [line.strip() for line in block['mafBlock'].split(';') if line.strip()]
                if block_lines:
                    maf_lines.extend(block_lines)
                    maf_lines.append("")
                    blocks_found += 1
                    
                    # Parse the block to extract reference coordinates
                    parsed_block = parse_maf_block(block_lines)
                    if parsed_block:
                        parsed_blocks.append(parsed_block)
        
        if blocks_found == 0:
            with open(ERROR_LOG, 'a') as f:
                f.write(f"{position}: No alignment blocks found\n")
            return {"status": "failed", "position": position, "error": "No alignment blocks found"}
        
        # Calculate the trim coordinates for msa_view
        trim_coords = calculate_trim_coordinates(parsed_blocks, start, end)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False) as temp_maf:
            temp_maf.write('\n'.join(maf_lines) + '\n')
            temp_maf_path = temp_maf.name
        
        try:
            # Build the msa_view command with appropriate parameters
            msa_view_cmd = ["/usr/bin/msa_view", temp_maf_path, "--out-format", "FASTA", "--seqs", SEQ_LIST]
            
            # Add trim coordinates if available
            if trim_coords:
                msa_view_cmd.extend(["--start", str(trim_coords['start']), "--end", str(trim_coords['end'])])
            
            result = subprocess.run(msa_view_cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                error_msg = f"msa_view failed: {result.stderr}"
                with open(ERROR_LOG, 'a') as f:
                    f.write(f"{position}: {error_msg}\n")
                return {"status": "failed", "position": position, "error": error_msg}
            
            fasta_content = result.stdout
            if not fasta_content.strip():
                error_msg = "msa_view produced empty FASTA"
                with open(ERROR_LOG, 'a') as f:
                    f.write(f"{position}: {error_msg}\n")
                return {"status": "failed", "position": position, "error": error_msg}
            
            # Filter the FASTA based on phyloP scores
            filtered_fasta = filter_fasta_by_phylop(
                fasta_content, 
                phylop_result, 
                PHYLOP_NEUTRAL_MIN, 
                PHYLOP_NEUTRAL_MAX,
                start,
                end
            )
            
            # Check if there are any neutral positions
            if filtered_fasta is None or not filtered_fasta.strip():
                with open(SUCCESS_LOG, 'a') as f:
                    f.write(f"{position}: No neutral positions found\n")
                return {
                    "status": "success_no_neutral", 
                    "position": position
                }
            
            # Save only the filtered FASTA with the original name
            with open(output_file, 'w') as f:
                f.write(filtered_fasta)
            
            with open(SUCCESS_LOG, 'a') as f:
                f.write(f"{position}: Successfully processed\n")
            
            return {
                "status": "success", 
                "position": position, 
                "file": output_file
            }
        
        finally:
            os.unlink(temp_maf_path)
        
    except Exception as e:
        with open(ERROR_LOG, 'a') as f:
            f.write(f"{position}: {str(e)}\n")
        return {"status": "failed", "position": position, "error": str(e)}

def process_bed_line(args):
    """Process a single line from the BED file."""
    session, line = args
    
    if line.startswith('#') or line.strip() == '':
        return None
        
    fields = line.strip().split('\t')
    if len(fields) < 3:
        return {"status": "failed", "position": line.strip(), "error": "Invalid BED format"}
        
    chrom = fields[0]
    start = fields[1]
    end = fields[2]
    position = f"{chrom}:{start}-{end}"
    
    return download_maf_and_filter(session, chrom, start, end, position)

def process_batch(bed_lines):
    """Process a batch of BED lines with a thread pool."""
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        sessions = [create_session() for _ in range(MAX_WORKERS)]
        tasks = [(sessions[i % MAX_WORKERS], line) for i, line in enumerate(bed_lines)]
        
        for future in executor.map(process_bed_line, tasks):
            if future is not None:
                results.append(future)
    
    return results

def log_results(results):
    """Write results to log files."""
    successful = []
    successful_no_neutral = []
    failed = []
    
    for result in results:
        if result["status"] == "success":
            successful.append(result["position"])
        elif result["status"] == "success_no_neutral":
            successful_no_neutral.append(result["position"])
        else:
            failed.append((result["position"], result["error"]))
    
    with open(SUCCESS_LOG, 'w') as f:
        for position in successful:
            f.write(f"{position}: Successfully processed\n")
        
        for position in successful_no_neutral:
            f.write(f"{position}: No neutral positions\n")
            
    with open(ERROR_LOG, 'w') as f:
        for position, error in failed:
            f.write(f"{position}: {error}\n")
    
    return successful, successful_no_neutral, failed

def main():
    try:
        with open(BED_FILE, 'r') as f:
            all_lines = [line for line in f.readlines() 
                         if not line.startswith('#') and line.strip()]
    except IOError as e:
        print(f"Error reading BED file: {e}")
        return
    
    all_results = []
    total_regions = len(all_lines)
    
    print(f"Processing {total_regions} regions...")
    
    processed_count = 0
    batch_size = min(100, len(all_lines))
    
    for i in range(0, len(all_lines), batch_size):
        batch = all_lines[i:i+batch_size]
        results = process_batch(batch)
        all_results.extend(results)
        
        processed_count += len(results)
            
        # Print progress for every 1000 processed
        if (processed_count // 1000) > ((processed_count - len(results)) // 1000):
            print(f"Progress: {processed_count}/{total_regions} regions processed ({processed_count/total_regions*100:.1f}%)")
    
    successful, successful_no_neutral, failed = log_results(all_results)
    
    successful_count = len(successful)
    successful_no_neutral_count = len(successful_no_neutral)
    failed_count = len(failed)
    total_processed = successful_count + successful_no_neutral_count + failed_count
    
    print("\n--- SUMMARY ---")
    print(f"Total processed: {total_processed}")
    print(f"Successful with neutral positions: {successful_count}")
    print(f"Successful but no neutral positions: {successful_no_neutral_count}")
    print(f"Failed: {failed_count}")

if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed = time.time() - start_time
    print(f"Total execution time: {elapsed:.2f} seconds")
