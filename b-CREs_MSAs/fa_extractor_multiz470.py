#!/usr/bin/env python3
"""
Script to download MAF data from UCSC API for regions in a BED file and convert to FASTA.
Optimized for speed while respecting UCSC guidelines.
Saves only FASTA files to minimize disk usage.
Ensures exact BED coordinates are extracted from MAF alignments.
Renames sequences in FASTA files to use scientific names.
"""

import requests
import concurrent.futures
import time
import os
import subprocess
import tempfile
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import re
import json
# Configuration - EDIT THESE SETTINGS
BED_FILE = "Intermediate_Files/Cactus_447_regions_below_60.bed"  # Path to your BED file
OUTPUT_DIR = "Intermediate_Files/Multiz470_Raw_MSAs"                # Directory to save FASTA files
GENOME = "hg38"                     # Genome assembly
TRACK = "multiz470way"             # Alignment track
MAX_WORKERS = min(32, os.cpu_count() * 4)  # Optimize for I/O bound tasks
MAX_RETRIES = 5                    # Number of retry attempts
TIMEOUT = 30                       # Request timeout in seconds
MIN_REQUEST_INTERVAL = 0.15        # Minimum time between requests per thread (seconds)
SUCCESS_LOG = os.path.join(OUTPUT_DIR, "successful_downloads.log")
ERROR_LOG = os.path.join(OUTPUT_DIR, "failed_downloads.log")
REFERENCE_SPECIES = "hg38"         # Reference species in the MAF (usually human)
with open('config.json') as f:
    config = json.load(f)
    MSA_VIEW_PATH = config['msa_view_path']
# Dictionary mapping UCSC names to scientific names
SPECIES_MAPPING = {
    "panPan3": "Pan_paniscus",
    "panTro6": "Pan_troglodytes",
    "gorGor6": "Gorilla_gorilla",
    "ponAbe3": "Pongo_abelii",
    "HLhylMol2": "Hylobates_muelleri",
    "HLnomLeu4": "Nomascus_leucogenys",
    "rheMac10": "Macaca_mulatta",
    "HLmacFus1": "Macaca_fuscata",
    "HLmacFas6": "Macaca_fascicularis",
    "macNem1": "Macaca_nemestrina",
    "HLpapAnu5": "Papio_anubis",
    "HLtheGel1": "Theropithecus_gelada",
    "cerAty1": "Cercocebus_atys",
    "HLmanSph1": "Mandrillus_sphinx",
    "manLeu1": "Mandrillus_leucophaeus",
    "HLcerMon1": "Cercopithecus_mona",
    "HLcerNeg1": "Cercopithecus_neglectus",
    "HLallNig1": "Allenopithecus_nigroviridis",
    "HLeryPat1": "Erythrocebus_patas",
    "chlSab2": "Chlorocebus_sabaeus",
    "HLpilTep2": "Piliocolobus_tephrosceles",
    "colAng1": "Colobus_angolensis",
    "HLsemEnt1": "Semnopithecus_entellus",
    "HLtraFra1": "Trachypithecus_francoisi",
    "nasLar1": "Nasalis_larvatus",
    "HLpygNem1": "Pygathrix_nemaeus",
    "HLrhiRox2": "Rhinopithecus_roxellana",
    "rhiBie1": "Rhinopithecus_bieti",
    "HLpitPit1": "Pithecia_pithecia",
    "HLateGeo1": "Ateles_geoffroyi",
    "HLaloPal1": "Alouatta_palliata",
    "saiBol1": "Saimiri_boliviensis",
    "HLsapApe1": "Sapajus_apella",
    "HLcebAlb1": "Cebus_albifrons",
    "aotNan1": "Aotus_nancymaae",
    "HLsagImp1": "Saguinus_imperator",
    "HLcalJac4": "Callithrix_jacchus",
    "HLcalPym1": "Cebuella_pygmaea",
    "tarSyr2": "Carlito_syrichta",
    "HLdauMad1": "Daubentonia_madagascariensis",
    "HLindInd1": "Indri_indri",
    "proCoq1": "Propithecus_coquerelli",
    "HLcheMed1": "Cheirogaleus_medius",
    "HLmirZaz1": "Mirza_zaza",
    "micMur3": "Microcebus_murinus",
    "HLlemCat1": "Lemur_catta",
    "HLproSim1": "Prolemur_simus",
    "HLeulFul1": "Eulemur_fulvus",
    "HLeulMon1": "Eulemur_mongoz",
    "eulFla1": "Eulemur_flavifrons",
    "eulMac1": "Eulemur_macaco",
    "HLnycCou1": "Nycticebus_coucang",
    "otoGar3": "Otolemur_garnettii",
}

# Read primates.txt for sequence list
try:
    with open("primates_multiz.txt", 'r') as f:
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

def rename_fasta_sequences(fasta_path):
    """
    Rename sequences in FASTA file to use scientific names.
    Returns True if successful, False otherwise.
    """
    try:
        # Read the original FASTA file
        with open(fasta_path, 'r') as f:
            content = f.read()
        
        # Create a temporary file for the new content
        temp_fasta_path = fasta_path + ".temp"
        renamed = False
        
        with open(temp_fasta_path, 'w') as f_temp:
            # Process each line
            for line in content.splitlines():
                if line.startswith('>'):
                    # This is a header line
                    new_line = line
                    
                    # Replace each UCSC species ID with its scientific name
                    for ucsc_id, scientific_name in SPECIES_MAPPING.items():
                        # Look for the species ID as a whole word in the header
                        # This handles various header formats like ">panPan3" or ">panPan3.chr1"
                        pattern = r'(^>|\s+)' + re.escape(ucsc_id) + r'(\.|$|\s+)'
                        match = re.search(pattern, new_line)
                        if match:
                            # Replace only the species ID part, preserving the rest
                            start, end = match.span(0)
                            prefix = new_line[start:start+match.group(1).__len__()]
                            suffix = new_line[end-match.group(2).__len__():end]
                            # Perform the replacement
                            new_line = new_line[:start] + prefix + scientific_name + suffix + new_line[end:]
                            renamed = True
                            # In case multiple replacements are needed, break after the first match
                            break
                    
                    f_temp.write(new_line + '\n')
                else:
                    # This is a sequence line, keep it as is
                    f_temp.write(line + '\n')
        
        # If we renamed at least one sequence, replace the original file with our new one
        if renamed:
            os.replace(temp_fasta_path, fasta_path)
        else:
            # No renaming was done, remove the temp file
            os.remove(temp_fasta_path)
        
        return renamed
    except Exception as e:
        print(f"Error in rename_fasta_sequences: {str(e)}")
        # Clean up temp file if it exists
        if os.path.exists(temp_fasta_path):
            try:
                os.remove(temp_fasta_path)
            except:
                pass
        return False
def download_maf(session, chrom, start, end, position):
    """Download MAF data, convert to FASTA, and save FASTA file with coordinate trimming."""
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
            msa_view_cmd = [MSA_VIEW_PATH, temp_maf_path, "--out-format", "FASTA", "--seqs", SEQ_LIST]
            
            # Add trim coordinates if available
            detail_msg = ""
            if trim_coords:
                msa_view_cmd.extend(["--start", str(trim_coords['start']), "--end", str(trim_coords['end'])])
                
                # For detailed logging, include both the MAF and BED coordinates
                maf_start = parsed_blocks[0]['start'] if parsed_blocks else "unknown"
                detail_msg = f" (MAF starts at {maf_start}, trimmed to columns {trim_coords['start']}-{trim_coords['end']})"
            else:
                detail_msg = " (untrimmed - using full MAF coordinates)"
            
            # Debug information to log
            debug_info = f"BED:{start}-{end}, "
            if parsed_blocks:
                ref_info = parsed_blocks[0]
                debug_info += f"MAF:{ref_info['start']}, SeqLen:{len(ref_info['seq'])}"
            else:
                debug_info += "No block info"
            
            result = subprocess.run(msa_view_cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                error_msg = f"msa_view failed: {result.stderr}"
                with open(ERROR_LOG, 'a') as f:
                    f.write(f"{position}: {error_msg} | {debug_info}\n")
                return {"status": "failed", "position": position, "error": error_msg}
            
            fasta_content = result.stdout
            if not fasta_content.strip():
                error_msg = "msa_view produced empty FASTA"
                with open(ERROR_LOG, 'a') as f:
                    f.write(f"{position}: {error_msg} | {debug_info}\n")
                return {"status": "failed", "position": position, "error": error_msg}
            
            # Write the FASTA content to file
            with open(output_file, 'w') as f:
                f.write(fasta_content)
            
            # Now rename the sequences in the FASTA file
            renamed = rename_fasta_sequences(output_file)
            
            with open(SUCCESS_LOG, 'a') as f:
                rename_msg = " (sequences renamed to scientific names)" if renamed else ""
                f.write(f"{position} (blocks: {blocks_found}){detail_msg}{rename_msg} | {debug_info}\n")
            
            return {
                "status": "success", 
                "position": position, 
                "file": output_file, 
                "blocks": blocks_found,
                "trimmed": trim_coords is not None,
                "renamed": renamed,
                "debug_info": debug_info
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
    
    return download_maf(session, chrom, start, end, position)

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

def log_results(successful, failed):
    """Write successful and failed downloads to log files."""
    with open(SUCCESS_LOG, 'w') as f:
        for position, blocks, trimmed, renamed, debug_info in successful:
            trim_status = "trimmed" if trimmed else "untrimmed"
            rename_status = "renamed" if renamed else "not renamed"
            f.write(f"{position} (blocks: {blocks}, {trim_status}, {rename_status}) | {debug_info}\n")
            
    with open(ERROR_LOG, 'w') as f:
        for position, error in failed:
            f.write(f"{position}: {error}\n")

def main():
    try:
        with open(BED_FILE, 'r') as f:
            all_lines = [line for line in f.readlines() 
                         if not line.startswith('#') and line.strip()]
    except IOError as e:
        print(f"Error reading BED file: {e}")
        return
    
    successful_positions = []
    failed_positions = []
    total_blocks = 0
    trimmed_count = 0
    renamed_count = 0
    total_regions = len(all_lines)
    
    print(f"Processing {total_regions} regions...")
    
    processed_count = 0
    batch_size = min(100, len(all_lines))
    
    for i in range(0, len(all_lines), batch_size):
        batch = all_lines[i:i+batch_size]
        results = process_batch(batch)
        
        for result in results:
            if result["status"] == "success":
                blocks = result.get("blocks", 0)
                trimmed = result.get("trimmed", False)
                renamed = result.get("renamed", False)
                debug_info = result.get("debug_info", "No debug info")
                if trimmed:
                    trimmed_count += 1
                if renamed:
                    renamed_count += 1
                total_blocks += blocks
                successful_positions.append((result["position"], blocks, trimmed, renamed, debug_info))
            else:
                failed_positions.append((result["position"], result["error"]))
            
            processed_count += 1
            
        # Print progress for every 1000 processed
        if (processed_count // 1000) > ((processed_count - len(results)) // 1000):
            print(f"Progress: {processed_count}/{total_regions} regions processed ({processed_count/total_regions*100:.1f}%)")
    
    log_results(successful_positions, failed_positions)
    
    successful = len(successful_positions)
    failed = len(failed_positions)
    
    print("\n--- SUMMARY ---")
    print(f"Total processed: {successful + failed}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Regions that required trimming: {trimmed_count} ({trimmed_count/successful*100:.1f}% of successful)")
    print(f"Regions with renamed sequences: {renamed_count} ({renamed_count/successful*100:.1f}% of successful)")
    
    if successful > 0:
        avg_blocks = total_blocks / successful
        print(f"Average alignment blocks per region: {avg_blocks:.2f}")

def rename_existing_fasta_files():
    """
    Utility function to rename sequences in all existing FASTA files in the output directory.
    Can be called separately if you already have downloaded FASTA files.
    """
    print(f"Renaming sequences in existing FASTA files in {OUTPUT_DIR}...")
    count = 0
    renamed_count = 0
    
    for filename in os.listdir(OUTPUT_DIR):
        if filename.endswith('.fa'):
            count += 1
            filepath = os.path.join(OUTPUT_DIR, filename)
            if rename_fasta_sequences(filepath):
                renamed_count += 1
                
            # Print progress for every 100 files
            if count % 100 == 0:
                print(f"Progress: {count} files processed, {renamed_count} renamed")
    
    print(f"Completed renaming. Processed {count} files, renamed sequences in {renamed_count} files.")

if __name__ == "__main__":
    # Choose one of these options:
    
    # Option 1: Run the full download and processing (including renaming)
    start_time = time.time()
    main()
    elapsed = time.time() - start_time
    print(f"Total execution time: {elapsed:.2f} seconds")
    
    # Option 2: If you already have FASTA files and just want to rename them
    # Uncomment the following line and comment out the main() call above
    # rename_existing_fasta_files()
