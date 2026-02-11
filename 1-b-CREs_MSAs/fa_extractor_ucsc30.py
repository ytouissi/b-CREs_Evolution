#!/usr/bin/env python3
"""
Script to download MAF data from UCSC and extract Multiple Sequence Alignment (MSA)
in FASTA format with species filtering and renaming.
"""

import requests
import concurrent.futures
import time
import os
import subprocess
import tempfile
import sys
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import json
# Configuration - EDIT THESE SETTINGS
# Get your UCSC API key from: https://genome.ucsc.edu/cgi-bin/hgHubConnect#hubDeveloper
# UCSC Browser: Track Hubs > Hub Developement > Generate  Key (at bottom of Hub Development page)
BED_FILE = "Intermediate_Files/Multiz470_regions_below_27.bed"           # Path to your BED file
OUTPUT_DIR = "Intermediate_Files/UCSC30_Raw_MSAs"                     # Directory to save FASTA files
MAX_WORKERS = min(32, os.cpu_count() * 4) # Optimize for I/O bound tasks
MAX_RETRIES = 5                           # Number of retry attempts
TIMEOUT = 30                              # Request timeout in seconds
MIN_REQUEST_INTERVAL = 0.15               # Minimum time between requests per thread (seconds)
SUCCESS_LOG = os.path.join(OUTPUT_DIR, "successful_downloads.log")
ERROR_LOG = os.path.join(OUTPUT_DIR, "failed_downloads.log")
with open('config.json') as f:
    config = json.load(f)
    MSA_VIEW_PATH = config['msa_view_path']
    UCSC_API_KEY = config['api_key']
# Species mapping for filtering and renaming
SPECIES_MAPPING = {
    "hg38": "hg38",
    "panPan2": "Pan_paniscus",
    "panTro5": "Pan_troglodytes",
    "gorGor5": "Gorilla_gorilla",
    "ponAbe2": "Pongo_abelii",
    "nomLeu3": "Nomascus_leucogenys",
    "rheMac8": "Macaca_mulatta",
    "macFas5": "Macaca_fascicularis",
    "macNem1": "Macaca_nemestrina",
    "cerAty1": "Cercocebus_atys",
    "papAnu3": "Papio_anubis",
    "chlSab2": "Chlorocebus_sabaeus",
    "manLeu1": "Mandrillus_leucophaeus",
    "nasLar1": "Nasalis_larvatus",
    "colAng1": "Colobus_angolensis",
    "rhiRox1": "Rhinopithecus_roxellana",
    "rhiBie1": "Rhinopithecus_bieti",
    "calJac3": "Callithrix_jacchus",
    "saiBol1": "Saimiri_boliviensis",
    "aotNan1": "Aotus_nancymaae",
    "tarSyr2": "Carlito_syrichta",
    "micMur3": "Microcebus_murinus",
    "proCoq1": "Propithecus_coquereli",
    "eulMac1": "Eulemur_macaco",
    "eulFla1": "Eulemur_flavifrons",
    "otoGar3": "Otolemur_garnettii",
}

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

def clean_maf_content(maf_content):
    """Remove error lines from MAF that break msa_view."""
    lines = maf_content.split('\n')
    cleaned_lines = []
    
    for line in lines:
        # Skip lines that are just dashes (error markers)
        if line.strip() and line.strip().replace('-', '') == '':
            continue
        # Skip error messages from UCSC
        if 'Bad start cookie' in line or 'freeing' in line:
            continue
        cleaned_lines.append(line)
    
    return '\n'.join(cleaned_lines)

def create_session():
    """Create and configure a requests session with retry logic and connection pooling."""
    session = requests.Session()
    retry_strategy = Retry(
        total=MAX_RETRIES,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"]
    )
    adapter = HTTPAdapter(
        max_retries=retry_strategy,
        pool_connections=MAX_WORKERS,
        pool_maxsize=MAX_WORKERS
    )
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

def download_maf(session, chrom, start, end, position):
    """Download MAF data and save as FASTA file."""
    url = "https://genome.ucsc.edu/cgi-bin/hgTables"
    output_file = os.path.join(OUTPUT_DIR, f"{chrom}_{start}_{end}.fa")
    
    # Complete form parameters for Table Browser to get MAF
    payload = {
        'hgta_track': 'cons30way',
        'hgta_table': 'multiz30way',
        'hgta_regionType': 'range',
        'position': f'{chrom}:{start}-{end}',
        'hgta_outputType': 'maf',
        'hgta_doTopSubmit': 'get output',
        'apiKey': UCSC_API_KEY
    }
    
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded'
    }
    
    try:
        # Use POST request with complete form data
        response = session.post(url, data=payload, headers=headers, timeout=TIMEOUT)
        response.raise_for_status()
        maf_content = response.text
        
        # Check if we got valid MAF content
        if "##maf" not in maf_content and "a score=" not in maf_content:
            with open(ERROR_LOG, 'a') as f:
                f.write(f"{position}: Invalid MAF response received\n")
            return {"status": "failed", "position": position, "error": "Invalid MAF response"}
        
        # Count alignment blocks
        blocks_found = maf_content.count('a score=')
        
        if blocks_found == 0:
            with open(ERROR_LOG, 'a') as f:
                f.write(f"{position}: No alignment blocks found\n")
            return {"status": "failed", "position": position, "error": "No alignment blocks found"}
        
        # Clean MAF content to remove error lines that break msa_view
        maf_content = clean_maf_content(maf_content)
        
        # Save the MAF content to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.maf', delete=False) as temp_maf:
            temp_maf.write(maf_content)
            temp_maf_path = temp_maf.name
        
        try:
            # Use msa_view to convert MAF to FASTA
            try:
                msa_view_cmd = [MSA_VIEW_PATH, temp_maf_path, "--out-format", "FASTA"]
                result = subprocess.run(msa_view_cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    raise Exception("msa_view failed")
                
                fasta_content = result.stdout
                
            except:
                # If msa_view fails, use Python-based extraction
                fasta_content = extract_fasta_from_maf(maf_content)
            
            # Check if we got valid FASTA content
            if not fasta_content.strip():
                fasta_content = extract_fasta_from_maf(maf_content)
                if not fasta_content.strip():
                    with open(ERROR_LOG, 'a') as f:
                        f.write(f"{position}: Empty FASTA output\n")
                    return {"status": "failed", "position": position, "error": "Empty FASTA output"}
            
            # Apply species filtering and renaming
            filtered_fasta = filter_and_rename_species(fasta_content)
            
            # Save the filtered FASTA file
            with open(output_file, 'w') as f:
                f.write(filtered_fasta)
            
            # Count species/sequences in filtered output
            species_count = filtered_fasta.count('>')
            
            with open(SUCCESS_LOG, 'a') as f:
                f.write(f"{position} (blocks: {blocks_found}, species: {species_count})\n")
            
            return {
                "status": "success", 
                "position": position, 
                "file": output_file, 
                "blocks": blocks_found,
                "species": species_count
            }
            
        finally:
            # Clean up temporary file
            os.unlink(temp_maf_path)
            
    except Exception as e:
        with open(ERROR_LOG, 'a') as f:
            f.write(f"{position}: {str(e)}\n")
        return {"status": "failed", "position": position, "error": str(e)}

def filter_and_rename_species(fasta_content):
    """Filter and rename species according to SPECIES_MAPPING."""
    filtered_lines = []
    current_species = None
    
    lines = fasta_content.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('>'):
            # Extract species identifier (removing '>' and any additional identifiers after space or dot)
            species_id = line[1:].split('.')[0].split()[0]
            
            # Check if this species is in our mapping
            if species_id in SPECIES_MAPPING:
                current_species = species_id
                # Add header with renamed species
                filtered_lines.append(f">{SPECIES_MAPPING[species_id]}")
            else:
                current_species = None
        elif current_species is not None:
            # Add sequence lines for the current species if it's in our mapping
            filtered_lines.append(line)
        
        i += 1
    
    return "\n".join(filtered_lines)

def extract_fasta_from_maf(maf_content):
    """Extract FASTA from MAF as a fallback if msa_view is not available."""
    fasta_content = []
    alignment_species = {}
    current_block = []
    in_block = False
    
    for line in maf_content.splitlines():
        line = line.strip()
        if not line:
            continue
            
        # New alignment block
        if line.startswith('a '):
            if current_block:
                # Process the previous block
                for seq_line in current_block:
                    if seq_line.startswith('s '):
                        parts = seq_line.split()
                        if len(parts) >= 7:
                            species = parts[1].split('.')[0]  # Get species from sequence name
                            sequence = parts[6]  # The sequence is in the 7th field
                            
                            # Add or append to species dictionary
                            if species in alignment_species:
                                alignment_species[species] += sequence
                            else:
                                alignment_species[species] = sequence
            
            # Start a new block
            current_block = []
            in_block = True
        elif in_block:
            current_block.append(line)
    
    # Process the last block if there is one
    if current_block:
        for seq_line in current_block:
            if seq_line.startswith('s '):
                parts = seq_line.split()
                if len(parts) >= 7:
                    species = parts[1].split('.')[0]
                    sequence = parts[6]
                    
                    if species in alignment_species:
                        alignment_species[species] += sequence
                    else:
                        alignment_species[species] = sequence
    
    # Convert alignment_species dictionary to FASTA format
    # Only include species in our mapping and use the mapped names
    for species, sequence in alignment_species.items():
        if species in SPECIES_MAPPING:
            fasta_content.append(f">{SPECIES_MAPPING[species]}")
            fasta_content.append(sequence)
    
    return "\n".join(fasta_content)

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
    
    # Add a small delay to avoid overwhelming the server
    time.sleep(MIN_REQUEST_INTERVAL)
    
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

def main():
    # Print species mapping info
    print(f"Species filtering enabled. Only keeping {len(SPECIES_MAPPING)} species with custom naming.")
    
    # Single region mode
    if len(sys.argv) > 1 and ':' in sys.argv[1] and '-' in sys.argv[1]:
        region = sys.argv[1]
        parts = region.split(':')
        chrom = parts[0]
        start, end = parts[1].split('-')
        
        print(f"Processing single region: {region}")
        session = create_session()
        result = download_maf(session, chrom, start, end, region)
        
        if result["status"] == "success":
            print(f"Success! FASTA data saved to: {result['file']}")
            print(f"Alignment blocks: {result['blocks']}")
            print(f"Species count: {result['species']}")
        else:
            print(f"Failed: {result['error']}")
        
        return
    
    # BED file mode
    bed_file = BED_FILE
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        bed_file = sys.argv[1]
    
    if not os.path.exists(bed_file):
        print(f"Error: BED file '{bed_file}' not found.")
        print("Usage options:")
        print(f"  {sys.argv[0]} chr11:71787967-71788313")
        print(f"  {sys.argv[0]} regions.bed")
        return
    
    try:
        with open(bed_file, 'r') as f:
            all_lines = [line for line in f.readlines() 
                         if not line.startswith('#') and line.strip()]
    except IOError as e:
        print(f"Error reading BED file: {e}")
        return
    
    successful_positions = []
    failed_positions = []
    total_blocks = 0
    total_species = 0
    total_regions = len(all_lines)
    
    print(f"Processing {total_regions} regions from {bed_file}...")
    
    processed_count = 0
    batch_size = min(100, len(all_lines))
    
    for i in range(0, len(all_lines), batch_size):
        batch = all_lines[i:i+batch_size]
        results = process_batch(batch)
        
        for result in results:
            if result["status"] == "success":
                blocks = result.get("blocks", 0)
                species = result.get("species", 0)
                total_blocks += blocks
                total_species += species
                successful_positions.append((result["position"], blocks, species))
            else:
                failed_positions.append((result["position"], result["error"]))
            
            processed_count += 1
            
        # Print progress
        print(f"Progress: {processed_count}/{total_regions} regions processed ({processed_count/total_regions*100:.1f}%)")
    
    # Write final logs
    with open(SUCCESS_LOG, 'w') as f:
        for position, blocks, species in successful_positions:
            f.write(f"{position} (blocks: {blocks}, species: {species})\n")
            
    with open(ERROR_LOG, 'w') as f:
        for position, error in failed_positions:
            f.write(f"{position}: {error}\n")
    
    successful = len(successful_positions)
    failed = len(failed_positions)
    
    print("\n--- SUMMARY ---")
    print(f"Total processed: {successful + failed}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    
    if successful > 0:
        avg_blocks = total_blocks / successful
        avg_species = total_species / successful
        print(f"Average alignment blocks per region: {avg_blocks:.2f}")
        print(f"Average species per region: {avg_species:.2f}")
        print(f"FASTA files saved to: {OUTPUT_DIR}")
        print(f"Species filtered and renamed according to mapping ({len(SPECIES_MAPPING)} species)")

if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed = time.time() - start_time
    print(f"Total execution time: {elapsed:.2f} seconds")