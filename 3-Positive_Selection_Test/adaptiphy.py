import os
import re
import random
import subprocess
import pandas as pd
from Bio import SeqIO
from ete3 import Tree
from scipy import stats
import multiprocessing as mp
import time
import json
from collections import defaultdict
import glob
import tempfile
import logging
import math

# Configuration
QUERY_DIR = "../1-b-CREs_MSAs/Adaptify_MSAs" # Path of the folder containing b-CREs multiple sequence alignement
REF_DIR = "../2-Neutral_MSAs/References" # Path of the folder containing neutral multiple sequence alignement
MASTER_TREE_FILE = "243Primates_Tree.nwk"

# Define foreground branches
HUMAN_FOREGROUND = 'hg38'
GREAT_APES_CLADE = "((Pongo_pygmaeus,Pongo_abelii),((Gorilla_beringei,Gorilla_gorilla),(hg38,(Pan_troglodytes,Pan_paniscus))))"
GREAT_APE_SPECIES = {'Pongo_pygmaeus', 'Pongo_abelii', 'Gorilla_beringei', 'Gorilla_gorilla',
                     'hg38', 'Pan_troglodytes', 'Pan_paniscus'}

MODELS = ['null', 'alt']
NUM_PROCESSES = max(1, mp.cpu_count() - 1)

# Main directory (current working directory)
MAIN_DIR = os.getcwd()

# Set up logging in the main directory
logging.basicConfig(filename=os.path.join(MAIN_DIR, "script_log.txt"), level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Progress tracking
progress_interval = 5  # Percentage interval for logging progress

# Single JSON file for all results
RESULTS_FILE = os.path.join(MAIN_DIR, "processed_bCREs.json")

def read_fasta_species(fasta_file):
    """Extract species names from a FASTA file"""
    species = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        species.add(record.id)
    return species

def get_base_name(name):
    """Remove trailing identifiers from species names"""
    return re.sub(r'[#$]\d+$', '', name).strip()

def clean_newick_format(newick_string):
    """Remove branch lengths and semicolon from Newick string"""
    clean_newick = re.sub(r':\d+(\.\d+)?', '', newick_string)
    if clean_newick.endswith(';'):
        clean_newick = clean_newick[:-1]
    return clean_newick

def identify_great_apes_in_alignment(species_list):
    """Identify great ape species present in the alignment"""
    present_great_apes = set()
    for species in species_list:
        base_name = get_base_name(species)
        if base_name in GREAT_APE_SPECIES:
            present_great_apes.add(base_name)
    return present_great_apes

def prune_tree_for_alignment(master_tree_file, alignment_file, output_tree_file):
    """Prune master tree to match alignment species"""
    try:
        tree = Tree(master_tree_file, format=1)
        species = read_fasta_species(alignment_file)
        name_mapping = {get_base_name(leaf.name): leaf.name for leaf in tree.iter_leaves()}
        species_to_keep = {get_base_name(sp) for sp in species if get_base_name(sp) in name_mapping}
        
        if HUMAN_FOREGROUND not in species_to_keep:
            logging.warning(f"{HUMAN_FOREGROUND} not in alignment {alignment_file}")
            return None, None
        
        leaves_to_keep = [name_mapping[name] for name in species_to_keep]
        if not leaves_to_keep:
            logging.warning(f"No matching species for {alignment_file}")
            return None, None
        
        tree.prune(leaves_to_keep, preserve_branch_length=False)
        newick = tree.write(format=1)
        clean_newick = clean_newick_format(newick)
        
        with open(output_tree_file, 'w') as f:
            f.write(clean_newick)
        
        present_great_apes = identify_great_apes_in_alignment(species_to_keep)
        return clean_newick, present_great_apes
    except Exception as e:
        logging.error(f"Error pruning tree: {str(e)}")
        return None, None

def run_hyphy_analysis(query_file, ref_file, model, tree_string, foreground_label, analysis_type, temp_dir):
    """Run HyPhy analysis and return result file path"""
    base_name = os.path.basename(query_file).split('.')[0]
    abs_query_file = os.path.abspath(query_file)
    abs_ref_file = os.path.abspath(ref_file)
    abs_template_file = os.path.abspath(f"{model}4-fgrnd_spec.bf")
    res_file = os.path.join(temp_dir, "results_folder", f"{base_name}.{analysis_type}.{model}.res")
    
    if not os.path.exists(abs_template_file):
        logging.error(f"Template file missing: {abs_template_file}")
        return None
    
    batch_file = os.path.join(temp_dir, f"{base_name}.{analysis_type}.{model}.bf")
    log_file = os.path.join(temp_dir, f"{base_name}.{analysis_type}.{model}.log")
    
    try:
        with open(batch_file, 'w') as f:
            random_seed = random.randint(1, 1000)
            f.write(f"random_seed={random_seed};\n")
            f.write(f"quer_seq_file= \"{abs_query_file}\";\n")
            f.write(f"ref_seq_file = \"{abs_ref_file}\";\n")
            f.write("fit_repl_count = 20;\n")
            f.write(f"tree= \"{tree_string}\";\n")
            f.write(f"fgrnd_branch_name = \"{foreground_label}\";\n")
            f.write(f"res_file = \"{res_file}\";\n")
            f.write(f"#include \"{abs_template_file}\";\n")
    except Exception as e:
        logging.error(f"Error creating batch file: {e}")
        return None
    
    try:
        cmd = f"hyphy {batch_file}"
        with open(log_file, 'w') as f:
            subprocess.run(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
        os.remove(batch_file)
        os.remove(log_file)
    except Exception as e:
        logging.error(f"Error running HyPhy: {e}")
        return None
    
    return res_file if os.path.exists(res_file) else None

def hyphy_task_wrapper(args):
    """Wrapper for parallel HyPhy execution"""
    bCRE, query_file, ref_file, tree_string, foreground_label, analysis_type, model, temp_dir = args
    res_file = run_hyphy_analysis(query_file, ref_file, model, tree_string, foreground_label, analysis_type, temp_dir)
    return {'bCRE': bCRE, 'analysis_type': analysis_type, 'model': model, 'res_file': res_file}

def extract_logl(res_file):
    """Extract log-likelihood from result file"""
    if not res_file or not os.path.exists(res_file):
        return None
    with open(res_file, 'r') as f:
        content = f.read()
        match = re.search(r"BEST LOG-L:\s*([-\d.]+)", content)
        if match:
            return float(match.group(1))
    return None

def calculate_lrt(null_logl, alt_logl):
    """Calculate LRT statistic"""
    if null_logl is None or alt_logl is None:
        return None
    return 2 * (alt_logl - null_logl)

def calculate_p_value(lrt_value):
    """Calculate p-value from LRT"""
    if lrt_value is None:
        return None
    return 1 - stats.chi2.cdf(lrt_value, df=1)

def generate_tasks(bCRE, query_file, ref_file, tree_string, present_great_apes, temp_dir):
    """Generate HyPhy tasks for a bCRE"""
    tasks = []
    tasks.append((bCRE, query_file, ref_file, tree_string, HUMAN_FOREGROUND, 'human', 'null', temp_dir))
    tasks.append((bCRE, query_file, ref_file, tree_string, HUMAN_FOREGROUND, 'human', 'alt', temp_dir))
    if len(present_great_apes) >= 2:
        tasks.append((bCRE, query_file, ref_file, tree_string, 'great_apes', 'great_apes', 'null', temp_dir))
        tasks.append((bCRE, query_file, ref_file, tree_string, 'great_apes', 'great_apes', 'alt', temp_dir))
    return tasks

def load_processed_bCREs():
    """Load processed bCREs and their results from the single JSON file"""
    if os.path.exists(RESULTS_FILE):
        try:
            with open(RESULTS_FILE, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError:
            logging.warning(f"Error loading {RESULTS_FILE}, creating new")
    return {"processed_ids": [], "results": []}

def save_processed_bCREs(data):
    """Save processed bCREs and their results to the single JSON file"""
    with open(RESULTS_FILE, 'w') as f:
        json.dump(data, f, indent=2)

def main():
    start_time = time.time()
    
    # Load existing processed_bCREs from the single JSON file
    processed_data = load_processed_bCREs()
    processed_bCREs = processed_data["processed_ids"]
    
    logging.info(f"Loaded {len(processed_bCREs)} already processed bCREs")
    
    query_files = [os.path.join(QUERY_DIR, f) for f in os.listdir(QUERY_DIR) if f.endswith('.fa')]
    if not query_files:
        logging.info(f"No query files in {QUERY_DIR}")
        return
    
    logging.info(f"Found {len(query_files)} query files")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Set directories inside temporary directory
        OUTPUT_DIR = os.path.join(temp_dir, "results_folder")
        HYPHY_DIR = os.path.join(temp_dir, "hyphy_output")
        PRUNED_TREE_DIR = os.path.join(temp_dir, "pruned_trees")
        for directory in [OUTPUT_DIR, HYPHY_DIR, PRUNED_TREE_DIR]:
            os.makedirs(directory, exist_ok=True)
        
        all_tasks = []
        bCRE_data = {}
        
        for query_file in query_files:
            bCRE = os.path.basename(query_file).split('.')[0]
            if bCRE in processed_bCREs:
                logging.info(f"Skipping processed bCRE: {bCRE}")
                continue
            
            ref_file = os.path.join(REF_DIR, f"{bCRE}.fa")
            if not os.path.exists(ref_file):
                logging.warning(f"Reference file missing for {bCRE}")
                continue
            
            pruned_tree_file = os.path.join(PRUNED_TREE_DIR, f"{bCRE}.nwk")
            tree_string, present_great_apes = prune_tree_for_alignment(MASTER_TREE_FILE, query_file, pruned_tree_file)
            if tree_string is None:
                logging.warning(f"Failed to prune tree for {bCRE}")
                continue
            
            bCRE_data[bCRE] = {'present_great_apes': present_great_apes}
            tasks = generate_tasks(bCRE, query_file, ref_file, tree_string, present_great_apes, temp_dir)
            all_tasks.extend(tasks)
        
        if all_tasks:
            total_tasks = len(all_tasks)
            logging.info(f"Running {total_tasks} HyPhy tasks across {len(bCRE_data)} bCREs with {NUM_PROCESSES} processes")
            
            # Calculate the number of tasks for each 5% interval
            interval_tasks = math.ceil(total_tasks * (progress_interval / 100))
            completed_tasks = 0
            last_logged_percent = 0
            
            results = []
            with mp.Pool(NUM_PROCESSES) as pool:
                for result in pool.imap_unordered(hyphy_task_wrapper, all_tasks):
                    results.append(result)
                    completed_tasks += 1
                    current_percent = math.floor((completed_tasks / total_tasks) * 100)
                    if current_percent >= last_logged_percent + progress_interval:
                        elapsed_time = time.time() - start_time
                        logging.info(f"Progress: {current_percent}% completed. Elapsed time: {elapsed_time:.2f} seconds")
                        last_logged_percent = current_percent - (current_percent % progress_interval)
            
            results_by_bCRE = defaultdict(list)
            for result in results:
                if result['res_file']:
                    results_by_bCRE[result['bCRE']].append(result)
            
            # Track newly processed bCREs
            new_results = []
            
            for bCRE, bCRE_results in results_by_bCRE.items():
                present_great_apes = bCRE_data[bCRE]['present_great_apes']
                
                human_null_res = next((r['res_file'] for r in bCRE_results if r['analysis_type'] == 'human' and r['model'] == 'null'), None)
                human_alt_res = next((r['res_file'] for r in bCRE_results if r['analysis_type'] == 'human' and r['model'] == 'alt'), None)
                great_ape_null_res = next((r['res_file'] for r in bCRE_results if r['analysis_type'] == 'great_apes' and r['model'] == 'null'), None)
                great_ape_alt_res = next((r['res_file'] for r in bCRE_results if r['analysis_type'] == 'great_apes' and r['model'] == 'alt'), None)
                
                human_null_logl = extract_logl(human_null_res)
                human_alt_logl = extract_logl(human_alt_res)
                great_ape_null_logl = extract_logl(great_ape_null_res)
                great_ape_alt_logl = extract_logl(great_ape_alt_res)
                
                human_lrt = calculate_lrt(human_null_logl, human_alt_logl)
                human_p_value = calculate_p_value(human_lrt)
                great_ape_lrt = calculate_lrt(great_ape_null_logl, great_ape_alt_logl)
                great_ape_p_value = calculate_p_value(great_ape_lrt)
                
                result = {
                    'bCRE': bCRE,
                    'human_branch': HUMAN_FOREGROUND,
                    'human_null_lnL': human_null_logl,
                    'human_alt_lnL': human_alt_logl,
                    'human_lrt': human_lrt,
                    'human_p_value': human_p_value,
                    'human_significant': 'Yes' if human_p_value and human_p_value < 0.05 else 'No',
                    'great_ape_branch': 'great_apes',
                    'great_ape_species_included': ", ".join(present_great_apes) if len(present_great_apes) >= 2 else None,
                    'great_ape_null_lnL': great_ape_null_logl,
                    'great_ape_alt_lnL': great_ape_alt_logl,
                    'great_ape_lrt': great_ape_lrt,
                    'great_ape_p_value': great_ape_p_value,
                    'great_ape_significant': 'Yes' if great_ape_p_value and great_ape_p_value < 0.05 else 'No'
                }
                
                # Add to new results
                new_results.append(result)
                
                # Add to processed_bCREs list
                processed_bCREs.append(bCRE)
                logging.info(f"Processed bCRE: {bCRE}")
            
            # Update the processed data and save
            if new_results:
                processed_data["results"].extend(new_results)
                processed_data["processed_ids"] = processed_bCREs
                save_processed_bCREs(processed_data)
                logging.info(f"Updated {RESULTS_FILE} with {len(new_results)} new results")
        
        # Create Excel output from the single JSON file
        if processed_data["results"]:
            df = pd.DataFrame(processed_data["results"])
            ordered_columns = [
                'bCRE', 'human_branch', 'human_null_lnL', 'human_alt_lnL', 'human_lrt', 'human_p_value', 'human_significant',
                'great_ape_branch', 'great_ape_species_included', 'great_ape_null_lnL', 'great_ape_alt_lnL', 'great_ape_lrt',
                'great_ape_p_value', 'great_ape_significant'
            ]
            df = df[ordered_columns]
            df.to_excel(os.path.join(MAIN_DIR, "evolutionary_analysis_results.xlsx"), sheet_name='Results', index=False)
            logging.info(f"Results saved to evolutionary_analysis_results.xlsx")
        else:
            logging.info("No results available for Excel output.")
    
    logging.info(f"Execution time: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
