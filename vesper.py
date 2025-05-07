import os
import pandas as pd
import glob
import gzip
import cstag
import re
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import pysam
import csv
from itertools import combinations
from collections import Counter
import subprocess
import warnings
import matplotlib.pyplot as plt
import time
import numpy as np
import yaml
import argparse
import shutil
import math

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

### Parse command line arguments ###
parser = argparse.ArgumentParser(description="Run vesper with specified config.")
parser.add_argument("-c", "--config", required=True, help="Path to vesper_config.yaml file")
args = parser.parse_args()

### Load the YAML config file ###
config_file = os.path.abspath(args.config)
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

### Assign variables from VEP_config.yaml file ###
max_workers = config["max_workers"]
cx_sam_dir = config["input_sam_dir"]
root_output_dir = config["root_output_dir"]
cDNA_file = config["cDNA_file"]
vep_path = config["vep_path"]
start_bp = config["start_bp"]
end_bp = config["end_bp"]
baseq_threshold = config["baseq_threshold"]
mapq_threshold = config["mapq_threshold"]
seq_error = config["seq_error"]
sample_names = config['sample_names']
read_length_filter = config["read_length_filter"]
ex = config["exon_number"]
transcript = config["transcript"]
sample_combinations = config["sample_combinations"]

### Define output directories ###
unzipped_sam_dir = os.path.join(root_output_dir, f'Ex{ex}_unzipped_sam')                                    ### to be deleted ###
filtered_sam_unzipped_dir = os.path.join(root_output_dir, f'Ex{ex}_unzipped_filtered_sam')                  ### to be deleted ###
filtered_sam_gzipped_dir = os.path.join(root_output_dir, f'Ex{ex}_filtered_sam_gz')                         ### final filtered SAM files to be processed further ###

intermediate_results_output_dir = os.path.join(root_output_dir, 'intermediate_files')
main_results_output_dir = os.path.join(root_output_dir, 'main_results') 

raw_output_dir = os.path.join(intermediate_results_output_dir, 'sample_rep_raw')                             ### to be kept in intermediate ###
variants_per_rep_dir = os.path.join(main_results_output_dir, 'variants_per_rep')                             ### to be kept in main_results ###
nm_histograms_dir = os.path.join(main_results_output_dir, 'nm_histograms')                                   ### to be kept in main_results ###
variants_per_sample_dir = os.path.join(main_results_output_dir, 'variants_per_sample')                       ### to be kept in main_results ###
HGVS_dir = os.path.join(intermediate_results_output_dir, 'HGVS')                                             ### to be kept in intermediate ###
vep_txt_dir = os.path.join(intermediate_results_output_dir, 'VEP', 'txt_files')                              ### to be kept in intermediate ###
vep_csv_dir = os.path.join(intermediate_results_output_dir, 'VEP', 'csv_files')                              ### to be kept in intermediate ###
vep_maf_dir = os.path.join(main_results_output_dir, 'VEP_MAF')                                               ### to be kept in main_results ###
variant_consequences_dir = os.path.join(main_results_output_dir, 'sample_variant_consequences')              ### to be kept in main_results ###
sample_variant_comparisons_dir = os.path.join(main_results_output_dir, 'sample_variant_comparisons')         ### to be kept in main_results ###
sample_SNV_counts_dir = os.path.join(main_results_output_dir, 'sample_SNV_counts')                           ### to be kept in main_results ###
sample_SNV_pairs_dir = os.path.join(main_results_output_dir, 'sample_SNV_pairs')                             ### to be kept in main_results ###

### Create directories ###
for directory in [unzipped_sam_dir, filtered_sam_unzipped_dir, filtered_sam_gzipped_dir, intermediate_results_output_dir, main_results_output_dir, raw_output_dir, variants_per_rep_dir, nm_histograms_dir, variants_per_sample_dir,
    HGVS_dir, vep_txt_dir, vep_csv_dir, vep_maf_dir, variant_consequences_dir, sample_variant_comparisons_dir, sample_SNV_counts_dir, sample_SNV_pairs_dir]:
    os.makedirs(directory, exist_ok=True)


####################################################################################################################################################################################################################
                                                                ### FILTER SAM.GZ FILES TO RETAIN ONLY PRIMARY READS WITH READ LENGTH = READ_LENGTH_FILTER  ###                                                     
####################################################################################################################################################################################################################

def decompress_cx_sam_gz_file(gz_file_path, output_path):
    """Decompress a *.sam.gz file"""
    with gzip.open(gz_file_path, 'rb') as f_in, open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def filter_cx_sam_file(sam_file_path, output_path):
    """Filters SAM file based on read length and mapping quality."""
    try:
        with pysam.AlignmentFile(sam_file_path, "r", check_sq=False) as infile, open(output_path, "w") as outfile:
            ### write headers ### 
            for header_line in infile.text.splitlines():
                outfile.write(header_line + "\n")

            ### apply filtering ###
            for read in infile.fetch(until_eof=True):
                if (
                    not read.is_secondary and
                    not read.is_supplementary and
                    not read.is_unmapped and
                    (len(read.query_sequence) == read_length_filter if read_length_filter else True) and
                    read.mapping_quality >= mapq_threshold
                ):
                    outfile.write(read.to_string() + "\n")
    except Exception as e:
        print(f"Error processing {sam_file_path}: {e}")

def compress_cx_sam_to_gz(sam_file_path, gz_output_path):
    """Compress a filtered SAM file into *.sam.gz format"""
    with open(sam_file_path, 'rb') as f_in, gzip.open(gz_output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def process_single_cx_sam(file_name):
    """Processes a single SAM or gzipped SAM file"""
    file_path = os.path.join(cx_sam_dir, file_name)
    
    ### step 1: unzip input SAM files if they aren't already ###
    if file_name.endswith('.gz'):
        sam_temp_path = os.path.join(unzipped_sam_dir, file_name.replace('.gz', ''))
        decompress_cx_sam_gz_file(file_path, sam_temp_path)
        input_sam_path = sam_temp_path
    else:
        input_sam_path = file_path

    filtered_sam_path = os.path.join(filtered_sam_unzipped_dir, os.path.basename(input_sam_path))
    filtered_gz_output_path = os.path.join(filtered_sam_gzipped_dir, os.path.basename(input_sam_path) + ".gz")

    ### step 2: filter the unzipped SAM files ###
    filter_cx_sam_file(input_sam_path, filtered_sam_path)

    ### step 3: gzip filtered SAM files ###
    compress_cx_sam_to_gz(filtered_sam_path, filtered_gz_output_path)

def process_cx_sam_files():
    """Processes SAM files in parallel using ProcessPoolExecutor"""
    files = [f for f in os.listdir(cx_sam_dir) if f.endswith('.gz') or f.endswith('.sam')]

    max_workers = min(len(files), os.cpu_count())                                                                          

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_cx_sam, file_name): file_name for file_name in files}

        for future in futures:
            try:
                future.result() 
            except Exception as e:
                print(f"Error processing {futures[future]}: {e}")


####################################################################################################################################################################################################################
                                                                                ### GET READ COUNTS TO CALCULATE VARIANT FREQUENCIES ###                                                                            
####################################################################################################################################################################################################################

def get_read_counts(file_path):
    """Process a single SAM file to calculate the average depth"""
    data = []
    
    ### open the filtered *.sam.gz files and process each read ###
    with pysam.AlignmentFile(file_path, "r") as samfile:
        for read in samfile:
            start_position = read.reference_start                                                                                   ### assign positions for each base in read ###
            if read.query_qualities is not None:
                for i, quality in enumerate(read.query_qualities):
                    position_value = start_position + i
                    score = 1 if quality >= baseq_threshold else 0                                                                  ### award score of 1 for position if base in that position has quality > baseq_threshold; else score of 0 awarded ###
                    data.append((read.query_name, position_value, score))                                                           ### append position and corresponding score of each read in data list ###
    
    ### create a df with columns "read_id", "position_value" and "score" ###
    df = pd.DataFrame(data, columns=["read_id", "position_value", "score"])
    
    ### remove duplicate rows based on read_id, position_value, and score ###
    df = df.drop_duplicates(subset=["read_id", "position_value", "score"])
    
    ### obtain total score for each position across all reads in SAM file ###
    grouped = df.groupby("position_value")["score"].sum()
    
    ### remove any position values that are smaller than start_bp and larger than end_bp i.e. outside of mutagenesis region ###
    grouped = grouped[(grouped.index >= start_bp) & (grouped.index <= end_bp)]
    
    ### calculate average depth per position by summing all scores in grouped and dividing by length of mutagenesis region ###
    total_score = grouped.sum()
    num_positions = len(grouped)
    average_score = total_score / (end_bp - start_bp + 1) if num_positions > 0 else 0
    
    ### round the average score to the nearest integer ###
    average_score = round(average_score)
    
    ### return the filename (without extension) and the average depth ###
    return os.path.basename(file_path).replace('.sam.gz', ''), average_score

def get_read_counts_parallel(filtered_sam_gzipped_dir):
    """Process SAM files in parallel using ProcessPoolExecutor to calculate average depth and returns a dictionary with filenames as keys and corresponding average depths as values"""
    sam_files = [os.path.join(filtered_sam_gzipped_dir, f) for f in os.listdir(filtered_sam_gzipped_dir) if f.endswith('.sam.gz')]
    read_counts = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(get_read_counts, file_path): file_path
            for file_path in sam_files
        }

        for future in as_completed(futures):
            file_path = futures[future]
            try:
                filename, average_score = future.result()
                read_counts[filename] = average_score
            except Exception as e:
                print(f"Failed to process {file_path}: {e}")

    return read_counts


######################################################################################################################################################################################################################
                                                                                        ### MAP VARIANTS TO CDNA POSITION ###                                                                                                                                                 
######################################################################################################################################################################################################################

### create dictionary to map bp to cDNA positions ###
def map_pos_to_cdna(cDNA_file):
    cdna_file_df = pd.read_csv(cDNA_file)
    cdna_file_df = cdna_file_df[["bp", "cDNA"]]
    cdna_file_df['cDNA'] = cdna_file_df['cDNA'].astype(str)
    cdna_file_df = cdna_file_df.drop_duplicates(subset=['bp']).reset_index(drop=True)
    bp_to_cdna = dict(zip(cdna_file_df['bp'], cdna_file_df['cDNA']))
    return bp_to_cdna
bp_to_cdna = map_pos_to_cdna(cDNA_file)

def convert_to_cdna(variants):
    """Converts each genomic variant separated by ', ' in cs_long_sub_baseq30 to cDNA positional variants using bp_to_cdna dictionary e.g. 168065*ag -> 2003*ag"""
    cdna_variants = []
    for variant in variants.split(', '):
        ### split genomic position and mutation part (e.g., "168065*ag" -> "168065" and "ag") ###
        try:
            pos, mutation = variant.split('*')
            pos = int(pos)                                                                                                          ### converts genomic position to an integer to look up corresponding cDNA position using bp_to_cdna ###
            cdna_pos = bp_to_cdna.get(pos, pos)                                                                                     ### use mapped cDNA position if available; otherwise keep original to prevent crashing script ###
            cdna_variants.append(f"{cdna_pos}*{mutation}")
        except ValueError:
            ### handle any unexpected format errors ###
            cdna_variants.append(variant)                                                                                           ### keep variant as is if it cannot be processed ###
    return ', '.join(cdna_variants)

def reformat_variant(variants):
    """Reformats each cDNA variant e.g. 2003*ag -> 2003A>G""" 
    formatted_variants = []
    for variant in variants.split(', '):
        try:
            cdna_position, bases = variant.split('*')
            ref_base, obs_base = bases[0].upper(), bases[1].upper()  # Capitalize the bases
            formatted_variants.append(f"{cdna_position}{ref_base}>{obs_base}")
        except ValueError:
            ### skip variant if there's a format error ###
            continue
    ### join reformatted variants back into a string ###
    return ', '.join(formatted_variants)


######################################################################################################################################################################################################################
                                                                                          ### DIRTY MATH VARIANT CALLING ###                                                                                         
######################################################################################################################################################################################################################

def open_sam_file(sam_file_path):
    """Open a SAM file, supporting both .sam and .sam.gz formats."""
    if sam_file_path.endswith('.gz'):
        return gzip.open(sam_file_path, 'rt')
    else:
        return open(sam_file_path, 'r')

def calculate_cs_short_sub(row):
    """Parses cs_short tag of each read to extract substitutions (denoted by '*') and their positions (assigned based on read's starting position)"""
    cs_tag = row['cs_short']
    current_position = row['read_start']
    substitutions = []

    ### regular expression to match each component in short cs_tag ###
    pattern = re.compile(r'(\d+)|(\*[a-z]{2})|(\+[a-z]+)|(-[a-z]+)')
    
    for match in pattern.finditer(cs_tag):
        if match.group(1):                                                                                                          ### ID exact matches (e.g., ":50") ###
            current_position += int(match.group(1))
        elif match.group(2):                                                                                                        ### ID substitutions (2 letters after '*') ###
            sub_base = match.group(2)[1:]                                                                                           ### extract base change substitutions after '*' ###
            substitutions.append(f"{current_position}*{sub_base}")
            current_position += 1                                                                                                   ### increment current position by 1 to move on to next position ###
        elif match.group(3):                                                                                                        ### ID insertions (denoted by '+') ###
            pass                                                                                                                    ### ignore inserted bases / pass ###
        elif match.group(4):                                                                                                        ### ID deletions (denoted by '-') ###
            del_length = len(match.group(4)) - 1                                                                                    ### compute number of bases deleted ###
            current_position += del_length                                                                                          ### increment reference position by adding back number of bases deleted ###

    ### join all substitutions for the row, separated by commas ###
    return ', '.join(substitutions)

def calculate_cs_long_sub(row):
    """Parses the `cs_masked` tag of each read to extract substitutions (denoted by '*') of qualities > baseq_threshold, together with their positions"""
    cs_masked = row['cs_masked']
    current_position = row['read_start']
    substitutions = []

    ### regular expression to match each component in cs_masked ###
    pattern = re.compile(r'(=[A-Za-z]+)|(\*[a-z]{2})|(\+[a-z]+)|(-[a-z]+)')
    
    for match in pattern.finditer(cs_masked):
        if match.group(1):                                                                                                          ### ID exact matches (e.g., "=CGATCG")
            ### extract matched segment and increment position for each base except '=' ###
            match_segment = match.group(1)[1:]                                                                                      ### ignore '=' ###
            current_position += len(match_segment)
        elif match.group(2):                                                                                                        ### ID substitutions (2 letters after '*') ###
            sub_base = match.group(2)[1:]                                                                                           ### ID base change substitutions after '*' ###
            substitutions.append(f"{current_position}*{sub_base}")
            current_position += 1                                                                                                   ### increment current position by 1 to move on to next position ###
        elif match.group(3):                                                                                                        ### ID insertions (denoted by '+') ###
            pass                                                                                                                    ### ignore inserted bases / pass ###
        elif match.group(4):                                                                                                        ### ID deletions (denoted by '-') ###
            del_length = len(match.group(4)) - 1                                                                                    ### compute number of bases deleted ###
            current_position += del_length                                                                                          ### increment reference position by adding back number of bases deleted ###

    ### join all substitutions for the row, separated by commas ###
    return ', '.join(substitutions)

# Function to remove variants with baseq < 30
def remove_n_variants(variants):
    """Remove variants with baseq < 30 as denoted by variant string containing 'n'"""
    filtered_variants = [variant for variant in variants.split(', ') if 'n' not in variant]                                         ### split variant list by ', ' and remove individual variants containing "n" ###
    ### re-join filtered variants back into a list separated by ', ' ###
    return ', '.join(filtered_variants)

def filter_variants(variants):
    """Parses through position value (pos) of each variant to remove those variants with pos outside of mutagenisis region"""
    filtered_variants = []
    for variant in variants.split(', '):
        ### split the position and mutation part ###
        try:
            pos, mutation = variant.split('*')
            pos = int(pos)
            ### keep variant only if position is within the specified range ###
            if start_bp <= pos <= end_bp:
                filtered_variants.append(f"{pos}*{mutation}")
        except ValueError:
            ### skip variant if there are format errors ###
            continue
    ### re-join the filtered variants back into a list separated by ', ' ###
    return ', '.join(filtered_variants)

def sam_to_dataframe(sam_file_path):
    """Convert a SAM file into a DataFrame with columns: 'read_id', 'flag', 'read_start', 'seq', 'base_qualities', 'md', 'cigar'"""
    data = []

    with pysam.AlignmentFile(sam_file_path, "rb") as samfile:
        for read in samfile:
            ### extract the relevant fields ###
            read_id = read.query_name
            flag = read.flag
            read_start = read.reference_start + 1  
            seq = read.query_sequence
            base_qualities = read.qual
            cigar = read.cigarstring
            length = read.query_alignment_length

            ### extract 'MD' tag if present ###
            md_tag = None
            for tag, value in read.get_tags():
                if tag == "MD":
                    md_tag = value
                    break
            
            ### append information to data list ###
            data.append([read_id, flag, read_start, seq, base_qualities, md_tag, cigar, length])
    
    ### create a df from the list of data ###
    df = pd.DataFrame(data, columns=['read_id', 'flag', 'read_start', 'seq', 'base_qualities', 'md', 'cigar', 'length'])

    ### sort by read_id ###
    df.sort_values(by='read_id', inplace=True)

    ### apply cstag.call to obtain 'cs_short' tag for each read ###
    df['cs_short'] = df.apply(
        lambda row: cstag.call(row['cigar'], row['md'], row['seq'])
        if row['cigar'] and row['md'] and row['seq'] else None,
        axis=1
    )

    ### apply cstag.call to obtain 'cs_long' tag for each read ###
    df['cs_long'] = df.apply(
        lambda row: cstag.call(row['cigar'], row['md'], row['seq'], long=True)
        if row['cigar'] and row['md'] and row['seq'] else None,
        axis=1
    )

    ### remove rows where 'cs_tag' contains no '*' i.e. reads with no subsitutions ###
    initial_row_count = len(df)
    df = df[df['cs_long'].str.contains(r'\*')]
    removed_row_count = initial_row_count - len(df)                                                                                 ### stores number of reads without substitutions as removed_row_count ###
    df.reset_index(drop=True, inplace=True)

    ### apply cstag.mask to mask bases in cs_long tags with base quality < baseq_threshold ###`
    df['cs_masked'] = df.apply(
        lambda row: cstag.mask(row['cs_long'], row['cigar'], row['base_qualities'], baseq_threshold)
        if row['cs_long'] and row['cigar'] and row['base_qualities'] else None,
        axis=1
    )

    df = df[["read_id", "flag", "read_start","cs_short", "cs_masked"]]
    df['cs_short_sub'] = df.apply(calculate_cs_short_sub, axis=1)                                                                   ### compute variants from cs_short tags ###
    df['cs_long_sub'] = df.apply(calculate_cs_long_sub, axis=1)                                                                     ### compute variants from cs_long tags ###
    df['cs_long_sub_baseq30'] = df['cs_long_sub'].apply(remove_n_variants)                                                          ### remove variants with baseq < 30 ###                   
    df['cs_long_sub_baseq30'] = df['cs_long_sub_baseq30'].apply(filter_variants)                                                    ### remove variants outside of position range ###


    ### create a copy of the DataFrame to process 'cs_long_sub' ###
    processed_df = df[["read_id", "flag", "cs_long_sub_baseq30"]].copy()

    ### remove rows where 'cs_long_sub' is NaN or an empty string after stripping whitespace ###
    processed_df['cs_long_sub_baseq30'] = processed_df['cs_long_sub_baseq30'].fillna('').str.strip() 
    processed_df = processed_df[processed_df['cs_long_sub_baseq30'] != ""]    

    processed_df['cdna_variants'] = processed_df['cs_long_sub_baseq30'].apply(convert_to_cdna)                                      ### convert to cDNA variants ###
    processed_df['cdna_variants'] = processed_df['cdna_variants'].apply(reformat_variant)                                           ### reformat cdna variants ###  


    ################ FIND COMMON VARIANTS SHARING THE SAME READ_ID ################
    
    ### drop 'flag' and 'cs_long_sub_baseq30' columns ###
    processed_df = processed_df[['read_id', 'cdna_variants']]

    def find_common_variants(variants_list):
        """ID common variants shared between paired reads; assumes readID for paired reads are the same"""
        ### split each string of variants by ", " and convert to a set to find the intersection ###
        sets = [set(variants.split(', ')) for variants in variants_list]
        common_variants = set.intersection(*sets) if sets else set()
        ### re-join the common variants back into a single comma-separated string ###
        return ', '.join(sorted(common_variants))

    processed_df_grouped = processed_df.groupby('read_id').agg({'cdna_variants': find_common_variants}).reset_index()

    raw_output_file = os.path.basename(sam_file_path).replace('.sam.gz', '.csv')
    raw_output_path = os.path.join(raw_output_dir, raw_output_file)
    processed_df_grouped.to_csv(raw_output_path, index=False)

    return processed_df_grouped, removed_row_count

def process_single_sam_file(sam_file):
    """Processes a single SAM file and returns the filename, processed DataFrame, and removed row count"""
    filename = os.path.basename(sam_file).split('.')[0]
    processed_df_grouped, removed_row_count = sam_to_dataframe(sam_file)
    return filename, processed_df_grouped, removed_row_count

def process_sam_cs_variants(filtered_sam_gzipped_dir):
    """Process each .sam.gz file in a directory, convert to a grouped DataFrame, and collect removed row counts"""
    sam_files = [os.path.join(filtered_sam_gzipped_dir, f) for f in os.listdir(filtered_sam_gzipped_dir) if f.endswith('.sam.gz')]
    dfs = []
    removed_counts = {}
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_sam_file, sam_file): sam_file for sam_file in sam_files}

        for future in as_completed(futures):
            sam_file = futures[future]
            try:
                filename, processed_df_grouped, removed_row_count = future.result()
                dfs.append(processed_df_grouped)
                removed_counts[filename] = removed_row_count
            except Exception as exc:
                print(f"{sam_file} generated an exception: {exc}")

    return dfs, removed_counts

######################################################################################################################################################################################################################
                                                                              ### GET COUNTS OF EACH UNIQUE CDNA VARIANT COMBINATION ###                                                                              
######################################################################################################################################################################################################################

def get_variants(file_path, variants_per_rep_dir, nm_histograms_dir, removed_counts, read_counts_dict):
    """
    Count number of times each unique 'cdna_variants' value was found in each csv file of raw_output_dir

    Parameters:
    - file_path: File path to each CSV file in raw_output_dir
    - variants_per_rep_dir: Path to the directory where processed files will be saved.
    - removed_counts: Dictionary containing removed row counts with sample IDs as keys.
    - nm_histograms_dir: Path to the directory where NM histograms will be saved as PDF files.
    """
    filename = os.path.basename(file_path)
    sample_ID = filename.split('.csv')[0]                                                                                           ### get sample_ex_rep e.g. iPSC_Ex9_Rep1, A_Ex11_Rep5 ###

    ### read each CSV file from raw_output_dir into a df and process ###
    df = pd.read_csv(file_path, usecols=["cdna_variants"])                                                                          ### csv --> df and keep only 'cdna_variants' ###
    processed_df = df.groupby("cdna_variants").size().reset_index(name='count')                                                     ### group by unique 'cdna_variants' and count occurrences of each unique cdna_variants value ###
    processed_df = processed_df.sort_values(by='count', ascending=False)

    ### add sample_ID to 'count' column name and use read_counts_dict to map corresponding depth ###
    processed_df.rename(columns={'count': f'count_{sample_ID}'}, inplace=True)                                          
    processed_df[f'depth_{sample_ID}'] = read_counts_dict[sample_ID] if sample_ID in read_counts_dict else print(f"{sample_ID} not found in read_counts_dict")

    output_file_path = os.path.join(variants_per_rep_dir, filename)
    processed_df.to_csv(output_file_path, index=False)

    ### calculate number of mutations (NM) for each cdna_variant ###
    processed_df['NM'] = processed_df['cdna_variants'].apply(lambda x: 1 + x.count(','))
    processed_df = processed_df.drop(columns=['cdna_variants', f'depth_{sample_ID}'])

    ### get counts of unique NMs per sample_ID ###
    grouped_df = processed_df.groupby('NM', as_index=False).agg({f'count_{sample_ID}': 'sum'})

    ### use removed_counts to add a row to group_df where NM=0 ###
    nm_zero_count = removed_counts.get(sample_ID, 0)
    nm_zero_row = pd.DataFrame({'NM': [0], f'count_{sample_ID}': [nm_zero_count]})
    grouped_df = pd.concat([grouped_df, nm_zero_row], ignore_index=True)

    ### plot histogram to visualise NM distribution per sample_rep ###
    plt.figure(figsize=(10, 6))
    plt.bar(grouped_df['NM'], grouped_df[f'count_{sample_ID}'], width=0.8, align='center', color='#DDD5F3')

    ### set y limits ###
    max_y_val = grouped_df[f'count_{sample_ID}'].max()
    y_buffer = int(max_y_val * 0.05)                                                                                                    ### set buffer to 5% of max value ###
    ylim = max_y_val + y_buffer
    plt.ylim(0, ylim)
    plt.xlabel("Number of Mutations per Variant", labelpad=10, fontsize=10)
    plt.ylabel("Count", labelpad=10, fontsize=10)
    plt.title(f"Distribution of Number of Mutations per Variant for {sample_ID}", fontsize=13, pad=13)
    plt.xticks(grouped_df['NM'], fontsize=8)

    plt.yticks(range(0, grouped_df[f'count_{sample_ID}'].max() + 999, 1000), fontsize=8)

    ### plot grid behind bars ###
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.gca().set_axisbelow(True)

    ### save histograms as pdf files ###
    histogram_path = os.path.join(nm_histograms_dir, f"{sample_ID}_NM_histogram.pdf")
    plt.savefig(histogram_path, format='pdf', bbox_inches='tight')
    plt.close()

def get_variants_parallel(raw_output_dir, variants_per_rep_dir, nm_histograms_dir, removed_counts, read_counts_dict):
    """
    Process CSV files in parallel using ProcessPoolExecutor.
    """
    csv_files = [os.path.join(raw_output_dir, f) for f in os.listdir(raw_output_dir) if f.endswith('.csv')]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(get_variants, file_path, variants_per_rep_dir, nm_histograms_dir, removed_counts, read_counts_dict): file_path for file_path in csv_files}

        for future in as_completed(futures):
            try:
                filename = future.result()
                # print(f"Processed {filename}")
            except Exception as e:
                print(f"Exception processing {futures[future]}: {e}")

def combine_variant_across_reps(group_name, files, variants_per_sample_dir):
    """ 
    Calculate total count, depth and frequency of each variant within a single sample and exon by merging the data of replicates.
    """
    merged_df = pd.DataFrame()

    for file in files:
        df = pd.read_csv(file)
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='cdna_variants', how='outer').fillna(0)

    ### calculate total count of each variant ###
    count_columns = [col for col in merged_df.columns if col.startswith('count_')]
    merged_df[f'total_count_{group_name}'] = merged_df[count_columns].sum(axis=1)

    ### calculate total depth of each variant ###
    depth_columns = [col for col in merged_df.columns if col.startswith('depth_')]
    non_zero_rows = merged_df[depth_columns].all(axis=1)

    ### calculate total_depth for rows where all values are non-zero ###
    if non_zero_rows.any():
        total_depth_value = merged_df.loc[non_zero_rows, depth_columns].sum(axis=1).iloc[0]
        merged_df[f'total_depth_{group_name}'] = total_depth_value
    else:                                                                                                                           ### if no rows meet condition, default to 0 for total depth ###
        merged_df[f'total_depth_{group_name}'] = 0  

    ### remove depth_* columns ###
    merged_df = merged_df.drop(columns=depth_columns)

    ### calculate frequency for each variant across reps (total count/ total depth) ###
    merged_df[f'total_frequency_{group_name}'] = merged_df[f'total_count_{group_name}'] / merged_df[f'total_depth_{group_name}']
    merged_df = merged_df[['cdna_variants'] + sorted(count_columns) + [f'total_count_{group_name}'] + [f'total_depth_{group_name}'] + [f'total_frequency_{group_name}']]

    ### save to csv file ###
    output_file = os.path.join(variants_per_sample_dir, f"{group_name}.csv")
    merged_df.to_csv(output_file, index=False)

    ### remove rows where total_frequency < 10^(-4) ###
    merged_MAF_df = merged_df[merged_df[f'total_frequency_{group_name}'] >= seq_error]
    output_MAF_file = os.path.join(variants_per_sample_dir, f"{group_name}_MAF.csv")
    merged_MAF_df.to_csv(output_MAF_file, index=False)

    return group_name

def combine_variant_across_reps_parallel(variants_per_rep_dir, variants_per_sample_dir):
    """
    Combine variant counts across replicates for each sample and exon, with parallel processing.
    """
    ### collect all csv files of replicates to be merged ###
    csv_files = glob.glob(os.path.join(variants_per_rep_dir, '*.csv'))

    ### group files by sample and exon ###
    grouped_files = {}
    for file in csv_files:
        filename = os.path.basename(file)
        sample_exon = filename.split('_Rep')[0]
        grouped_files.setdefault(sample_exon, []).append(file)

    ### use ProcessPoolExecutor to process groups in parallel ###
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(combine_variant_across_reps, group_name, files, variants_per_sample_dir): group_name
            for group_name, files in grouped_files.items()
        }

        for future in as_completed(futures):
            group_name = futures[future]
            try:
                future.result()
                # print(f"Processed group: {group_name}")
            except Exception as exc:
                print(f"Exception processing group {group_name}: {exc}")



######################################################################################################################################################################################################################
                                                                                      ### FORMAT VARIANTS INTO HGVS NOTATION ###                                                                                     
######################################################################################################################################################################################################################

def format_refdf(cDNA_file):
    refdf = pd.read_csv(cDNA_file)
    refdf['cDNA'] = refdf['cDNA'].astype(str)
                         
    ### extract numerical values from 'cDNA' ###
    def extract_numerical_value(cdna):
        match = re.search(r'\d+[-]?\d*', cdna)
        return match.group(0) if match else None
    
    refdf['cDNA'] = refdf['cDNA'].apply(extract_numerical_value)
    refdf = refdf[["bp", "cDNA", "reference", "observed_nucleotide", "AA_Consequence"]]

    return refdf

def HGVS_prep(variants_per_sample_dir, cDNA_file, HGVS_dir):
    """
    Reformat variants to prep for transformation into HGVS notation by determining whether the variants fall in the same codon number

    Parameters:
    variants_per_sample_dir (str): 
    cDNA_file (str): The path to the CSV file containing 'bp', 'cDNA', 'reference' and 'AA_Consequence' data
    HGVS_dir: output dir of files containing HGVS formatted variants per sample

    Workflow:
    * Reads cdna_file into df, retaining only the "cDNA", "reference" and "AA_Consequence" columns relevant for HGVS notation 
    * Extracts the AA positions from "AA_Consequence" and sets "cDNA" as df index for efficient lookups
    * Calls process_changes_generalized function to process each row of variants, grouping them by "AA_position" and formatting variants as SNV (different codon) or MNV (same codon)
    * Calls prepend_hgvs_notation function to format the variants into HGVS notation for VEP processing
    """
    ### Read cdna_file into a DataFrame ###
    refdf = format_refdf(cDNA_file)
    refdf = refdf.loc[:, refdf.columns.intersection(['cDNA','reference', 'AA_Consequence'])]
    refdf['AA_Position'] = refdf['AA_Consequence'].str.extract(r'(\d+-?\d*)')
    refdf.drop(columns=['AA_Consequence'], inplace=True)
    refdf.drop_duplicates(inplace=True)
    refdf.reset_index(drop=True, inplace=True)
    refdf = refdf.set_index("cDNA")
    refdf["AA_Position"] = refdf["AA_Position"].astype(str)

    ### Inner function 1 ###
    def process_changes_generalized(row, refdf):
        """
        Splits base changes of each variant to identify the variant positions and base changes separately, groups and sorts positions, 
        determine if they are SNV or MNV based on the mapped AA_position, before returning final HGVS notation string

        Parameters:
        Row from each row of sample_exon.csv from variants_per_sample_dir and reference df ("refdf")
        """
        changes = row.split(", ")
        positions = [change.split(">")[0][:-1] for change in changes]
        mutations = [change.split(">")[1] for change in changes]

        aa_positions = refdf.loc[positions, "AA_Position"].tolist()
        grouped_positions = {}

        for pos, aa_pos, mut in zip(positions, aa_positions, mutations):
            if aa_pos not in grouped_positions:
                grouped_positions[aa_pos] = []
            grouped_positions[aa_pos].append((pos, mut))

        new_changes = []
        for aa_pos, pos_mut_pairs in grouped_positions.items():
            pos_mut_pairs.sort(key=lambda x: int(x[0].split('-')[0]))

            ### construct variant strings ###
            if len(pos_mut_pairs) == 1:
                pos, mut = pos_mut_pairs[0]
                change = f"{pos}{refdf.loc[pos, 'reference']}>{mut}"
            else:
                ### check for consecutive positions and construct the string accordingly ###
                consecutive = True
                refs = ""
                for i in range(len(pos_mut_pairs) - 1):
                    current_pos, current_mut = pos_mut_pairs[i]
                    next_pos, next_mut = pos_mut_pairs[i + 1]
                    refs += current_mut
                    if int(next_pos.split('-')[0]) - int(current_pos.split('-')[0]) > 1:
                        consecutive = False
                        missing_refs = ''.join(refdf.loc[str(pos), 'reference'] for pos in range(int(current_pos.split('-')[0]) + 1, int(next_pos.split('-')[0])))
                        refs += missing_refs
                ### add mutation of last position ###
                refs += next_mut

                if consecutive:
                    change = f"{pos_mut_pairs[0][0]}_{pos_mut_pairs[-1][0]}delins{refs}"
                else:
                    change = f"{pos_mut_pairs[0][0]}_{pos_mut_pairs[-1][0]}delins{refs}"

            new_changes.append(change)
        return ", ".join(new_changes)

    ### Inner function 2 ###
    def prepend_hgvs_notation(int_var_str):
        """
        Prepends transcript to mutation to obtain HGVS notation for each mutation
        """
        mutations = int_var_str.split(", ")                                                                                         ### split the string by comma to process each mutation individually ###
        hgvs_mutations = [transcript + mutation for mutation in mutations]                                                          ### prepend the HGVS notation to each mutation ###
        return ", ".join(hgvs_mutations)   

    ### Inner function 3 ###
    def format_int_var_hgvs_to_txt(sample_var_df, HGVS_txt_path):
        """
        Writes 'Int_var_HGVS' from each row of sample_var_df into a new line in a "*.txt" file, where each set of mutations (variants) is separated by empty line

        Parameters:
        Takes sample_var_df as input file and outputs txt file to HGVS_txt_path
        """
        with open(HGVS_txt_path, 'w') as f:
            for _, row in sample_var_df.iterrows():
                int_var_hgvs = row['Int_var_HGVS']
                hgvs_mutations = int_var_hgvs.split(", ")
                for mutation in hgvs_mutations:
                    f.write(mutation + '\n')
                f.write('\n')

    for file_path in glob.glob(os.path.join(variants_per_sample_dir, "*_MAF.csv")):                                                 ### loop through each *_MAF.csv file in variants_per_sample_dir ###
        sample_var_df = pd.read_csv(file_path)
        sample_ID = os.path.basename(file_path).replace("_MAF.csv", "")                                                             ### get the sample ID from the file name ###

        sample_var_df = sample_var_df[["cdna_variants", f'total_count_{sample_ID}', f'total_depth_{sample_ID}', f'total_frequency_{sample_ID}']]

        sample_var_df['Int_var'] = sample_var_df['cdna_variants'].apply(lambda x: process_changes_generalized(x, refdf))            ### processes nm_df with 'process_changes_generalized' and 'prepend_hgvs_notation' functions column wise ###
        sample_var_df['Int_var_HGVS'] = sample_var_df['Int_var'].apply(prepend_hgvs_notation)

        HGVS_csv_path = os.path.join(HGVS_dir, f"{sample_ID}_HGVS.csv")
        sample_var_df.to_csv(HGVS_csv_path, index=False, encoding='utf-8-sig')

        HGVS_txt_path = os.path.join(HGVS_dir, f"{sample_ID}_HGVS.txt")
        format_int_var_hgvs_to_txt(sample_var_df, HGVS_txt_path)                                                                    ### calls 'format_int_var_hgvs_to_txt' function to write the 'Int_var_HGVS' column from nm_df into a "_HGVS_for_ensembl.txt" file in output *_HGVS_for_ensembl directory for VEP ###



#####################################################################################################################################################################################################################
                                                                                               ### RUN ENSEMBL VEP ###                                                                                               
#####################################################################################################################################################################################################################

### VEP options ###
vep_options = "--database --everything --no_stats --quiet"

### Define function to run VEP ###
def run_vep(HGVS_txt_file, vep_txt_file, vep_path, vep_options):
    """
    Execute VEP on input *_HGVS.txt file using 'singularity' command.

    Parameters:
    HGVS_txt_file (str): Path of input *_HGVS.txt file
    vep_txt_file (str): Path of VEP output file
    vep_path (str): Path to vep.sif
    vep_options (str): String containing any additional options to pass to the VEP command
    """
    command = f"singularity exec {vep_path} vep {vep_options} --input_file {HGVS_txt_file} --output_file {vep_txt_file}"
    try:
        subprocess.run(command, shell=True, check=True)
        # print(f"VEP output saved to: {vep_txt_file}")
    except Exception as e:
        print(f"Error running VEP for {HGVS_txt_file}: {e}")

### Parallel processing for multiple files ###
def run_vep_in_parallel(HGVS_dir, vep_txt_dir, vep_path, vep_options):
    """
    Run VEP in parallel for all *_HGVS.txt files in the input directory.

    Parameters:
    HGVS_dir (str): Directory containing *_HGVS.txt files.
    vep_txt_dir (str): Directory to store VEP output txt files.
    vep_path (str): Path to vep.sif
    vep_options (str): VEP options to be passed to the command
    """
    ### collect all *_HGVS.txt input files ###
    hgvs_input_files = [os.path.join(HGVS_dir, f) for f in os.listdir(HGVS_dir) if f.endswith('_HGVS.txt')]

    ### define output file paths ###
    tasks = []
    for hgvs_file in hgvs_input_files:
        base_name = os.path.basename(hgvs_file).replace('_HGVS.txt', '_VEP.txt')
        vep_txt_file = os.path.join(vep_txt_dir, base_name)
        tasks.append((hgvs_file, vep_txt_file))

    ### use ProcessPoolExecutor for parallel processing ###
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(run_vep, task[0], task[1], vep_path, vep_options): task
            for task in tasks
        }

        for future in as_completed(future_to_task):
            hgvs_file, vep_txt_file = future_to_task[future]
            try:
                future.result()  
            except Exception as e:
                print(f"Failed to process {hgvs_file}: {e}")


######################################################################################################################################################################################################################
                                                                                              ### PROCESS VEP OUTPUT ###                                                                                              
######################################################################################################################################################################################################################

### Define function to process each row of VEP output txt file ###
def process_vep_data_row(row):
    """
    Splits additional data from VEP output into constituent columns and check for necessary fields 

    Parameters:
    row (list): A row from VEP result

    Returns:
    tuple: The row without additional data, additional data as a dict, and a flag indicating if row should be skipped based on the absence of "MANE_SELECT" column
    """
    extra_data = row[-1]
    extra_columns = dict(item.split('=', 1) for item in extra_data.split(';') if '=' in item)
    ### determine if the "MANE_SELECT" key is missing in the dictionary ###
    skip_row = ('MANE_SELECT' not in extra_columns)
    return row[:-1], extra_columns, skip_row

### Define function to convert VEP output txt file into CSV file for easier visualisation and downstream data analysis ###
def convert_vep_txt_to_csv(HGVS_txt_file, vep_txt_file, vep_csv_file, header_lines=73):
    """
    Converts *_VEP.txt file to *_VEP.csv file

    Parameters:
    HGVS_txt_file (str): Path of input *_HGVS.txt file
    vep_output_txt (str): Path of input *_VEP.txt file 
    vep_csv_file (str): Output file path for *_VEP.csv files
    header_lines (int, optional): Number of header lines in the VEP file; defaults to 73

    Workflow:
    * Reads empty line positions from the input *_HGVS.txt file
    * Iterates over the  *_VEP_output.txt file, processes each row, and keeps only necessary data
    * Creates a CSV writer and writes headers
    * Writes processed data and additional columns to output *_VEP_output.csv file
    
    If any exception occurs during this process, it is caught and printed.
    """
    try:
        empty_line_positions = []
        with open(HGVS_txt_file, 'r') as hgvs_file:
            for line in hgvs_file:
                empty_line_positions.append(not line.strip())

        all_keys = set()
        processed_data = []
        with open(vep_txt_file, 'r') as vep_infile:
            vep_reader = csv.reader(vep_infile, delimiter='\t')
            for i, row in enumerate(vep_reader):
                if i <= header_lines:
                    continue
                processed_row, extra_columns, skip_row = process_vep_data_row(row)
                if skip_row:
                    continue
                for key in extra_columns.keys():
                    all_keys.add(key)
                processed_data.append((processed_row, extra_columns))

        with open(vep_csv_file, 'w', newline='') as vep_outfile:
            vep_writer = csv.writer(vep_outfile)
            headers = ["Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position",
                       "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation"] + sorted(all_keys)
            vep_writer.writerow(headers)
            filtered_row_index = 0
            for is_empty_line in empty_line_positions:
                    if not is_empty_line and filtered_row_index < len(processed_data):
                        row_data, extra_data = processed_data[filtered_row_index]
                        vep_writer.writerow(row_data + [extra_data.get(key, '') for key in sorted(all_keys)])
                        filtered_row_index += 1
                    else:
                        vep_writer.writerow([])
    except Exception as e:
        print(f"Error in convert_txt_to_csv: {e}")

def process_single_vep_file(vep_txt_file, HGVS_dir, vep_csv_dir):
    """
    Process a single *_VEP.txt file and generate a *_VEP.csv file.

    Parameters:
    vep_txt_file (str): Path to the VEP output txt file.
    HGVS_dir (str): Directory containing the *_HGVS.txt files.
    vep_csv_dir (str): Directory to store the converted *_VEP.csv files.
    """
    base_name = os.path.basename(vep_txt_file).replace("_VEP.txt", "")
    HGVS_txt_file = os.path.join(HGVS_dir, f"{base_name}_HGVS.txt")
    vep_csv_file = os.path.join(vep_csv_dir, f"{base_name}_VEP.csv")

    if not os.path.exists(HGVS_txt_file):
        print(f"Missing corresponding HGVS file for {vep_txt_file}. Skipping.")
        return

    try:
        convert_vep_txt_to_csv(HGVS_txt_file, vep_txt_file, vep_csv_file)
        # print(f"Converted {vep_txt_file} to {vep_csv_file}.")
    except Exception as e:
        print(f"Error processing {vep_txt_file}: {e}")

### Define a function to process all *_VEP.txt files in parallel ###
def process_vep_files_in_parallel(vep_txt_dir, HGVS_dir, vep_csv_dir):
    """
    Process all *_VEP.txt files in the specified directory in parallel.

    Parameters:
    vep_txt_dir (str): Directory containing *_VEP.txt files.
    HGVS_dir (str): Directory containing *_HGVS.txt files.
    vep_csv_dir (str): Directory to store the converted *_VEP.csv files.
    max_workers (int): Number of parallel processes to use.
    """
    ### collect all *_VEP.txt files ###
    vep_txt_files = [os.path.join(vep_txt_dir, f) for f in os.listdir(vep_txt_dir) if f.endswith("_VEP.txt")]

    ### use ProcessPoolExecutor for parallel processing ###
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_single_vep_file, vep_txt_file, HGVS_dir, vep_csv_dir): vep_txt_file
            for vep_txt_file in vep_txt_files
        }

        for future in as_completed(futures):
            vep_txt_file = futures[future]
            try:
                future.result() 
            except Exception as e:
                print(f"Failed to process {vep_txt_file}: {e}")


######################################################################################################################################################################################################################
                                                                               ### ADD SAMPLE COUNT, DEPTH AND FREQUENCY TO VEP CSV ###                                                                               
######################################################################################################################################################################################################################

def append_maf_to_vep(vep_csv_dir, variants_per_sample_dir, vep_maf_dir):
    """
    Appends Total_Depth, Total_Count, and Total_Frequency values from {sample}_Ex{ex}_MAF.csv to {sample}_Ex{exon}_VEP.csv.

    Parameters:
    vep_csv_dir (str): Directory containing {sample}_Ex{ex}_VEP.csv files.
    variants_per_sample_dir (str): Directory containing {sample}_Ex{ex}_MAF.csv files.
    vep_maf_dir (str): Directory to store the resulting {sample}_Ex{ex}_VEP_MAF.csv files.
    """

    vep_csv_files = glob.glob(os.path.join(vep_csv_dir, '*.csv'))
    maf_files = glob.glob(os.path.join(variants_per_sample_dir, '*.csv'))

    for vep_csv_file in vep_csv_files:
        sample_exon = os.path.basename(vep_csv_file).split('_VEP.csv')[0]
        matching_maf_file = os.path.join(variants_per_sample_dir, f"{sample_exon}_MAF.csv")

        if matching_maf_file in maf_files:
            # print(f"Processing: {vep_csv_file} with {matching_maf_file}")

            ### read VEP and MAF CSV files of same sample_exon together ###
            vep_df = pd.read_csv(vep_csv_file, skip_blank_lines=False)
            maf_df = pd.read_csv(matching_maf_file)

            ### extract total_depth value (constant throughout the MAF file) ###
            total_depth = maf_df[f"total_depth_{sample_exon}"].iloc[0]

            ### identify chunk sizes in maf_df (i.e. variants per read in SAM file) ###
            chunk_sizes = []
            for variants in maf_df['cdna_variants']:
                num_variants = len(variants.split(', '))
                chunk_sizes.append(num_variants)

            ### add new columns for total count, depth and frequency to VEP dataframe ###
            vep_df["cdna_variants"] = pd.NA
            vep_df[f"Total_Depth_{sample_exon}"] = pd.NA
            vep_df[f"Total_Count_{sample_exon}"] = pd.NA
            vep_df[f"Total_Frequency_{sample_exon}"] = pd.NA

            ### identify chunks in vep_df (i.e. variants per read in SAM file) ###
            chunks = []
            current_chunk = []
            for index, row in vep_df.iterrows():
                if pd.isna(row['Uploaded_variation']):
                    if current_chunk:
                        chunks.append(current_chunk)
                        current_chunk = []
                else:
                    current_chunk.append(index)
            if current_chunk:
                chunks.append(current_chunk)

            ### ensure the number of chunks matches the number of rows in maf_df ###
            if len(chunks) != len(maf_df):
                print(f"Warning: Number of chunks in {vep_csv_file} does not match number of rows in {matching_maf_file}.")
                continue

            ### assign values to chunks ###
            for chunk_indices, (_, maf_row) in zip(chunks, maf_df.iterrows()):
                total_count = maf_row[f"total_count_{sample_exon}"]
                total_frequency = maf_row[f"total_frequency_{sample_exon}"]

                vep_df.loc[chunk_indices, "cdna_variants"] = maf_row['cdna_variants']
                vep_df.loc[chunk_indices, f"Total_Depth_{sample_exon}"] = total_depth
                vep_df.loc[chunk_indices, f"Total_Count_{sample_exon}"] = total_count
                vep_df.loc[chunk_indices, f"Total_Frequency_{sample_exon}"] = total_frequency

            ### reorder columns to make "cdna_variants" the first column ###
            columns_order = ["cdna_variants"] + [col for col in vep_df.columns if col != "cdna_variants"]
            vep_df = vep_df[columns_order]

            ### save the modified DataFrame to *_VEP_MAF.csv ###
            output_file = os.path.join(vep_maf_dir, f"{sample_exon}_VEP_MAF.csv")
            vep_df.to_csv(output_file, index=False, na_rep='', encoding='utf-8-sig')
            # print(f"Saved: {output_file}")
        else:
            print(f"No matching MAF file for {vep_csv_file}")


######################################################################################################################################################################################################################
                                                                          ### GET NUMBER OF MUTATIONS AND AA CHANGES PER VARIANT (COMBO) ###                                                                          
######################################################################################################################################################################################################################

### create dictionary to map protein position to original AA ###
def map_pp_to_aa(cDNA_file):
    cdna_file_df = pd.read_csv(cDNA_file)
    cdna_file_df = cdna_file_df.loc[cdna_file_df['reference'] == cdna_file_df['observed_nucleotide']]
    cdna_file_df = cdna_file_df.reset_index(drop=True)
    cdna_file_df['Original_AA'] = cdna_file_df['AA_Consequence'].str.split('.').str[1].str[0:1]
    cdna_file_df['Protein_position'] = cdna_file_df['AA_Consequence'].str.split('.').str[1].str[1:]
    cdna_file_df = cdna_file_df[["Protein_position", "Original_AA"]]
    cdna_file_df = cdna_file_df.drop_duplicates(subset=['Protein_position', 'Original_AA']).reset_index(drop=True)
    cdna_file_df['Protein_position'] = cdna_file_df['Protein_position'].astype(float).astype('Int64')
    pp_aa_map = cdna_file_df.set_index('Protein_position')['Original_AA'].to_dict()
    return pp_aa_map

def get_variant_aa_consequences(vep_maf_dir, variant_consequences_dir, pp_aa_map):
    """
    Add amino acid information to *_VEP_MAF.csv files 

    Parameters:
    vep_maf_dir (str): Directory containing *_VEP_MAF.csv files.
    variant_consequences_dir (str): Directory to save processed output files.
    pp_aa_map (dict): A dictionary mapping Protein_position to Original_AA.
    """

    ### find all *_VEP_MAF.csv files in the input directory ###
    maf_files = glob.glob(os.path.join(vep_maf_dir, '*_VEP_MAF.csv'))
    
    for maf_file in maf_files:
        ### Extract sample and exon information ###
        sample_exon = os.path.basename(maf_file).split('_VEP')[0]

        ### read the MAF file ###
        maf_df = pd.read_csv(maf_file, skip_blank_lines=False)

        ### process the file ###
        maf_df['Protein_position'] = maf_df['Protein_position'].astype(float).astype('Int64')
        maf_df['Original_AA'] = maf_df['Protein_position'].map(pp_aa_map)                                                                                                   ### use pp_aa_map to obtain reference amino acid of each protein position ###
        maf_df['Observed_AA'] = maf_df['Amino_acids'].apply(lambda x: x.split('/')[-1] if isinstance(x, str) and '/' in x else x if isinstance(x, str) else None)           ### determine observed amino acid ###
        maf_df["AA_change"] = maf_df["Original_AA"] + maf_df["Protein_position"].astype(str) + maf_df["Observed_AA"]                                                        ### construct AA_change string ###
        maf_df['AA_change'] = maf_df.apply(lambda row: '-' if row['Original_AA'] == row['Observed_AA'] else row['AA_change'], axis=1)

        maf_df = maf_df[["cdna_variants", "Uploaded_variation", "Consequence", "AA_change", f"Total_Depth_{sample_exon}", f"Total_Count_{sample_exon}", f"Total_Frequency_{sample_exon}"]]

        ### rename Total_Depth, Total_Count, and Total_Frequency columns ###
        maf_df = maf_df.rename(columns={f"Total_Depth_{sample_exon}": "Total_Depth", f"Total_Count_{sample_exon}": "Total_Count", f"Total_Frequency_{sample_exon}": "Total_Frequency"})

        ### reformat values in "Uploaded_variation" and "Consequence" columns for cleaner output ###
        maf_df["Uploaded_variation"] = maf_df["Uploaded_variation"].str.replace(transcript, "")                                     ### remove gencode transcript; retain only variant info ###
        maf_df["Consequence"] = maf_df["Consequence"].str.replace(",", ", ")                                                        ### split variant consequences with additional whitespace after comma ###

        ### join the values of other columns together into a string separated by ', ' for each unique 'cdna_variants' value ###
        maf_df = maf_df.groupby('cdna_variants').agg(
            {
                "Uploaded_variation": ", ".join,
                "Consequence": ", ".join,
                "AA_change": ", ".join,
                "Total_Depth": lambda x: ", ".join(map(str, x)),
                "Total_Count": lambda x: ", ".join(map(str, x)),
                "Total_Frequency": lambda x: ", ".join(map(str, x))
            }
        ).reset_index()
       
        ### remove duplicated values from each joined string e.g. [A, A, B, C] --> [A, B, C]
        maf_df["Consequence"] = maf_df["Consequence"].apply(lambda x: ", ".join(dict.fromkeys(x.split(", "))))
        maf_df["AA_change"] = maf_df["AA_change"].apply(lambda x: ", ".join(dict.fromkeys(x.split(", "))))
        maf_df["Total_Depth"] = maf_df["Total_Depth"].apply(lambda x: ", ".join(dict.fromkeys(x.split(", "))))
        maf_df["Total_Count"] = maf_df["Total_Count"].apply(lambda x: ", ".join(dict.fromkeys(x.split(", "))))
        maf_df["Total_Frequency"] = maf_df["Total_Frequency"].apply(lambda x: ", ".join(dict.fromkeys(x.split(", "))))

        ### clean up values in AA_change ###
        maf_df["AA_change"] = maf_df["AA_change"].apply(lambda x: x if x.strip() == "-" else x.replace("-, ", "").replace(", -", ""))

        ### calculate the number of mutations for each cdna_variant ###
        maf_df['Number_of_SNVs'] = maf_df['cdna_variants'].apply(lambda x: len(x.split(", ")))

        ### calculate the number of amino acid changes for each cdna_variant ###
        maf_df['Number_of_AA_changes'] = maf_df['AA_change'].apply(lambda x: len(x.split(", ")))
        maf_df['Number_of_AA_changes'] = maf_df.apply(lambda row: 0 if row['AA_change'] == "-" else row['Number_of_AA_changes'], axis=1)

        ### add back transcript to 'Uploaded_variation' values ###
        maf_df["Uploaded_variation"] = maf_df["Uploaded_variation"].apply(lambda x: f"{transcript}{x}")

        ### reorder columns ###
        maf_df = maf_df[["cdna_variants", "Uploaded_variation", "AA_change", "Total_Count", "Total_Depth", "Total_Frequency", "Number_of_SNVs", "Number_of_AA_changes", "Consequence"]]

        ### rename "Uploaded_variation" to "Variant(s)" ###
        maf_df = maf_df.rename(columns={"Uploaded_variation": "Variant(s)"})

        ### save the processed file as csv files in variant_consequences_dir ###
        output_file = os.path.join(variant_consequences_dir, f"{sample_exon}_variant_consequences.csv")
        maf_df.to_csv(output_file, index=False, na_rep='', encoding='utf-8-sig')
        # print(f"Processed and saved: {output_file}")


#####################################################################################################################################################################################################################
                                                                         ### COMPARE VARIANT COUNTS ACROSS DIFFERENT SAMPLE COMBINATIONS ###                                                                         
#####################################################################################################################################################################################################################

def create_variations_dict(variant_consequences_dir):
    """
    Converts *_variant_consequences.csv into dfs, combines them and stores only 'Variant(s)', 'AA_change', 'Number_of_SNVs', 'Number_of_AA_changes' & 'Consequence' information into variations_dict dictionary
    """
    csv_files = glob.glob(os.path.join(variant_consequences_dir, '*_variant_consequences.csv'))
    combined_df = pd.DataFrame()

    for file in csv_files:
        df = pd.read_csv(file)
        combined_df = pd.concat([combined_df, df[['cdna_variants', 'Variant(s)', 'AA_change', 'Number_of_SNVs', 'Number_of_AA_changes', 'Consequence']]], ignore_index=True)

    ### remove duplicates based on combined_df["Variation"] ###
    combined_df = combined_df.drop_duplicates(subset='Variant(s)')
    variations_dict = combined_df.set_index('Variant(s)').to_dict('index')
    return variations_dict

def combine_and_calculate_counts(list_df_primary, list_df_secondary, primary_label, secondary_label, variation_cols=['Variant(s)', 'cdna_variants']):
    """
    Combine primary and secondary DataFrames, group by both 'Variant(s)' and 'cdna_variants', calculate counts, and merge results.
    """

    ### combine primary dfs, rename columns and calculate sum ###
    df_primary_combined = pd.concat(list_df_primary).groupby(variation_cols, as_index=False).sum()
    df_primary_combined = df_primary_combined.rename(columns={'Total_Depth': f'Total_Depth_{primary_label}', 'Total_Count': f'Total_Count_{primary_label}', 'Total_Frequency': f'Total_Frequency_{primary_label}'})

    ### remove Total_Frequency column ###
    df_primary_combined = df_primary_combined.drop(columns=[f'Total_Frequency_{primary_label}'])
    df_primary_combined[f'Total_Frequency_{primary_label}'] = df_primary_combined[f'Total_Count_{primary_label}'] / df_primary_combined[f'Total_Depth_{primary_label}']
    
    ### combine secondary dfs, rename columns and calculate sum ###
    df_secondary_combined = pd.concat(list_df_secondary).groupby(variation_cols, as_index=False).sum()
    df_secondary_combined = df_secondary_combined.rename(columns={'Total_Depth': f'Total_Depth_{secondary_label}','Total_Count': f'Total_Count_{secondary_label}', 'Total_Frequency': f'Total_Frequency_{secondary_label}'})

    ### remove Total_Frequency column ###
    df_secondary_combined = df_secondary_combined.drop(columns=[f'Total_Frequency_{secondary_label}'])
    df_secondary_combined[f'Total_Frequency_{secondary_label}'] = df_secondary_combined[f'Total_Count_{secondary_label}'] / df_secondary_combined[f'Total_Depth_{secondary_label}']

    ### merge combined primary and combined secondary dfs ###
    merged_df = pd.merge(df_primary_combined, df_secondary_combined, on=variation_cols, how='outer')

    ### fill NaN values with zero for all relevant columns ###
    relevant_columns = [f'Total_Depth_{primary_label}', f'Total_Depth_{secondary_label}',
                        f'Total_Count_{primary_label}', f'Total_Count_{secondary_label}',
                        f'Total_Frequency_{primary_label}', f'Total_Frequency_{secondary_label}']
    for column in relevant_columns:
        merged_df[column] = merged_df[column].fillna(0)

    return merged_df

def process_combinations(sample_variant_comparisons_dir):
    """
    Processes combinations of data according to defined {sample}s to be compared
        
    Parameters:
        * sample_variant_comparisons_dir (str): sample_variant_comparisons_dir pre-defined globally to store all output files
    
    Workflow:
        * Extract and store {sample} from *_variant_consequences.csv filenames and in a dictionary
        * Data Combination and Analysis:
            * Calls the combine_and_calculate_counts function to combine data from primary and secondary dfs for each combination defined in the combinations list
            * Appends variation-specific data in "AA_change", "Number_of_SNVs", "Number_of_AA_changes" & "Consequence" columns to final result_df 
            * Save result_df to output {primary}_vs{secondary}_Ex{ex}.csv files
    """
    ### creates variations_dict where key is Variant(s) and values are 'AA_change', 'Number_of_SNVs', 'Number_of_AA_changes' & 'Consequence' ###
    variations_dict = create_variations_dict(variant_consequences_dir)

    ### list of sample_exon dataframes to process ###
    df_names = [f"{sample}_Ex{ex}" for sample in sample_names]

    ### generate all combinations to process according to group1 and group2 labels in config.yaml file ###
    combinations = []
    for combo in sample_combinations:
        group1 = [item + f'Ex{ex}' for item in combo['group1']]
        group2 = [item + f'Ex{ex}' for item in combo['group2']]
        group1_label = combo['group1_label']
        group2_label = combo['group2_label']
        combinations.append((group1, group2, group1_label, group2_label, ex))

    dfs = {sample_exon: pd.read_csv(os.path.join(variant_consequences_dir, f"{sample_exon}_variant_consequences.csv")) for sample_exon in df_names}
    for primaries, secondaries, primary_label, secondary_label, experiment in combinations:
        list_df_primary = [dfs[sample_exon] for sample_exon in primaries]
        list_df_secondary = [dfs[sample_exon] for sample_exon in secondaries]
        result_df = combine_and_calculate_counts(list_df_primary, list_df_secondary, primary_label, secondary_label)

        ### append variation-specific data from the dictionary only if not already in the df ###
        for column in ['AA_change', 'Number_of_SNVs', 'Number_of_AA_changes', 'Consequence']:
            if column not in result_df.columns:
                result_df[column] = result_df['Variant(s)'].apply(lambda x: variations_dict[x][column] if x in variations_dict else None)

        columns_order = ['cdna_variants', 'Variant(s)', 'AA_change', f'Total_Count_{primary_label}', f'Total_Depth_{primary_label}', f'Total_Frequency_{primary_label}',
                         f'Total_Count_{secondary_label}', f'Total_Depth_{secondary_label}', f'Total_Frequency_{secondary_label}', 'Number_of_SNVs', 'Number_of_AA_changes', 'Consequence']
        result_df = result_df[columns_order]
        for column in result_df.columns:
            if 'Total_Depth' in column:
                ### find the maximum total depth value ###
                max_value = result_df[result_df[column] > 0][column].max()
                ### replace all values in the column with the maximum total depth value found ###
                result_df[column] = max_value

        output_file = os.path.join(sample_variant_comparisons_dir, f"{primary_label}_vs_{secondary_label}_Ex{experiment}.csv")
        result_df.to_csv(output_file, index=False, float_format='%.9f')


######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
                                                                                  ### ADDITIONAL ANALYSIS: SNV VARIANT EXPANSION ###                                                                                
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################

def SNV_expansion(variant_consequences_dir, sample_SNV_counts_dir):
    """
    Expands compound SNV entries (MNVs) in 'cdna_variants' into individual SNVs,
    sums their Total_Count across all appearances in each input file, 
    and saves per-sample SNV-level counts and frequencies into separate CSVs.
    """
    ### list all CSV files in the input directory ###
    csv_files = [f for f in os.listdir(variant_consequences_dir) if f.endswith('.csv')]
    
    for csv_file in csv_files:
        ### read each csv file into df ###
        df = pd.read_csv(os.path.join(variant_consequences_dir, csv_file))
        
        ### extract the sample_exon from the filename ###
        sample_exon = csv_file.split('_variant_consequences.csv')[0]
        
        ### create a dictionary to store summed counts for each unique variant ###
        variant_counts = {}
        
        ### iterate over each row in the DataFrame ###
        for _, row in df.iterrows():
            ### split 'cdna_variants' by ', ' to get individual variants ###
            variants = row['cdna_variants'].split(', ')
            count = row['Total_Count']
            total_depth = row['Total_Depth']
            
            ### sum the count for each unique variant ###
            for variant in variants:
                if variant not in variant_counts:
                    variant_counts[variant] = count
                else:
                    variant_counts[variant] += count
        
        ### create a new DataFrame from the summed variant counts ###
        result_df = pd.DataFrame({
            'cdna_variants': list(variant_counts.keys()),
            'Total_Count': list(variant_counts.values()),
            'Total_Depth': total_depth  # total_depth is constant
        })
        
        ### rename columns ###
        result_df.rename(columns={'cdna_variants': 'cdna_SNV', 'Total_Count':f'Total_Count_{sample_exon}', 'Total_Depth':f'Total_Depth_{sample_exon}'}, inplace=True)

        ### calculate Total_Frequency ###
        result_df[f'Total_Frequency_{sample_exon}'] = result_df[f'Total_Count_{sample_exon}'] / result_df[f'Total_Depth_{sample_exon}']

        ### save to csv ###
        output_file = os.path.join(sample_SNV_counts_dir, f"{sample_exon}_SNV_counts.csv")
        result_df.to_csv(output_file, index=False)

def combine_SNV_sample(sample_SNV_counts_dir):
    """
    Combines all per-sample SNV count CSVs into a single comparison table, normalises depth columns, and outputs multiple summary files.
    """
    ### find all {prefix}_condensed.csv excluding {prefix}_SNV_filtered_condensed.csv files in the directory ###
    snv_files = glob.glob(os.path.join(sample_SNV_counts_dir, '*.csv'))

    ### combine all snv_files into a single DataFrame ###
    combined_df = pd.DataFrame()

    for file_path in snv_files:
        ### extract {sample}_Ex{ex_number} from the filename ###
        prefix = os.path.basename(file_path).split('_SNV_counts.csv')[0]

        ### read csv into df ###
        df = pd.read_csv(file_path)

        ### initialise with current df if combined_df is empty ###
        if combined_df.empty:
            combined_df = df

        ### merge with existing df if combined_df is not empty ###
        else:
            combined_df = pd.merge(combined_df, df, on='cdna_SNV', how='outer')

    ### fill missing values with 0 in the combined df ###
    combined_df.fillna(0, inplace=True)

    ### fill each Total_Depth column with its maximum value ###
    total_depth_cols = [col for col in combined_df.columns if col.startswith('Total_Depth_')]
    for col in total_depth_cols:
        combined_df[col] = combined_df[col].max()

    ### define the desired sample order and column types ###
    sample_order = [f"{sample}_Ex{ex}" for sample in sample_names]
    metric_order = ['Total_Count', 'Total_Depth', 'Total_Frequency']

    ### adjust the sorting function ###
    def custom_sort(column_name):
        ### ensure 'cdna_SNV' comes first ###
        if column_name == 'cdna_SNV':
            return (0, "", "")
        
        ### loop through each sample and metric to determine its position ###
        for sample_idx, sample in enumerate(sample_order, start=1):
            for metric_idx, metric in enumerate(metric_order, start=1):
                if column_name.startswith(f"{metric}_{sample}"):
                    return (sample_idx, metric_idx, sample)                                                                         ### sorting priority: sample --> metric ###

    ### apply the custom sorting function ###
    sorted_columns = sorted(combined_df.columns, key=custom_sort)

    ### reorder DataFrame columns ###
    combined_df = combined_df[sorted_columns]

    ### save to csv file ###
    output_file = os.path.join(sample_SNV_counts_dir, f"All_Samples_SNV_Counts_Ex{ex}.csv")
    combined_df.to_csv(output_file, index=False)

    ### create df with only 'Variant', 'Total_Count', 'Total_Depth' and 'Total_Frequency' columns for cleaner output ###
    freq_cols = [col for col in combined_df.columns if col.startswith('Total_Frequency')]
    count_cols = [col for col in combined_df.columns if col.startswith('Total_Count')]
    depth_cols = [col for col in combined_df.columns if col.startswith('Total_Depth')]
    combined_df = combined_df.drop(columns=freq_cols)
    combined_df.rename(columns={'cdna_SNV': 'Variant'}, inplace=True)
    combined_df["Total_Count"] = combined_df[count_cols].sum(axis=1)
    combined_df["Total_Depth"] = combined_df[depth_cols].sum(axis=1)
    combined_df["Total_Frequency"] = combined_df["Total_Count"] / combined_df["Total_Depth"]
    combined_df.columns = combined_df.columns.str.replace(f"_Ex{ex}", "")
    
    ### save as csv file ###
    condensed_output_file = os.path.join(sample_SNV_counts_dir, f"Condensed_All_Samples_SNV_Counts_Ex{ex}.csv")
    combined_df.to_csv(condensed_output_file, index=False)

    ### Sub Counts of 0 with 1 to avoid inf error in downstream stats analyses ###
    sub_1_df = combined_df.drop(columns=['Total_Count', 'Total_Depth', 'Total_Frequency'])
    sub_1_count_cols = [col for col in sub_1_df.columns if col.startswith('Total_Count')]
    sub_1_depth_cols = [col for col in sub_1_df.columns if col.startswith('Total_Depth')]
    sub_1_df[sub_1_count_cols] = sub_1_df[sub_1_count_cols].replace(0, 1)
    sub_1_df["Total_Count"] = sub_1_df[sub_1_count_cols].sum(axis=1)
    sub_1_df["Total_Depth"] = sub_1_df[sub_1_depth_cols].sum(axis=1)
    sub_1_df["Total_Frequency"] = sub_1_df["Total_Count"] / sub_1_df["Total_Depth"]

    ### save as csv file ###
    sub_1_output_file = os.path.join(sample_SNV_counts_dir, f"Sub_1_Condensed_All_Samples_SNV_Counts_Ex{ex}.csv")
    sub_1_df.to_csv(sub_1_output_file, index=False)


######################################################################################################################################################################################################################
######################################################################################################################################################################################################################
                                                                                      ### ADDITIONAL ANALYSIS: GET SNV PAIRS ###                                                                                    
######################################################################################################################################################################################################################
######################################################################################################################################################################################################################

def get_variant_pairs(variant_consequences_dir, sample_SNV_pairs_dir):
    """
    Processes per-sample variant consequence CSV files to extract co-occurring cDNA SNV pairs and calculates their joint counts and frequencies.
    """
    ### list all CSV files in the input directory ###
    csv_files = [f for f in os.listdir(variant_consequences_dir) if f.endswith('.csv')]
    
    for csv_file in csv_files:
        ### extract the sample_exon from the filename ###
        sample_exon = csv_file.split('_variant_consequences.csv')[0]

        ### read csv into df ###
        df = pd.read_csv(os.path.join(variant_consequences_dir, csv_file))
        
        ### initialise pair_counts dictionary to store variant pairs and their counts ###
        pair_counts = {}

        ### extract the column names for count and total reads dynamically ###
        count_col = [col for col in df.columns if col.startswith("Total_Count")][0]
        depth_col = [col for col in df.columns if col.startswith("Total_Depth")][0]
        
        depth = df[depth_col].iloc[0]                                                                                               ### assumes total_reads is constant across all rows ###

        ### iterate over each row in df ###
        for _, row in df.iterrows():
            ### split the cdna_variants by ', ' into individual variants ###
            variants = row['cdna_variants'].split(', ')
            count = row[count_col]

            ### parse '-' as variant 2 if there's only 1 variant ###
            if len(variants) == 1:
                pair = (variants[0], '-')
                if pair not in pair_counts:
                    pair_counts[pair] = count
                else:
                    pair_counts[pair] += count
            else:
                ### generate all unique pairs (combinations) of variants and count occurrences ###
                for var1, var2 in combinations(variants, 2):
                    ### ensure consistent order of variant pairs ###
                    pair = tuple(sorted((var1, var2)))
                    if pair not in pair_counts:
                        pair_counts[pair] = count
                    else:
                        pair_counts[pair] += count

        ### create df with 'Variant1', 'Variant2', 'Total_Count', and 'Total_Depth' columns ###
        result_data = {
            'Variant1': [pair[0] for pair in pair_counts.keys()],
            'Variant2': [pair[1] for pair in pair_counts.keys()],
            count_col: list(pair_counts.values()),
            depth_col: [depth] * len(pair_counts)
        }

        result_df = pd.DataFrame(result_data)

        ### remove rows with '-' in either 'Variant1' or 'Variant2' ###
        result_df = result_df[(result_df['Variant1'] != '-') & (result_df['Variant2'] != '-')]

        ### rename columns ###
        result_df.rename(columns={'Total_Count': f'Total_Count_{sample_exon}', 'Total_Depth': f'Total_Depth_{sample_exon}'}, inplace=True)

        ### add 'Total_Frequency' column ###
        result_df[f'Total_Frequency_{sample_exon}'] = result_df[f'Total_Count_{sample_exon}'] / result_df[f'Total_Depth_{sample_exon}']

        ### save as csv files ### 
        output_file = os.path.join(sample_SNV_pairs_dir, csv_file.replace('_variant_consequences.csv', '_SNV_pair.csv'))
        result_df.to_csv(output_file, index=False)

def all_sample_SNV_pairs(sample_SNV_pairs_dir):
    """
    Combines SNV pair data from all samples for a given exon into a single df, processes it into multiple views (full, condensed, grouped comparison), and saves each view as a CSV file.
    """
    ### finds all *_SNV_pair.csv files in the directory ###
    snv_files = glob.glob(os.path.join(sample_SNV_pairs_dir, '*_SNV_pair.csv'))

    combined_df = pd.DataFrame()

    for file_path in snv_files:
        ### extract sample and exon information from the filename ###
        file_name = os.path.basename(file_path)
        sample_exon = file_name.split('_SNV_pair.csv')[0]
        count_col = f'Total_Count_{sample_exon}'
        depth_col = f'Total_Depth_{sample_exon}'
        freq_col = f'Total_Frequency_{sample_exon}'

        ### read csv into df ###
        df = pd.read_csv(file_path)
        df = df[["Variant1", "Variant2", count_col, depth_col, freq_col]]

        ### ensure that variant pairs are consistent by sorting each pair ###
        df[['Variant1', 'Variant2']] = df.apply(lambda row: sorted([row['Variant1'], row['Variant2']]), axis=1, result_type='expand')

        ### merge with the combined df, filling missing pairs with 0 for counts and max value for depth ###
        if combined_df.empty:
            combined_df = df
        else:
            combined_df = pd.merge(
                combined_df, df,
                on=['Variant1', 'Variant2'],
                how='outer'
            )

    ### fll missing values in count columns with 0 and in depth columns with the max value for that column ###
    for col in combined_df.columns:
        if col.startswith("Total_Count"):
            combined_df[col].fillna(0, inplace=True)
        elif col.startswith("Total_Depth"):
            combined_df[col].fillna(combined_df[col].max(), inplace=True)

    ### generate desired order of columns ###
    base_columns = ["Variant1", "Variant2"]
    metric_order = ["Total_Count", "Total_Depth", "Total_Frequency"]
    sample_order = sample_names

    ordered_columns = base_columns + [
        f"{metric}_{sample}_Ex{ex}"
        for sample in sample_order
        for metric in metric_order
    ]

    ### reorder the DataFrame columns ###
    combined_df = combined_df[ordered_columns]

    count_cols = [col for col in combined_df.columns if col.startswith("Total_Count")]
    depth_cols = [col for col in combined_df.columns if col.startswith("Total_Depth")]
    freq_cols = [col for col in combined_df.columns if col.startswith("Total_Frequency")]

    ### fill empty freq cols with 0 ###
    combined_df[freq_cols] = combined_df[freq_cols].fillna(0)

    ### save to csv file ###
    combined_sample_file_path = os.path.join(sample_SNV_pairs_dir, f"All_Samples_SNV_Pairs_Ex{ex}.csv")
    combined_df.to_csv(combined_sample_file_path, index=False)

    ### make copy for condensed version ###
    condensed_df = combined_df.copy()

    ### drop depth and frequency columns for cleaner output ###
    condensed_df = condensed_df.drop(columns=depth_cols + freq_cols)

    ### rename count columns ###
    condensed_df.columns = condensed_df.columns.str.replace("Total_Count", "Count")
    condensed_df.columns = condensed_df.columns.str.replace(f"_Ex{ex}", "")
    condensed_df["Count_iPSC+A"] = condensed_df[[f"Count_iPSC", f"Count_A"]].sum(axis=1)
    condensed_df["Count_B+C"] = condensed_df[[f"Count_B", f"Count_C"]].sum(axis=1)
    condensed_df["Total_Count"] = condensed_df[[f"Count_iPSC+A", f"Count_B+C"]].sum(axis=1)

    ### save condensed file to csv ###
    condensed_sample_file_path = os.path.join(sample_SNV_pairs_dir, f"All_Samples_SNV_Pairs_Condensed_Ex{ex}.csv")
    condensed_df.to_csv(condensed_sample_file_path, index=False)

    ### make copy of the combined DataFrame for comparison between iPSC+A vs B+C ###
    compare_df = combined_df.copy()

    compare_df["Total_Counts_iPSC+A"] = combined_df[[f"Total_Count_iPSC_Ex{ex}", f"Total_Count_A_Ex{ex}"]].sum(axis=1)
    compare_df["Total_Depths_iPSC+A"] = combined_df[[f"Total_Depth_iPSC_Ex{ex}", f"Total_Depth_A_Ex{ex}"]].sum(axis=1)
    compare_df["Total_Frequencies_iPSC+A"] = compare_df["Total_Counts_iPSC+A"] / compare_df["Total_Depths_iPSC+A"]

    compare_df["Total_Counts_B+C"] = combined_df[[f"Total_Count_B_Ex{ex}", f"Total_Count_C_Ex{ex}"]].sum(axis=1)
    compare_df["Total_Depths_B+C"] = combined_df[[f"Total_Depth_B_Ex{ex}", f"Total_Depth_C_Ex{ex}"]].sum(axis=1)
    compare_df["Total_Frequencies_B+C"] = compare_df["Total_Counts_B+C"] / compare_df["Total_Depths_B+C"]

    compare_df = compare_df[["Variant1", "Variant2", "Total_Counts_iPSC+A", "Total_Depths_iPSC+A", "Total_Frequencies_iPSC+A", "Total_Counts_B+C", "Total_Depths_B+C", "Total_Frequencies_B+C"]]

    ### save combined file to csv ###
    compare_sample_file_path = os.path.join(sample_SNV_pairs_dir, f"iPSC+A_vs_B+C_SNV_Pairs_Ex{ex}.csv")
    compare_df.to_csv(compare_sample_file_path, index=False)


def main():
    start_time = time.time()

    process_cx_sam_files_parallel_start_time = time.time()
    process_cx_sam_files()
    process_cx_sam_files_parallel_end_time = time.time()
    print(f"Time elapsed for process_cx_sam_files_parallel: {process_cx_sam_files_parallel_end_time - process_cx_sam_files_parallel_start_time:.2f} seconds")

    read_counts_dict_start_time = time.time()
    read_counts_dict = get_read_counts_parallel(filtered_sam_gzipped_dir)
    read_counts_dict_end_time = time.time()
    print(f"Time elapsed for read_counts_dict: {read_counts_dict_end_time - read_counts_dict_start_time:.2f} seconds")

    process_sam_cs_variants_start_time = time.time()
    processed_dfs, removed_counts = process_sam_cs_variants(filtered_sam_gzipped_dir)
    process_sam_cs_variants_end_time = time.time()
    print(f"Time elapsed for process_sam_cs_variants: {process_sam_cs_variants_end_time - process_sam_cs_variants_start_time:.2f} seconds")

    get_variants_start_time = time.time()
    get_variants_parallel(raw_output_dir, variants_per_rep_dir, nm_histograms_dir, removed_counts, read_counts_dict)
    get_variants_end_time = time.time()
    print(f"Time elapsed for get_variants: {get_variants_end_time - get_variants_start_time:.2f} seconds")

    combine_variant_across_reps_start_time = time.time()
    combine_variant_across_reps_parallel(variants_per_rep_dir, variants_per_sample_dir)
    combine_variant_across_reps_end_time = time.time()
    print(f"Time elapsed for combine_variant_across_reps: {combine_variant_across_reps_end_time - combine_variant_across_reps_start_time:.2f} seconds")

    format_HGVS_variants_start_time = time.time()
    HGVS_prep(variants_per_sample_dir, cDNA_file, HGVS_dir)
    format_HGVS_variants_end_time = time.time()
    print(f"Time elapsed for format_HGVS_variants: {format_HGVS_variants_end_time - format_HGVS_variants_start_time:.2f} seconds")

    run_vep_start_time = time.time()
    run_vep_in_parallel(HGVS_dir, vep_txt_dir, vep_path, vep_options)
    run_vep_end_time = time.time()
    print(f"Time elapsed for run_vep: {run_vep_end_time - run_vep_start_time:.2f} seconds")

    process_vep_files_start_time = time.time()
    process_vep_files_in_parallel(vep_txt_dir, HGVS_dir, vep_csv_dir)
    process_vep_files_end_time = time.time()
    print(f"Time elapsed for process_vep_files: {process_vep_files_end_time - process_vep_files_start_time:.2f} seconds")

    append_maf_to_vep_start_time = time.time()
    append_maf_to_vep(vep_csv_dir, variants_per_sample_dir, vep_maf_dir)
    append_maf_to_vep_end_time = time.time()
    print(f"Time elapsed for append_maf_to_vep: {append_maf_to_vep_end_time - append_maf_to_vep_start_time:.2f} seconds")

    get_variant_aa_consequences_start_time = time.time()
    pp_aa_map = map_pp_to_aa(cDNA_file)
    get_variant_aa_consequences(vep_maf_dir, variant_consequences_dir, pp_aa_map)
    get_variant_aa_consequences_end_time = time.time()
    print(f"Time elapsed for get_variant_aa_consequences: {get_variant_aa_consequences_end_time - get_variant_aa_consequences_start_time:.2f} seconds")

    process_combinations_start_time = time.time()
    process_combinations(sample_variant_comparisons_dir)
    process_combinations_end_time = time.time()
    print(f"Time elapsed for process_combinations: {process_combinations_end_time - process_combinations_start_time:.2f} seconds")

    SNV_expansion_start_time = time.time()
    SNV_expansion(variant_consequences_dir, sample_SNV_counts_dir)
    SNV_expansion_end_time = time.time()
    print(f"Time elapsed for SNV_expansion: {SNV_expansion_end_time - SNV_expansion_start_time:.2f} seconds")

    combine_SNV_sample_start_time = time.time()
    combine_SNV_sample(sample_SNV_counts_dir)
    combine_SNV_sample_end_time = time.time()
    print(f"Time elapsed for combine_SNV_sample: {combine_SNV_sample_end_time - combine_SNV_sample_start_time:.2f} seconds")

    get_variant_pairs_start_time = time.time()
    get_variant_pairs(variant_consequences_dir, sample_SNV_pairs_dir)
    get_variant_pairs_end_time = time.time()
    print(f"Time elapsed for get_variant_pairs: {get_variant_pairs_end_time - get_variant_pairs_start_time:.2f} seconds")

    all_sample_SNV_pairs_start_time = time.time()
    all_sample_SNV_pairs(sample_SNV_pairs_dir)
    all_sample_SNV_pairs_end_time = time.time()
    print(f"Time elapsed for all_sample_SNV_pairs: {all_sample_SNV_pairs_end_time - all_sample_SNV_pairs_start_time:.2f} seconds")

    end_time = time.time()

    print(f"Total time elapsed: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
