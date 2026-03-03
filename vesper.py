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
from pathlib import Path
import seaborn as sns

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
start_bp = config["start_bp"]
end_bp = config["end_bp"]
baseq_threshold = config["baseq_threshold"]
mapq_threshold = config["mapq_threshold"]
seq_error = config["seq_error"]
read_length_filter = config["read_length_filter"]
ex = config["exon_number"]
samples_to_analyse = config["samples_to_analyse"]
mpileup_file = config["mpileup_file_path"]

### Define output directories ###
unzipped_sam_dir = os.path.join(root_output_dir, f'Ex{ex}_unzipped_sam')                                
filtered_sam_unzipped_dir = os.path.join(root_output_dir, f'Ex{ex}_unzipped_filtered_sam')                
filtered_sam_gzipped_dir = os.path.join(root_output_dir, f'Ex{ex}_filtered_sam_gz')                  

sample_variants_dir = os.path.join(root_output_dir, 'sample_rep_raw')                             
variants_per_rep_dir = os.path.join(root_output_dir, 'variants_per_rep')                     
nm_histograms_dir = os.path.join(root_output_dir, 'nm_histograms')      
sample_snv_counts_dir = os.path.join(root_output_dir, "sample_snv_counts")
vep_mpileup_merged_dir = os.path.join(root_output_dir, "vep_mpileup_merged")
sample_SNV_pairs_dir = os.path.join(root_output_dir, "sample_snv_pairs")                   

### Create directories ###
for directory in [unzipped_sam_dir, filtered_sam_unzipped_dir, filtered_sam_gzipped_dir, sample_variants_dir, variants_per_rep_dir, nm_histograms_dir]:
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
    raw_output_path = os.path.join(sample_variants_dir, raw_output_file)
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
    Count number of times each unique 'cdna_variants' value was found in each csv file of sample_variants_dir

    Parameters:
    - file_path: File path to each CSV file in sample_variants_dir
    - variants_per_rep_dir: Path to the directory where processed files will be saved.
    - removed_counts: Dictionary containing removed row counts with sample IDs as keys.
    - nm_histograms_dir: Path to the directory where NM histograms will be saved as PDF files.
    """
    filename = os.path.basename(file_path)
    sample_ID = filename.split('.csv')[0]                                                                                           ### get sample_ex_rep e.g. iPSC_Ex9_Rep1, A_Ex11_Rep5 ###

    ### read each CSV file from sample_variants_dir into a df and process ###
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

def get_variants_parallel(sample_variants_dir, variants_per_rep_dir, nm_histograms_dir, removed_counts, read_counts_dict):
    """
    Process CSV files in parallel using ProcessPoolExecutor.
    """
    csv_files = [os.path.join(sample_variants_dir, f) for f in os.listdir(sample_variants_dir) if f.endswith('.csv')]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(get_variants, file_path, variants_per_rep_dir, nm_histograms_dir, removed_counts, read_counts_dict): file_path for file_path in csv_files}

        for future in as_completed(futures):
            try:
                filename = future.result()
                # print(f"Processed {filename}")
            except Exception as e:
                print(f"Exception processing {futures[future]}: {e}")
                

def cdna_variants_per_sample(sample_variants_dir, col="cdna_variants"):
    """
    For files named {sample}_Ex{ex}_Rep{rep}.csv in sample_variants_dir:
      - group files by {sample}_Ex{ex}
      - read/concat all reps in each group
      - group by cdna_variants and count occurrences -> 'count'
      - write {sample}_Ex{ex}_cdna_variant_counts.csv to output_dir
    """
    sample_variants_dir = Path(sample_variants_dir)

    pat = re.compile(r'^(?P<sample>.+)_Ex(?P<ex>\d+)_Rep(?P<rep>\d+)\.csv$')
    buckets = {}

    for p in sorted(sample_variants_dir.glob("*.csv")):
        m = pat.match(p.name)
        if not m:
            continue
        key = f"{m['sample']}_Ex{m['ex']}"
        df = pd.read_csv(p)

        if col not in df.columns:
            raise ValueError(f"{p.name} is missing required column '{col}'")

        buckets.setdefault(key, []).append(df[[col]].copy())

    outputs = []
    for key, parts in buckets.items():
        if not parts:
            continue
        cat = pd.concat(parts, ignore_index=True)

        # drop NAs; remove this line if you want to count NaNs as a category
        cat = cat.dropna(subset=[col])

        counts = (
            cat.groupby(col, as_index=False)
               .size()
               .rename(columns={"size": "count"})
               .sort_values("count", ascending=False)
        )

        out_path = sample_variants_dir / f"{key}.csv"
        counts.to_csv(out_path, index=False)
        outputs.append(str(out_path))

    if not outputs:
        raise FileNotFoundError("No files matched pattern {sample}_Ex{ex}_Rep{rep}.csv")

    return outputs

def count_snvs(sample_variants_dir, sample_snv_counts_dir, pattern=f"*_Ex{ex}.csv", col="cdna_variants",count_col="count"):
    """
    For every CSV matching `pattern` in `input_dir`:
      - split comma-separated variant lists in `col`
      - explode to one row per variant
      - sum `count_col` per variant
      - write <stem>_flat.csv to `output_dir` (or next to inputs if None)
    Returns: list of written file paths.
    """
    sample_variants_dir = Path(sample_variants_dir)
    sample_snv_counts_dir = Path(sample_snv_counts_dir)
    sample_snv_counts_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(sample_variants_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {sample_variants_dir}")

    written = []
    for p in files:
        df = pd.read_csv(p)

        if col not in df.columns or count_col not in df.columns:
            raise ValueError(f"{p.name} missing required columns: {col!r} and/or {count_col!r}")

        # Ensure numeric counts
        df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(0)

        # Split "A, B" -> ["A","B"], explode, clean
        df[col] = df[col].astype(str).str.split(r"\s*,\s*")
        flat = df.explode(col)
        flat = flat[flat[col].notna() & (flat[col] != "")]
        flat[col] = flat[col].str.strip()

        out = (flat.groupby(col, as_index=False)[count_col]
                    .sum()
                    .sort_values(count_col, ascending=False))

        out_path = sample_snv_counts_dir / f"{p.stem}.csv"
        out['pos'] = out['cdna_variants'].str.extract(r'(?P<pos>-?\d+)').astype(int)
        out = out.sort_values('pos').reset_index(drop=True)
        out = out.drop(columns=['pos'])  # Drop pos column if not needed in output
        # rename 'cdna_variants' column to 'snv_variant'
        out = out.rename(columns={'cdna_variants': 'snv_variant'})
        out.to_csv(out_path, index=False)
        written.append(str(out_path))

    return written

def merge_vep_mpileup(mpileup_csv, sample_snv_counts_dir, vep_mpileup_merged_dir):
    sample_snv_counts_dir = Path(sample_snv_counts_dir)
    out_dir = Path(vep_mpileup_merged_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mp = pd.read_csv(mpileup_csv)
    mp.rename(columns={'Variant': 'snv_variant'}, inplace=True)

    for c in ['Total_Count_iPSC','Total_Depth_iPSC',
              'Total_Count_A','Total_Depth_A',
              'Total_Count_B','Total_Depth_B',
              'Total_Count_C','Total_Depth_C']:
        mp[c] = pd.to_numeric(mp[c], errors='coerce')

    # 2) merge into each per-sample file and plot count vs *_counts
    colmap_counts_only = {
        'A': 'Total_Count_A',
        'B': 'Total_Count_B',
        'C': 'Total_Count_C',
        'iPSC': 'Total_Count_iPSC',
    }

    written = []
    plotted = []

    for p in sorted(sample_snv_counts_dir.glob("*_Ex*.csv")):
        stem = p.stem       
        sample = stem.split('_Ex', 1)[0]
        if sample not in colmap_counts_only:
            continue

        df = pd.read_csv(p)
        if 'snv_variant' not in df.columns or 'count' not in df.columns:
            # must have both for merge and plotting
            continue

        ycol = colmap_counts_only[sample]
        sub = mp[['snv_variant', ycol]].copy()
        merged = df.merge(sub, on='snv_variant', how='left')

        # write merged
        out_csv = out_dir / p.name
        merged.to_csv(out_csv, index=False)
        written.append(str(out_csv))

        # plot: x = pooled count, y = *_counts_mpileup
        x = pd.to_numeric(merged['count'], errors='coerce')
        y = pd.to_numeric(merged[ycol], errors='coerce')
        m = x.notna() & y.notna()
        if m.any():
            lo = float(min(x[m].min(), y[m].min()))
            hi = float(max(x[m].max(), y[m].max()))
            plt.figure(figsize=(6,6))
            plt.scatter(x[m], y[m], s=16, alpha=0.7)
            plt.plot([lo, hi], [lo, hi], linestyle='--', linewidth=1)
            plt.xlabel('count')
            plt.ylabel(ycol)
            plt.title(f"{stem}: {ycol} vs count (n={m.sum()})")
            plt.tight_layout()
            out_png = out_dir / f"{stem}_{ycol}_vs_count.png"
            plt.savefig(out_png, dpi=200)
            plt.close()
            plotted.append(str(out_png))

def get_variant_pairs(sample_variants_dir, sample_SNV_pairs_dir):
    """
    Processes per-sample variant consequence CSV files to extract co-occurring cDNA SNV pairs and calculates their joint counts and frequencies.
    """
    # create output directory if it doesn't exist
    os.makedirs(sample_SNV_pairs_dir, exist_ok=True)
    ### list all CSV files in the input directory ###
    csv_files = [f for f in os.listdir(sample_variants_dir) if f.endswith(f'Ex{ex}.csv')]                     

    for csv_file in csv_files:
        ### extract the sample_exon from the filename ###
        sample_exon = csv_file.split('.csv')[0]

        ### read csv into df ###
        df = pd.read_csv(os.path.join(sample_variants_dir, csv_file))                                   
        
        ### initialise pair_counts dictionary to store variant pairs and their counts ###                    
        pair_counts = {}

        ### extract the column names for count and total reads dynamically ###                               
        count_col = [col for col in df.columns if col.startswith("count")][0]                                                                    ### assumes total_reads is constant across all rows ###

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
            count_col: list(pair_counts.values())
        }

        result_df = pd.DataFrame(result_data)

        ### rename columns ###
        result_df.rename(columns={'count': f'Paired_Count_{sample_exon}'}, inplace=True)

        ### save as csv files ### 
        output_file = os.path.join(sample_SNV_pairs_dir, csv_file.replace('.csv', '_SNV_pair.csv'))
        result_df.to_csv(output_file, index=False)


def annotate_snv_pairs(sample_SNV_pairs_dir, sample_snv_counts_dir):
    pairs_dir = Path(sample_SNV_pairs_dir)
    counts_dir = Path(sample_snv_counts_dir)

    # filename patterns
    pat_pairs  = re.compile(r'^(?P<sample>[^_]+)_Ex(?P<ex>\d+)_SNV_pair\.csv$')
    pat_counts = re.compile(r'^(?P<sample>[^_]+)_Ex(?P<ex>\d+)\.csv$')

    # preload counts files into a dict: key -> {variant: total_count}
    counts_maps = {}
    for p in counts_dir.glob("*.csv"):
        m = pat_counts.match(p.name)
        if not m:
            continue
        key = f"{m['sample']}_Ex{m['ex']}"
        df = pd.read_csv(p)
        if not {'snv_variant','count'}.issubset(df.columns):
            continue
        # normalize
        df = df[['snv_variant','count']].copy()
        df['snv_variant'] = df['snv_variant'].astype(str).str.strip()
        df['count'] = pd.to_numeric(df['count'], errors='coerce')
        counts_maps[key] = dict(df.dropna(subset=['snv_variant'])[['snv_variant','count']].values)

    written = []
    for p in sorted(pairs_dir.glob("*.csv")):
        m = pat_pairs.match(p.name)
        if not m:
            continue
        key = f"{m['sample']}_Ex{m['ex']}"
        if key not in counts_maps:
            print(f"Skip {p.name}: no matching counts file for {key}")
            continue
        var2count = counts_maps[key]

        df = pd.read_csv(p)

        # locate the paired-count column (e.g., "Paired_Count_{sample}_Ex{ex}")
        paired_cols = [c for c in df.columns if c.startswith("Paired_Count")]
        if not paired_cols:
            # fallback: assume the 3rd column is the paired count
            if df.shape[1] < 3:
                print(f"Skip {p.name}: cannot find Paired_Count column")
                continue
            paired_col = df.columns[2]
        else:
            paired_col = paired_cols[0]

        # make sure needed columns exist
        for col in ['Variant1','Variant2', paired_col]:
            if col not in df.columns:
                raise ValueError(f"{p.name} is missing required column '{col}'")

        # numeric paired counts
        df[paired_col] = pd.to_numeric(df[paired_col], errors='coerce')

        # helpers
        def pct_for_variant(variant, paired):
            if isinstance(variant, str) and variant.strip() == '-':
                return '-'  # as requested
            v = str(variant).strip() if pd.notna(variant) else None
            if not v or v not in var2count:
                return '-'  # no total count available
            total = var2count[v]
            if pd.isna(total) or total == 0:
                return '-'  # avoid div-by-zero / NaN
            return round(float(paired) / float(total) * 100, 2) if pd.notna(paired) else '-'

        df['%_rt_Variant1'] = [pct_for_variant(v1, pc) for v1, pc in zip(df['Variant1'], df[paired_col])]
        df['%_rt_Variant2'] = [pct_for_variant(v2, pc) for v2, pc in zip(df['Variant2'], df[paired_col])]

        out_path = pairs_dir / f"{m['sample']}_Ex{m['ex']}_SNV_pair_pct.csv"
        df.to_csv(out_path, index=False)
        written.append(str(out_path))

    if not written:
        raise FileNotFoundError("No matching *_SNV_pair.csv files processed.")
    return written

def calculate_mpileup_counts(mpileup_csv, sample_SNV_pairs_dir, pattern="*_SNV_pair_pct.csv"):
    pairs_dir = Path(sample_SNV_pairs_dir)
    sample_SNV_pairs_dir = Path(sample_SNV_pairs_dir)

    mp = pd.read_csv(mpileup_csv)

    for c in ['Total_Count_iPSC', 'Total_Count_A', 'Total_Count_B', 'Total_Count_C']:
        mp[c] = pd.to_numeric(mp[c], errors='coerce')

    mp = mp[['Variant', 'Total_Count_iPSC', 'Total_Count_A', 'Total_Count_B', 'Total_Count_C']]

    # Helper: pick the right counts column for a sample, with fallbacks
    def pick_counts_col(sample):
        preferred = {
            "A":    ["Total_Count_A", "AlphaActNeg_counts", "A_counts"],
            "B":    ["Total_Count_B", "normalBNP_counts",  "B_counts"],
            "C":    ["Total_Count_C", "highBNP_counts",    "C_counts"],
            "iPSC": ["Total_Count_iPSC", "iPSC_counts"],
        }
        for col in preferred.get(sample, []):
            if col in mp.columns:
                return col
        raise ValueError(f"No counts column in mpileup_csv for sample '{sample}'. "
                         f"Looked for: {preferred.get(sample, [])}")

    # Process each *_SNV_pair_pct.csv
    pat = re.compile(r"^(?P<sample>[^_]+)_Ex(?P<ex>\d+)_SNV_pair_pct\.csv$")
    written = []
    for p in sorted(pairs_dir.glob(pattern)):
        m = pat.match(p.name)
        if not m:
            continue
        sample = m["sample"]
        counts_col = pick_counts_col(sample)

        df = pd.read_csv(p)
        for need in ["Variant1", "Variant2"]:
            if need not in df.columns:
                raise ValueError(f"{p.name} missing required column '{need}'")

        # normalize variants and bring in total counts for Variant1/Variant2
        v1 = df["Variant1"].astype(str).str.strip()
        v2 = df["Variant2"].astype(str).str.strip()

        # Map variant -> total mpileup count for this sample
        var2count = dict(
            mp[["Variant", counts_col]].assign(
                **{counts_col: pd.to_numeric(mp[counts_col], errors="coerce")}
            ).dropna().values
        )

        total1 = v1.map(var2count)       # may be NaN if missing/'-'
        total2 = v2.map(var2count)

        # Percent columns: coerce to numeric (NaN if '-' or empty)
        pct1 = pd.to_numeric(df.get("%_rt_Variant1", np.nan), errors="coerce")
        pct2 = pd.to_numeric(df.get("%_rt_Variant2", np.nan), errors="coerce")

        # Masks where we can compute numbers (variant not '-', pct not NaN, total not NaN/zero)
        ok1 = (v1 != "-") & pct1.notna() & total1.notna() & (total1 != 0)
        ok2 = (v2 != "-") & pct2.notna() & total2.notna() & (total2 != 0)

        # Initialize with '-' strings, then fill numeric where ok
        df[f"calculated_Variant1_mpileup_count_{sample}"] = "-"
        df[f"calculated_Variant2_mpileup_count_{sample}"] = "-"

        df.loc[ok1, f"calculated_Variant1_mpileup_count_{sample}"] = (
            (pct1[ok1] / 100.0) * total1[ok1]
        ).round()

        df.loc[ok2, f"calculated_Variant2_mpileup_count_{sample}"] = (
            (pct2[ok2] / 100.0) * total2[ok2]
        ).round()

        v1 = pd.to_numeric(df[f'calculated_Variant1_mpileup_count_{sample}'], errors='coerce')
        v2 = pd.to_numeric(df[f'calculated_Variant2_mpileup_count_{sample}'], errors='coerce')

        out_path = sample_SNV_pairs_dir / f"{sample}_SNV_pair_stats.csv"
        df.to_csv(out_path, index=False)
        written.append(str(out_path))

    if not written:
        raise FileNotFoundError(f"No files matching {pattern} in {pairs_dir}")
    return written

def remove_files(sample_SNV_pairs_dir):
    """keep only *_SNV_pair_stats.csv files, remove others"""
    pairs_dir = Path(sample_SNV_pairs_dir)
    for p in pairs_dir.glob("*_SNV_pair.csv"):
        p.unlink()  # remove file
    for p in pairs_dir.glob("*_SNV_pair_pct.csv"):
        p.unlink()  # remove file
    for p in pairs_dir.glob("*_SNV_pair_stats.csv"):
        pass
    
def combine_snv_pair_stats(sample_paired_dir):
    """
    Reads {sample}_SNV_pair_stats.csv files (each with:
      Variant1, Variant2, calculated_Variant1_mpileup_count_<sample>,
      calculated_Variant2_mpileup_count_<sample> ),
    concatenates them, and returns a wide DataFrame with one row per (Variant1, Variant2)
    and two columns per sample:
      calculated_Variant1_mpileup_count_<sample>, calculated_Variant2_mpileup_count_<sample>.
    Missing combos per sample are filled with 0.
    """
    sample_paired_dir = Path(sample_paired_dir)
    pattern = "*_SNV_pair_stats.csv"
    files = sorted(sample_paired_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {sample_paired_dir}")

    # filename: "<sample>_SNV_pair_stats.csv"
    pat = re.compile(r"^(?P<sample>[^_]+)_SNV_pair_stats\.csv$")

    long_parts = []
    for p in files:
        m = pat.match(p.name)
        if not m:
            # skip files that don't match naming convention
            continue
        sample = m["sample"]
        # exclude samples not in the analysis list i.e. A
        if sample not in samples_to_analyse:
            continue

        df = pd.read_csv(p)
        if not {"Variant1", "Variant2"}.issubset(df.columns):
            raise ValueError(f"{p.name} must have 'Variant1' and 'Variant2' columns.")

        # Try to find the two calculated columns for this sample
        c1_default = f"calculated_Variant1_mpileup_count_{sample}"
        c2_default = f"calculated_Variant2_mpileup_count_{sample}"
        # fallback: pick first columns that start with the prefixes
        c1 = c1_default if c1_default in df.columns else next(
            (c for c in df.columns if c.startswith("calculated_Variant1_mpileup_count_")), None)
        c2 = c2_default if c2_default in df.columns else next(
            (c for c in df.columns if c.startswith("calculated_Variant2_mpileup_count_")), None)

        if c1 is None or c2 is None:
            raise ValueError(f"{p.name}: could not find calculated columns for sample '{sample}'")

        # coerce to numeric (turn '-' or blanks into 0)
        df[c1] = pd.to_numeric(df[c1], errors="coerce").fillna(0)
        df[c2] = pd.to_numeric(df[c2], errors="coerce").fillna(0)

        # build long rows for variant1 and variant2 calculated columns
        part1 = df[["Variant1", "Variant2", c1]].rename(columns={c1: "value"})
        part1["metric"] = "calculated_Variant1_mpileup_count"
        part1["sample"] = sample

        part2 = df[["Variant1", "Variant2", c2]].rename(columns={c2: "value"})
        part2["metric"] = "calculated_Variant2_mpileup_count"
        part2["sample"] = sample

        long_parts.append(part1)
        long_parts.append(part2)

    if not long_parts:
        raise RuntimeError("No valid *_SNV_pair_stats.csv files parsed.")

    long = pd.concat(long_parts, ignore_index=True)

    # Pivot to wide: one row per (Variant1, Variant2); columns = metric x sample
    wide = (
        long.pivot_table(
            index=["Variant1", "Variant2"],
            columns=["metric", "sample"],
            values="value",
            aggfunc="sum",   # sum duplicates within a file if any
            fill_value=0
        )
        .sort_index()
    )

    # Flatten MultiIndex columns to "<metric>_<sample>"
    wide.columns = [f"{m}_{s}" for m, s in wide.columns]
    wide = wide.reset_index()

    combined_csv = sample_paired_dir / "combined_snv_pair_stats.csv"
    # remove rows if all calculated* columns are 0
    wide = wide[(wide.filter(like="calculated_").sum(axis=1) != 0)]
    # save the wide DataFrame to CSV
    wide.to_csv(combined_csv, index=False)

    # choose only rows where both Variant1 and Variant2 are not '-'
    co_travelling_vars = wide[(wide['Variant1'] != '-') & (wide['Variant2'] != '-')]
    co_travelling_csv = sample_paired_dir / f"exon{ex}_co_travelling_snv_pairs.csv"
    co_travelling_vars.to_csv(co_travelling_csv, index=False)
    
    lone_travelling_vars = wide[(wide['Variant1'] == '-') | (wide['Variant2'] == '-')]
    lone_travelling_csv = sample_paired_dir / f"exon{ex}_lone_travelling_snv.csv"
    lone_travelling_vars.to_csv(lone_travelling_csv, index=False)



def main():
    start_time = time.time()

    process_cx_sam_files()
    read_counts_dict = get_read_counts_parallel(filtered_sam_gzipped_dir)
    processed_dfs, removed_counts = process_sam_cs_variants(filtered_sam_gzipped_dir)
    get_variants_parallel(sample_variants_dir, variants_per_rep_dir, nm_histograms_dir, removed_counts, read_counts_dict)
    cdna_variants_per_sample(sample_variants_dir)
    count_snvs(sample_variants_dir, sample_snv_counts_dir)
    merge_vep_mpileup(mpileup_file, sample_snv_counts_dir, vep_mpileup_merged_dir)
    get_variant_pairs(sample_variants_dir, sample_SNV_pairs_dir)
    annotated_files = annotate_snv_pairs(sample_SNV_pairs_dir, sample_snv_counts_dir)
    calculate_mpileup_counts(mpileup_file, sample_SNV_pairs_dir)
    remove_files(sample_SNV_pairs_dir)
    combined_df = combine_snv_pair_stats(sample_SNV_pairs_dir)

    end_time = time.time()
    print(f"Total time elapsed: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()