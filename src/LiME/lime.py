import argparse
import sys
import time
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from pywaffle import Waffle
import re

EOL = '\n'

def SNV(vcf_file):
    # extract the SNV positions from a VCF file
    POS_list, QUAL_list, DP_list, DP_threshold = [], [], [], 0
    motif = re.compile(r'DP=(?P<deepth>[\d]+);')

    with open(vcf_file, 'r') as f:
        for line in f:
            if "#" not in line:
                POS_list.append(int(line.split("\t")[1]))
                QUAL_list.append(float(line.split("\t")[5]))
                match = motif.search(line.rstrip().split("\t")[7])
                DP_list.append(int(match['deepth']))
    
    data = {'POS':POS_list, 'QUAL':QUAL_list, 'DP':DP_list}
    df = pd.DataFrame(data)

    DP_threshold = 0.75* max(DP_list)
    df_filtered = df.query('QUAL>100 & DP>@DP_threshold')
    df_filtered = df_filtered.reset_index(drop=True)
    return df_filtered['POS'].tolist()

# Divide the positions given into small subgroups with 5 or less than 5 positions.
# ------
# Parameters:
# positions: list
#   list of positions given
# ------
# Returns:
# chunks: nested list
#   contains inner-lists of five positions.
#   If the last inner-list has only 1 position,
#   the last two inner-lists will be combined and re-split into 
#   two inner-lists with 3 positions in each of them.
def chunking(positions):
    i = 0
    chunks = []
    if len(positions) % 5 == 1:
        while i < (len(positions)-6):
            chunks.append(positions[i:(i + 5)])
            i += 5
        while i < len(positions):
            chunks.append(positions[i:(i + 3)])
            i += 3
    else:
        while i < len(positions):
            chunks.append(positions[i:(i + 5)])
            i += 5
    return chunks

# Collect the reads with specific nucleotides at each position and store this information in dictionary format.
# ------
# Parameters:
# Reads_input: pysam.libcalignmentfile.AlignmentFile
#   use BAM file using pysam
# positions: list
#   list of positions
# ------
# Returns:
# infotree: dictionary
#   key: str
#       position + nucleotide
#   value: list
#       list of read_ID harboring that "position + nucleotide"
def pos_nu_reads(Reads_input, positions):
    positions = [i-1 for i in positions]
    pos_i = 0
    nu_i = 0
    read_list = []
    infotree = {}
    nucleotides = ['A', 'T', 'G', 'C']
    number_of_target_positions = len(positions)

    while pos_i < number_of_target_positions:
        while nu_i < 4:
            for pileupcolumn in Reads_input.pileup():
                if pileupcolumn.pos == positions[pos_i]:
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if pileupread.alignment.query_sequence[pileupread.query_position] == nucleotides[nu_i]:
                                read_list.append(pileupread.alignment.query_name)
                    infotree[str(positions[pos_i]+1) + nucleotides[nu_i]] = read_list
            nu_i += 1
            read_list = []
        pos_i += 1
        nu_i = 0
    return(infotree)

# for each subgroups, collect all patterns detected
# return {patter: reads} and {read list included in this subgroup}.

# Get all possible nucleotide combination on the positions given (patterns), and collect all reads with each pattern
# ------
# Parameters:
# infotree: dictionary
#   key: str
#       position + nucleotide
#   value: list
#       list of read_ID harboring that "position + nucleotide"
# positions: list
#   list of positions
# ------
# Returns:
# intersection_tree: dictionary
#   key:str
#       pattern
#   value:
#       ID of all reads with that pattern
# o_readID: list
#   list of the ID of all the reads with any pattern detected for the positions of interest
def pattern_generation(infotree, positions):
    i = 0
    nucleotides = ['A', 'T', 'G', 'C']
    n = len(positions)
    p = [0] * n
    combination = [0] * n
    name_list = ['0'] * n
    name_list_i = 0
    intersection_name = ''
    intersection_tree = {}

    o_pattern = []
    o_readID = []

    while p[i] < 4:
        if i < (n - 1):
            combination[i] = p[i]
            p[i] += 1
            for x in range(i + 1, n):
                p[x] = 0
            i += 1
        else:
            combination[i] = p[i]  # get all possible nucleotide combination on the positions given (patterns)
            pos_i = 0
            intersection = set(infotree[str(positions[pos_i] + 1) + nucleotides[combination[pos_i]]]) & \
                           set(infotree[str(positions[pos_i + 1] + 1) + nucleotides[combination[pos_i + 1]]])
            name_list[pos_i] = str(positions[pos_i] + 1) + nucleotides[combination[pos_i]]
            name_list[pos_i + 1] = str(positions[pos_i + 1] + 1) + nucleotides[combination[pos_i + 1]]
            pos_i = 2
            while pos_i < n:
                intersection = intersection & set(infotree[str(positions[pos_i] + 1) + nucleotides[combination[pos_i]]])
                name_list[pos_i] = str(positions[pos_i] + 1) + nucleotides[combination[pos_i]]
                pos_i += 1        # for each pattern, collect the reads having this pattern
            if len(intersection) != 0:   # transform the list of 'position + nucleotide' into a string
                while name_list_i < n:
                    if name_list_i != (n - 1):
                        intersection_name = intersection_name + name_list[name_list_i] + ' '
                        name_list_i += 1
                    else:
                        intersection_name = intersection_name + name_list[name_list_i]
                        break
                intersection_tree[intersection_name] = intersection  # use the string as key to form a dictionary
                o_pattern.append(intersection_name)                  # store all the patterns in a list
                o_readID.append(intersection)                        # store all the reads in a nested list
                name_list_i = 0
                intersection_name = ''
            p[i] += 1
            while p[i] == 4 and i > 0:
                i -= 1
    return intersection_tree, list(o_readID)

# search for specific item in inner nested-dictionary

# return the the key of that inner dictionary
# store the keys found as a string

# Given the inner value in a nested dictionary, find the inner key corresponding to it.
# ------
# Parameters:
# nester_dictionary: dictionary
#   target nested dictionary
# target: str
#   the inner value
# ------
# Notes:
#   (In LiME only) the outer key is the positions in each subgroup, the outer value is the inner key and the inner value
#   the inner key is the pattern for those positions, the inner value is the ID of all the reads with that pattern.
def nested_dic_item_searching(nested_dictionary, target):
    pattern = ""
    for subgroup, subgroup_info in nested_dictionary.items():
        for key in subgroup_info:
            if target in subgroup_info[key]:
                pattern = pattern + " " + key
    return(pattern)
    
def lime(args):
    
    # Read the indexed BAM file as the input.
    test_dataset = pysam.AlignmentFile(args.i, "rb")

    # import positions of interest
    if args.vcf != None:
        pos_list = SNV(args.vcf)
    elif args.txt != None:
        with open(args.txt,'r') as f:
            pos_list = [int(i) for i in f.readline().rstrip().rsplit(', ')]
    else:
        print(f'Error: missing positions of interest.{EOL}')
        sys.exit()
        
    # LiME takes at least 2 positions.
    if len(pos_list) < 2:
        sys.exit(f"Error: The number of positions should be bigger than 1.{EOL}")
    start_time = time.time() # start timing
    
    chunks = chunking(pos_list)
    subgroup_read_tree = {}
    subgroup_info_tree = {}
    
    for item in chunks:
        subgroup_info_tree[str(item)], subgroup_read_tree[str(item)] = (
            pattern_generation(pos_nu_reads(test_dataset, item), item)
        )

    # get the list of reads shared in all subgroups
    shared_read = []
    if len(chunks) > 1:  # skip the loop while only one subgroup exits
        chunk_i = 0
        shared_read = set(subgroup_read_tree[str(chunks[chunk_i])]) & set(subgroup_read_tree[str(chunks[chunk_i + 1])])
        chunk_i = 2
        while chunk_i < len(chunks):
            shared_read = shared_read & set(subgroup_read_tree[str(chunks[chunk_i])])
            chunk_i += 1
    else:
        shared_read = set(subgroup_read_tree[str(chunks[0])])

    # Form big patterns using information from all subgroups
    # Generate a dictionary to store the information {unique_big_pattern: reads}
    pool = []
    pattern = ""
    pattern_tree = {}
    tem_list = []

    for item in shared_read:   # loop the reads with pattern detected
        pattern = nested_dic_item_searching(subgroup_info_tree, item)  # get the pattern in the read from the dictionary
        if pattern not in pool:  # create a dictionary with all unique pattern as key, reads with specific pattern as value
            pool.append(pattern)
            pattern_tree[pattern] = item
        else:
            tem_list.append(pattern_tree[pattern])
            tem_list.append(item)
            pattern_tree[pattern] = ",".join(tem_list)
            tem_list = []


    final_pattern_list = []
    final_read_list = []
    final_read_number = []
    for key, value in pattern_tree.items():  # get the columns to form DataFrame
        final_pattern_list.append(key)
        final_read_list.append(value.split(","))
        final_read_number.append(len(value.split(",")))

    data = {'Pattern': final_pattern_list,
            'Read_name': final_read_list,
            'Read_number': final_read_number}
    Pattern_Reads = pd.DataFrame(data, columns=['Pattern', 'Read_name', 'Read_number'])  # form the DataFrame
    if args.f:
        Pattern_Reads.to_csv(args.f, index=False)  # output the DataFrame in a CSV file

    # generate a waffle plot reflecting the pattern diversity in the input reads
    if args.p:
        fig = plt.figure(  
            FigureClass=Waffle,
            columns=60,
            values=Pattern_Reads.Read_number,
            cmap_name="tab20c"
            )

        fig.savefig(args.p, bbox_inches='tight')  

    finish_time = time.time()
    print(EOL.join([
        f"Program Finished in {float(finish_time - start_time)} seconds.",
        f'-------------'
        ])
         )