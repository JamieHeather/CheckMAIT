# CheckMAIT v2
# Jamie Heather, UCL, August 2016

# Searches through DCR/CDR3 data, looking for recombinations that in other settings have been recorded to belong to MAIT cells, or other T-cell subsets bearing invariant TCR alpha chains
# Requires as input a .dcrcdr3 file, i.e. the output of CDR3translator with the '-do' flag set to True

# Given the sparsity (and relative variability) of beta chain invariant sequences, this is only really recommended for alpha chain analyses

from __future__ import division 
import gzip
import argparse
import sys, os
import collections as coll
import os
import urllib2

# GET invariant SEQUENCES UP ONLINE
  # FIRST THING = IMPORT THEM IN

# FIX want option to check for a) CDR3 matches b) VJ matches and c) both matched (default)

def args():
    """args(): Obtains command line arguments which dictate the script's behaviour"""

    # Help flag
    parser = argparse.ArgumentParser(
        description='Identify rearranged TCR sequences in suitably demultiplexed FASTQ files, producing using the ligation TCRseq protocol.')
    # Add arguments
    parser.add_argument(
        '-in', '--infile', type=str, help='Name of .dcrcdr3 file to check for invariant sequences', required=True)
    parser.add_argument(
        '-o', '--output', type=bool, help='Output all detected invariant cell sequences to a separate output file', default=False)  
    parser.add_argument(
        '-dz', '--dontgzip', type=bool, help='Stop the output file files automatically being compressed with gzip (True/False)', required=False, default=False)
    parser.add_argument(
        '-pop', '--population', type=str, help='Specify which invariant TCR bearing population you wish to look for: \'mait\' (default), \'inkt\' or \'gem\'', required=False, default='mait')
    parser.add_argument(
        '-s', '--suppresssummary', type=bool, help='Output summary data (True/False)', required=False, default=False)    
    
    return parser.parse_args() 


def import_invariant_seqs():
    ''' Read in published invariant VJ-CDR3 sequences '''
    
    lines = []
    v_j_combs = []
    with open("alpha_"+ inputargs['population'].upper() +".vjcdr3") as infile:
      for line in infile:
        lines.append(tuple(line.rstrip().split(", ")))
        counts['invariant_sequences_in_file'] += 1
        v_j = tuple(line.rstrip().split(", ")[:2])
        if v_j not in v_j_combs:
          v_j_combs.append(v_j)
    return lines, v_j_combs


# GO THROUGH A GIVEN DCRCDR3 FILE
  # CHECK IF FILE SEQUENCE IS IN invariant LIST
  # IF SO, CAN DO 1 OR BOTH OF 2 OPTIONS
    # RECORD THE NUMBER OF DIFF AND CUMULATIVE PROP FREQ OF invariantS
    # OUTPUT THE ACTUAL invariant SEQUENCES TO NEW FILE
    

# OPTIONS TO LOOK FOR NKT/GEM



if __name__ == '__main__':
    
    inputargs = vars(args())
    # check file is alpha and .dcrcdr3 FIXXXXXXXXXXX
    
    counts = coll.Counter()
    
    # Read in invariant sequences and V-J combinations used to make them
    invariant_seqs, v_j_combs = import_invariant_seqs()

    # Read through infile and check for invariant cell sequences
    # FIX add check for presence of infile
    if inputargs['output'] == True:
      outfile_name = inputargs['infile'].split(".")[0] +'.'+ inputargs['population']
      outfile = open(outfile_name, 'w')

    if inputargs['infile'].endswith('.gz'):
      opener = gzip.open
    else:
      opener = open
      
    with opener(inputargs['infile']) as infile:
      for line in infile:
        
        bits = [x.split(", ") for x in line.rstrip().split(":")]
        v = bits[0][0]
        j = bits[0][1]
        cdr3 = bits[1][0]
        vj_cdr3 = tuple([v, j, cdr3])
        freq = int(bits[1][1])
        #sys.exit()
        
        counts['lines_in_unique'] += 1
        counts['lines_in_total'] += freq
        
        if tuple([v,j]) in v_j_combs:
          counts['v_j_match_unique'] += 1
          counts['v_j_match_total'] += freq         
        
          if vj_cdr3 in invariant_seqs:
            counts['invariant_line_unique'] += 1
            counts['invariant_line_total'] += freq
          
            if inputargs['output'] == True:
              outfile.write(line)
      
    if inputargs['output'] == True:
      outfile.close()
      if inputargs['dontgzip'] == False:
        print "Compressing output file..."
       
        with open(outfile_name) as infile, gzip.open(outfile_name + '.gz', 'wb') as zippedfile:
            zippedfile.writelines(infile)
        os.unlink(outfile_name)



# NEED TO ADD DECOMBINATOR OPTION TO SET URL FOR WHERE TO SEARCH FOR FILES?
 # AND FIX CHAIN BEHAViOUR!


# 


counts['pc_vjmatch_unique'] =  counts['v_j_match_unique'] / counts['lines_in_unique'] * 100
counts['pc_vjmatch_total'] = counts['v_j_match_total'] / counts['lines_in_total'] * 100

counts['pc_invariant_unique'] =  counts['invariant_line_unique'] / counts['lines_in_unique'] * 100
counts['pc_invariant_total'] = counts['invariant_line_total'] / counts['lines_in_total'] * 100


print counts['pc_invariant_unique'], "% all unique lines in detected as", inputargs['population'].upper()
print counts['pc_invariant_total'], "% all total lines in detected as", inputargs['population'].upper()




