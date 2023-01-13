##### Command line script for running CRISPR-SURF #####
##### Author: Jonathan Y. Hsu (jyhsu@mit.edu) #####
##### Developed in the Pinello and Joung Labs at MGH/Harvard #####

# Import packages
import os
import sys
import time
import logging
import argparse
import numpy as np
import pandas as pd

# Import CRISPR-SURF functions
from CRISPR_SURF_Counts import crispr_surf_counts

# Time stamp analysis start
timestr = time.strftime("%Y-%m-%d_%H.%M.%S")

##### Argument handeling
parser = argparse.ArgumentParser(description = 'CRISPR-SURF command line tool to analyze CRISPR tiling screen data across various screening modalities.')
parser.add_argument('-f', '--sgRNA_library', type = str, required = True, help = 'Input sgRNA library file with the following required columns: Chr, Start, Stop, sgRNA_Sequence, Strand, sgRNA_Type. Additional sgRNA count columns required if no FASTQs provided: Replicate1_Control_Count, Replicate1_Sample_Count, etc.')
parser.add_argument('-control_fastqs', '--control_fastqs', default = None, nargs = '+', help = 'List of control FASTQs (before selection) for sgRNA counting. Example: Rep1_Control.fastq Rep2_Control.fastq Rep3_Control.fastq')
parser.add_argument('-sample_fastqs', '--sample_fastqs', default = None, nargs = '+', help = 'List of sample FASTQs (after selection) for sgRNA counting. Example: Rep1_Sample.fastq Rep2_Sample.fastq Rep3_Sample.fastq')
parser.add_argument('-nuclease', '--nuclease', default = 'cas9', type = str, choices = ['cas9', 'cpf1'], help = 'Nuclease used in the CRISPR tiling screen experiment. This information is used to determine the cleavage index if indels are specified as the perturbation.')
parser.add_argument('-pert', '--perturbation', default = 'indel', type = str, choices = ['indel', 'crispri', 'crispra', 'be'], help = 'Perturbation type used in the CRISPR tiling screen experiment. This information is used to determine the perturbation index for a given sgRNA.')
parser.add_argument('-norm', '--normalization', default = 'median', type = str, choices = ['none', 'total', 'median'], help = 'Normalization method between sequencing libraries.')
parser.add_argument('-count_method', '--count_method', default = 'tracrRNA', type = str, choices = ['tracrRNA', 'index'], help = 'Counting method for sgRNA counting from FASTQ. The tracrRNA option aligns a consensus sequence directly downstream of the sgRNA. The index option uses provided indices to grab sgRNA sequence from read.')
parser.add_argument('-tracrRNA', '--tracrRNA', type = str, default = 'GTTTTAG', help = 'If -count_method == tracrRNA. The consensus tracrRNA sequence directly downstream of the sgRNA for counting from FASTQ.')
parser.add_argument('-sgRNA_index', '--sgRNA_index', default = [0, 20], nargs = '+', help = 'If -count_method == index. The sgRNA start and stop indices (0-index) within the FASTQ reads. Example: 0 20')
parser.add_argument('-count_min', '--count_minimum', type = int, default = 50, help = 'The minimum number of counts for a given sgRNA in each control sample.')
parser.add_argument('-dropout', '--dropout_penalty', type = str, default = 'true', choices = ['t', 'true', 'yes', 'f', 'false', 'no'], help = 'The dropout penalty removes sgRNAs that have a 0 count in any of the control/sample replicates.')
parser.add_argument('-TTTT', '--TTTT_penalty', type = str, default = 'true', choices = ['t', 'true', 'yes', 'f', 'false', 'no'], help = 'The TTTT penalty removes sgRNAs that have a homopolymer stretch of Ts >= 4.')
parser.add_argument('-sgRNA_length', '--sgRNA_length', type = int, default = 20, help = 'Length of sgRNAs used in the CRISPR tiling screen experiment. This must match the sgRNA length provided in the sgRNA library file.')
parser.add_argument('-reverse', '--reverse_score', type = str, default = 'false', choices = ['t', 'true', 'yes', 'f', 'false', 'no'], help = 'Reverse the enrichment score. Generally applied to depletion screens where a positive score is associated with depletion of a sgRNA.')
parser.add_argument('-out_dir', '--out_directory', type = str, default = '.', help = 'The output directory for CRISPR-SURF counts.')
args = parser.parse_args()

##### Initialize arguments
sgRNA_library = args.sgRNA_library
control_fastqs = args.control_fastqs
sample_fastqs = args.sample_fastqs
nuclease = args.nuclease
perturbation = args.perturbation
normalization = args.normalization
count_method = args.count_method
tracrRNA = args.tracrRNA
sgRNA_index = args.sgRNA_index
count_minimum = args.count_minimum
dropout_penalty = args.dropout_penalty
TTTT_penalty = args.TTTT_penalty
sgRNA_length = args.sgRNA_length
reverse_score = args.reverse_score
out_dir = args.out_directory

##### Establish logging mechanism
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

fh = logging.FileHandler(out_dir + '/crispr-surf-counts.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

##### Run CRISPR-SURF Counts
if control_fastqs is not None:
	try:
		control_fastqs = [str(x) for x in control_fastqs]
	except:
		logger.error('The input for control_fastqs produced an error. Example: Rep1_Control.fastq Rep2_Control.fastq Rep3_Control.fastq')
		sys.exit('The input for control_fastqs produced an error. Example: Rep1_Control.fastq Rep2_Control.fastq Rep3_Control.fastq')

if sample_fastqs is not None:
	try:
		sample_fastqs = [str(x) for x in sample_fastqs]
	except:
		logger.error('The input for sample_fastqs produced an error. Example: Rep1_Sample.fastq Rep2_Sample.fastq Rep3_Sample.fastq')
		sys.exit('The input for sample_fastqs produced an error. Example: Rep1_Sample.fastq Rep2_Sample.fastq Rep3_Sample.fastq')

try:
	sgRNA_start = int(sgRNA_index[0])
	sgRNA_stop = int(sgRNA_index[1])
except:
	logger.error('The sgRNA_index input produced an error. Example: 0 20')
	sys.exit('The sgRNA_index input produced an error. Example: 0 20')

crispr_surf_counts(sgRNA_library = sgRNA_library, control_fastqs = control_fastqs, sample_fastqs = sample_fastqs, nuclease = nuclease, perturbation = perturbation, normalization = normalization, count_method = count_method, tracrRNA = tracrRNA, sgRNA_start = sgRNA_start, sgRNA_stop = sgRNA_stop, count_minimum = count_minimum, dropout_penalty = dropout_penalty, TTTT_penalty = TTTT_penalty, sgRNA_length = sgRNA_length, reverse_score = reverse_score, out_dir = out_dir)