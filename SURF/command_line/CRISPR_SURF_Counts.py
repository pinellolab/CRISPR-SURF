#!/usr/bin/env python
import math
import gzip
import sys
import time
import logging
import numpy as np
import subprocess as sb

from CRISPR_SURF_Support_Functions import total_count_normalization, median_normalization, normalize_sgRNA_counts, str2bool

# # Time stamp analysis start
# timestr = time.strftime("%Y-%m-%d_%H.%M.%S")

# ##### Establish logging mechanism
# logger = logging.getLogger()
# logger.setLevel(logging.DEBUG)

# formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

# fh = logging.FileHandler(out_dir + '/crispr-surf.log')
# fh.setLevel(logging.DEBUG)
# fh.setFormatter(formatter)
# logger.addHandler(fh)

# ch = logging.StreamHandler()
# ch.setLevel(logging.DEBUG)
# ch.setFormatter(formatter)
# logger.addHandler(ch)

def crispr_surf_counts(sgRNA_library, t0_fastqs = None, tend_fastqs = None, normalization = 'median', counting_method = 'tracrRNA', tracrRNA = 'GTTTTAG', sgRNA_start = 0, sgRNA_stop = 20, count_minimum = 100, dropout_penalty = 'true', sgRNA_length = 20, nuclease = 'cas9', perturbation = 'indel', reverse_score = 'false', out_dir = '.', counted = 'false'):

	"""
	Function to convert pairs of FASTQs (pre- and post- enrichment) into sgRNA log2fc scores for each replicate.
	FASTQ files and sgRNA library file are required.

	if counted == False:
		sgRNA library file format: chr start stop sequence strand sgRNA_type

	elif counted == True:
		sgRNA library file format: chr start stop sequence strand sgRNA_type rep1_control_count rep2_control_count (... more control rep counts ...) rep1_sample_count rep2_sample_count (... more sample rep counts ...)

	Directly feeds as input into CRISPR-SURF analysis pipeline.
	"""

	# Store sgRNA library data
	sgRNA_library_data = {}
	negative_controls_list = []
	positive_controls_list = []
	t0_headers = []
	tend_headers = []
	sgRNA_counts = {}
	tmp_count = 0
	with open(sgRNA_library, 'r') as f:

		for line in f:

			line = line.strip().split(',')

			if not str2bool(counted):

				try:
					chrom, start_tmp, stop_tmp, sgRNA_sequence, strand, sgRNA_type = str(line[0]), str(line[1]), str(line[2]), str(line[3]), str(line[4]), str(line[5])

				except:
					sys.exit('ERROR: There needs to be 6 columns in the input sgRNA library file - chr start stop sequence strand sgRNA_type ...')

			else:

				if len(line) <= 7:
					sys.exit('ERROR: There needs to be at least 8 columns in the input sgRNA library file - chr start stop sequence strand sgRNA_type replicate1_control_count replicate1_sample_count ...')

				chrom, start_tmp, stop_tmp, sgRNA_sequence, strand, sgRNA_type = str(line[0]), str(line[1]), str(line[2]), str(line[3]), str(line[4]), str(line[5])

				try:
					sgRNA_counts[sgRNA_sequence] = [int(x) for x in line[6:]]

				except:
					sgRNA_counts[sgRNA_sequence] = [0 for x in range(len(line[6:]))]

				# Construct headers once
				if tmp_count == 0:

					replicates = int(len(line[6:])/2)

					t0_fastqs, tend_fastqs = list(range(replicates)), list(range(replicates))

					for i in range(1, replicates + 1):
						t0_headers.append('Replicate%s_Control_Count' % str(i))
						tend_headers.append('Replicate%s_Sample_Count' % str(i))

					tmp_count += 1

			if sgRNA_type == 'negative_control':
				negative_controls_list.append(sgRNA_sequence.upper())

			elif sgRNA_type == 'positive_control':
				positive_controls_list.append(sgRNA_sequence.upper())

			else:
				sgRNA_type == 'observation'

			try:
				start = min([int(start_tmp), int(stop_tmp)])
				stop = max([int(start_tmp), int(stop_tmp)])

				if perturbation == 'indel':

					if nuclease == 'cas9':

						if strand == '+':
							cleavage_pos = start + 17

						elif strand == '-':
							cleavage_pos = stop - 17

						else:
							cleavage_pos = int((start + stop)/2)

					elif nuclease == 'cpf1':

						if strand == '+':
							cleavage_pos = stop - 1

						elif strand == '-':
							cleavage_pos = start + 1

						else:
							cleavage_pos = int((start + stop)/2)

					else:
						cleavage_pos = int((start + stop)/2)

				elif perturbation == 'crispri':
					cleavage_pos = int((start + stop)/2)

				elif perturbation == 'crispra':
					cleavage_pos = int((start + stop)/2)

				else:
					cleavage_pos = int((start + stop)/2)

			except:
				start = start_tmp
				stop = stop_tmp
				cleavage_pos = 'NA'

			sgRNA_library_data[sgRNA_sequence] = [chrom, start, stop, cleavage_pos, sgRNA_sequence, strand, sgRNA_type]

	# Files to track problematic sgRNAs
	sgRNAs_missing_out = open(out_dir + '/sgRNAs_missing.csv', 'w')

	if not str2bool(counted):

		if (t0_fastqs == None) or (tend_fastqs == None):
			sys.exit('ERROR: The counted parameter is False and no FASTQs are given for either control or sample to perform counting procedure ...')

		# logger.info('Counting sgRNAs ...')

		# List control and sample FASTQs
		t0_fastqs = list(t0_fastqs)
		tend_fastqs = list(tend_fastqs)

		# Dictionary for sgRNA counts
		sgRNA_counts = {}
		count_stats = {}

		for i in range(len(t0_fastqs)):

			count_stats[t0_fastqs[i]] = [0, 0]
			count_stats[tend_fastqs[i]] = [0, 0]

		# Count control FASTQs
		for t0_n in range(len(t0_fastqs)):

			# logger.info('Counting sgRNAs in T0 FASTQ Replicate%s ...' % (t0_n + 1))

			t0_headers.append('Replicate%s_T0_Count' % str(int(t0_n) + 1))

			# Check if zipped
			if t0_fastqs[t0_n].endswith('.gz'):
				infile = gzip.open(t0_fastqs[t0_n], 'r')

			else:
				infile = open(t0_fastqs[t0_n] , 'r')

			# Loop through control FASTQ
			total_reads = 0

			while infile.readline():

				read_sequence = infile.readline().strip()
				infile.readline()
				infile.readline()

				total_reads += 1

				if counting_method == 'index':
					sgRNA_sequence = read_sequence[int(sgRNA_start):int(sgRNA_stop)]

				elif counting_method == 'tracrRNA':
					tracrRNA_index = read_sequence.rfind(tracrRNA)

					if tracrRNA_index != -1:

						if tracrRNA_index - sgRNA_length >= 0:
							sgRNA_sequence = read_sequence[tracrRNA_index - sgRNA_length:tracrRNA_index]

						else:
							sgRNA_sequence = 'G' + read_sequence[:tracrRNA_index]

				if sgRNA_sequence not in sgRNA_counts:
					sgRNA_counts[sgRNA_sequence] = [0]*(len(t0_fastqs) + len(tend_fastqs))

				sgRNA_counts[sgRNA_sequence][t0_n] += 1

			count_stats[t0_fastqs[t0_n]][1] = total_reads

		# Count sample FASTQs
		for tend_n in range(len(tend_fastqs)):

			# logger.info('Counting sgRNAs in Tend FASTQ Replicate%s ...' % (tend_n + 1))

			tend_headers.append('Replicate%s_Sample_Count' % str(int(tend_n) + 1))

			if tend_fastqs[tend_n].endswith('.gz'):
				infile = gzip.open(tend_fastqs[tend_n], 'r')

			else:
				infile = open(tend_fastqs[tend_n], 'r')

			# Loop through sample FASTQ
			total_reads = 0

			while infile.readline():

				read_sequence = infile.readline().strip()
				infile.readline()
				infile.readline()

				total_reads += 1

				if counting_method == 'index':
					sgRNA_sequence = read_sequence[int(sgRNA_start):int(sgRNA_stop)]

				elif counting_method == 'tracrRNA':
					tracrRNA_index = read_sequence.rfind(tracrRNA)

					if tracrRNA_index != -1:

						if tracrRNA_index - sgRNA_length >= 0:
							sgRNA_sequence = read_sequence[tracrRNA_index - sgRNA_length:tracrRNA_index]

						else:
							sgRNA_sequence = 'G' + read_sequence[:tracrRNA_index]

				if sgRNA_sequence in sgRNA_counts:
					sgRNA_counts[sgRNA_sequence][len(t0_fastqs) + tend_n] += 1

			count_stats[tend_fastqs[tend_n]][1] = total_reads

		# logger.info('Filtering sgRNAs ...')

		# Filter sgRNAs from provided library
		if len(t0_headers + tend_headers) > 0:
			sgRNAs_missing_out.write('sgRNA_Sequence' + ',' + ','.join(map(str, t0_headers + tend_headers)) + ',' + 'Reason' + '\n')

		else:
			sgRNAs_missing_out.write('sgRNA_Sequence' + ',' + 'Reason' + '\n')

		sgRNA_library_list = set(list(sgRNA_library_data.keys()) + list(negative_controls_list) + list(positive_controls_list))
		sgRNA_count_list = set(sgRNA_counts.keys())
		sgRNAs_missing_list = [x for x in sgRNA_library_list if x not in sgRNA_count_list]
		sgRNAs_delete_list = [x for x in sgRNA_count_list if x not in sgRNA_library_list]

		 # Output missing sgRNAs from provided list
		for sgRNA_sequence in sgRNAs_missing_list:
			del sgRNA_library_data[sgRNA_sequence]
			sgRNAs_missing_out.write(sgRNA_sequence + ',' + ','.join(map(str, [0]*(len(t0_fastqs)*2))) + ',' + 'The sgRNA is not present in FASTQ control files.' + '\n')

		# Delete sgRNAs that aren't in provided list
		for sgRNA_sequence in sgRNAs_delete_list:
			del sgRNA_counts[sgRNA_sequence]

	# Apply minimum sgRNA count for control sample
	replicate_replacement = {}
	usable_replicates = {}
	top_replicate = {}
	for sgRNA_sequence in list(sgRNA_counts):

		# Eliminate sgRNAs with homopolymer T greater than length 3
		if sum([1 if 'T'*x in sgRNA_sequence else 0 for x in range(4, len(sgRNA_sequence))]) > 0:
			sgRNAs_missing_out.write(sgRNA_sequence + ',' + ','.join(map(str, sgRNA_counts[sgRNA_sequence])) + ',' + 'This sgRNA contains a homopolymer T sequence greater than length 3.' + '\n')
			del sgRNA_counts[sgRNA_sequence]
			del sgRNA_library_data[sgRNA_sequence]

		else:
			# Impose minimum sgRNA count and dropout penalty
			sgRNA_control_count = [x for x in sgRNA_counts[sgRNA_sequence][:len(t0_fastqs)]]

			if str2bool(dropout_penalty):
				dropout_combined = [1 if x.count(0) > 0 else 0 for x in zip(sgRNA_counts[sgRNA_sequence][:len(t0_fastqs)], sgRNA_counts[sgRNA_sequence][len(t0_fastqs):])]
			else:
				dropout_combined = [0]*len(sgRNA_control_count)

			replicate_replacement[sgRNA_sequence] = [index for (index, count, dropout) in zip(range(len(sgRNA_control_count)), sgRNA_control_count, dropout_combined) if ((count < count_minimum) or (dropout == 1))]

			if len(t0_fastqs) <= 3:

				if len(replicate_replacement[sgRNA_sequence]) > 0:
					sgRNAs_missing_out.write(sgRNA_sequence + ',' + ','.join(map(str, sgRNA_counts[sgRNA_sequence])) + ',' + 'More than half of the replicates failed minimum control read count requirement and/or had sgRNA dropout.' + '\n')
					del sgRNA_counts[sgRNA_sequence]
					del sgRNA_library_data[sgRNA_sequence]
					del replicate_replacement[sgRNA_sequence]

				else:
					usable_replicates[sgRNA_sequence] = [index for (index, count, dropout) in zip(range(len(sgRNA_control_count)), sgRNA_control_count, dropout_combined) if ((count >= count_minimum) and (dropout == 0))]

			else:

				if len(replicate_replacement[sgRNA_sequence]) > int(math.ceil(len(t0_fastqs)/2.0)) - 1:
					sgRNAs_missing_out.write(sgRNA_sequence + ',' + ','.join(map(str, sgRNA_counts[sgRNA_sequence])) + ',' + 'More than half of the replicates failed minimum control read count requirement and/or had sgRNA dropout.' + '\n')
					del sgRNA_counts[sgRNA_sequence]
					del sgRNA_library_data[sgRNA_sequence]
					del replicate_replacement[sgRNA_sequence]

				else:
					usable_replicates[sgRNA_sequence] = [index for (index, count, dropout) in zip(range(len(sgRNA_control_count)), sgRNA_control_count, dropout_combined) if ((count >= count_minimum) and (dropout == 0))]

	sgRNAs_missing_out.close()

	# logger.info('Normalizing and computing sgRNA LFC scores ...')

	# Calculate counted reads in FASTQs

	if not str2bool(counted):

		for sgRNA_sequence in sgRNA_counts:

			for i in range(len(t0_fastqs)):
				count_stats[t0_fastqs[i]][0] += sgRNA_counts[sgRNA_sequence][i]
				count_stats[tend_fastqs[i]][0] += sgRNA_counts[sgRNA_sequence][len(t0_fastqs) + i]

	# Normalize sgRNA counts
	if normalization == 'median':
		sgRNA_norm_counts = normalize_sgRNA_counts(sgRNA_counts, method='median')

	elif normalization == 'total':
		sgRNA_norm_counts = normalize_sgRNA_counts(sgRNA_counts, method='total')

	elif normalization == 'none':
		sgRNA_norm_counts = normalize_sgRNA_counts(sgRNA_counts, method='none')

	# Compute log2FC score for FASTQ pairs
	sgRNA_lfc = {}
	for sgRNA_sequence in sgRNA_norm_counts:

		sgRNA_lfc[sgRNA_sequence] = []

		for i in range(len(t0_fastqs)):

			control_value = sgRNA_norm_counts[sgRNA_sequence][i]
			sample_value = sgRNA_norm_counts[sgRNA_sequence][len(t0_fastqs) + i]

			ratio = float(sample_value + 1.0)/float(control_value + 1.0)

			if str2bool(reverse_score):

				try:
					log2fc = -math.log(ratio, 2)

				except:
					log2fc = 'NaN'

			else:

				try:
					log2fc = math.log(ratio, 2)

				except:
					log2fc = 'NaN'

			sgRNA_lfc[sgRNA_sequence].append(log2fc)

	# Replace replicate lfc scores below count minimum and dropout using replicate with most supporting reads
	for sgRNA_sequence in replicate_replacement:
		for replicate in replicate_replacement[sgRNA_sequence]:
			sgRNA_lfc[sgRNA_sequence][replicate] = min([sgRNA_lfc[sgRNA_sequence][x] for x in usable_replicates[sgRNA_sequence]], key = abs)
			# sgRNA_lfc[sgRNA_sequence][replicate] = np.mean([sgRNA_lfc[sgRNA_sequence][x] for x in usable_replicates[sgRNA_sequence]])

	# logger.info('Outputting sgRNAs summary table ...')

	# Update sgRNA library dictionary with raw/normalized counts and LFC scores
	for sgRNA_sequence in sgRNA_counts:

		sgRNA_library_data[sgRNA_sequence] += sgRNA_counts[sgRNA_sequence]
		sgRNA_library_data[sgRNA_sequence] += sgRNA_norm_counts[sgRNA_sequence]
		sgRNA_library_data[sgRNA_sequence] += sgRNA_lfc[sgRNA_sequence]

	# Create header
	summary_header = ['Chr', 'Start', 'Stop', 'Perturbation_Index', 'sgRNA_Sequence', 'Strand', 'sgRNA_Type'] + t0_headers + tend_headers + [x.replace('Count','NormCount') for x in t0_headers] + [x.replace('Count','NormCount') for x in tend_headers]

	for i in range(len(t0_fastqs)):
		summary_header += ['Log2FC_Replicate%s' % (i + 1)]

	# Output summary sgRNA table
	with open(out_dir + '/sgRNAs_summary_table.csv', 'w') as f:

		for sgRNA_sequence in sgRNA_library_data:
			f.write(','.join(map(str, sgRNA_library_data[sgRNA_sequence])) + '\n')

	# Sort output bed for DBSCAN
	sb.call('sort -k1,1 -k2,2n %s > %s/tmp' % (out_dir + '/sgRNAs_summary_table.csv', out_dir), shell = True)
	sb.call('mv %s/tmp %s' % (out_dir, out_dir + '/sgRNAs_summary_table.csv'), shell = True)
	sb.call('echo %s | cat - %s > temp && mv temp %s' % (','.join(map(str, summary_header)), out_dir + '/sgRNAs_summary_table.csv', out_dir + '/sgRNAs_summary_table.csv'), shell = True)

	# Output FASTQ counting statistics
	if not str2bool(counted):
		with open(out_dir + '/FASTQ_count_statistics.csv', 'w') as f:
			f.write(','.join(map(str, ['FASTQ', 'Reads_Counted', 'Reads_Total', 'Percentage_Counted'])) + '\n')

			for fastq in count_stats:
				f.write(str(fastq.split('/')[-1]) + ',' + str(count_stats[fastq][0]) + ',' + str(count_stats[fastq][1]) + ',' + str(float(count_stats[fastq][0])/float(count_stats[fastq][1])) + '\n')
