#!/usr/bin/env python
import math
import gzip
import sys
import time
import logging
import numpy as np
import pandas as pd
import subprocess as sb

from CRISPR_SURF_Support_Functions import total_count_normalization, median_normalization, normalize_sgRNA_counts, str2bool

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def crispr_surf_counts(sgRNA_library, control_fastqs = None, sample_fastqs = None, nuclease = 'cas9', perturbation = 'indel', normalization = 'median', count_method = 'tracrRNA', tracrRNA = 'GTTTTAG', sgRNA_start = 0, sgRNA_stop = 20, count_minimum = 50, dropout_penalty = 'true', TTTT_penalty = 'true', sgRNA_length = 20, reverse_score = 'false', out_dir = '.'):

	"""
	Function to convert pairs of FASTQs (pre- and post- enrichment) into sgRNA log2fc scores for each replicate.
	FASTQ files and sgRNA library file are required.

	sgRNA library file format W/ FASTQ: chr start stop sequence strand sgRNA_type
	sgRNA library file format W/O FASTQ: chr start stop sequence strand sgRNA_type rep1_control_count rep2_control_count (... more control rep counts ...) rep1_sample_count rep2_sample_count (... more sample rep counts ...)

	Directly feeds as input into CRISPR-SURF analysis pipeline.
	"""

	# Read in sgRNA library file
	if '.csv' in sgRNA_library:
		df = pd.read_table(sgRNA_library, sep = ',', compression = 'infer')

	elif '.tsv' in sgRNA_library:
		df = pd.read_table(sgRNA_library, sep = '\t', compression = 'infer')

	else:
		logger.error("Please use either .csv or .tsv file format for the sgRNA library file ...")
		sys.exit("Please use either .csv or .tsv file format for the sgRNA library file ...")

	# Check necessary column names
	headers = ['chr', 'start', 'stop', 'sgrna_sequence', 'strand', 'sgrna_type']
	headers_out = ['Chr', 'Start', 'Stop', 'sgRNA_Sequence', 'Strand', 'sgRNA_Type']
	column_names = [x.lower() for x in df.columns]
	df.columns = column_names

	if all(x in column_names for x in headers):
		logger.info("Successfully read in sgRNA library file ...")

	else:
		logger.error("Necessary column names are missing (%s) ..." % ' '.join(map(str,headers)))
		sys.exit("Necessary column names are missing (%s) ..." % ' '.join(map(str,headers)))

	# Make sure all columns have correct formatting
	df['chr'] = df['chr'].str.lower()
	df['sgrna_sequence'] = df['sgrna_sequence'].str.upper()
	df['sgrna_type'] = df['sgrna_type'].str.lower()

	# Make sure sgRNA type is correct
	type_not_recognized = [x for x in df['sgrna_type'].unique() if x not in ['negative_control','positive_control', 'observation']]
	if len(type_not_recognized) > 0:
		logger.error("The sgRNA_Type column contains unknown annotations with %s. Only observation, negative_control, and positive_control values are allowed ..." % ','.join(map(str, type_not_recognized)))
		sys.exit("The sgRNA_Type column contains unknown annotations with %s. Only observation, negative_control, and positive_control values are allowed ..." % ','.join(map(str, type_not_recognized)))

	# Male sure start and stop columns are correct
	if np.issubdtype(df['start'].dtype, np.number) and np.issubdtype(df['stop'].dtype, np.number):
		df['start'] = df['start'].fillna(-1).astype(int).replace(-1, np.nan)
		df['stop'] = df['stop'].fillna(-1).astype(int).replace(-1, np.nan)

	else:
		logger.error("Please make sure all Start and Stop values are numeric or NaN (for non-targeting sgRNAs) ...")
		sys.exit("Please make sure all Start and Stop values are numeric or NaN (for non-targeting sgRNAs) ...")

	# Determine whether or not to count sgRNAs based on input sgRNA library file
	if len(df.columns) == 6:
		counted = False

	elif len(df.columns)%2 == 0:

		# Construct headers for new sgRNA summary table
		n_replicates = int((len(df.columns) - 6)/2)
		count_headers = []
		for n in range(1, n_replicates + 1):
			headers.append('replicate%s_control_count' % str(n))
			headers_out.append('Replicate%s_Control_Count' % str(n))
			count_headers.append('replicate%s_control_count' % str(n))

		for n in range(1, n_replicates + 1):
			headers.append('replicate%s_sample_count' % str(n))
			headers_out.append('Replicate%s_Sample_Count' % str(n))
			count_headers.append('replicate%s_sample_count' % str(n))

		if all(x in column_names for x in headers):

			for count_header in count_headers:
				if np.issubdtype(df[count_header].dtype, np.number):
					df[count_header] = df[count_header].fillna(0).astype(int)

				else:
					logger.error("The column %s contained non-numeric values. Please make sure all values are numeric ..." % count_header)
					sys.exit("The column %s contained non-numeric values. Please make sure all values are numeric ..." % count_header)

			sgRNA_counts = {}
			for index, row in df.iterrows():
				sgRNA_counts[row['sgrna_sequence']] = [row[x] for x in count_headers]

			counted = True

		else:
			logger.error("Necessary column names are missing assuming sgRNA counting is done (%s) ..." % ' '.join(map(str,headers)))
			sys.exit("Necessary column names are missing assuming sgRNA counting is done (%s) ..." % ' '.join(map(str,headers)))

	else:
		logger.error("There should be an even number of sgRNA count columns. Please modify the sgRNA library file to have even column number ...")
		sys.exit("There should be an even number of sgRNA count columns. Please modify the sgRNA library file to have even column number ...")

	# Add perturbation index column
	perturbation_indices = []
	if 'indel' in perturbation.lower():
		if 'cas9' in nuclease.lower():
			for index, row in df.iterrows():
				try:
					if row['strand'] == '+':
						perturbation_indices.append(int(min(row['start'], row['stop']) + 17))
					elif row['strand'] == '-':
						perturbation_indices.append(int(max(row['start'], row['stop']) - 17))
					else:
						perturbation_indices.append(int((row['start'] + row['stop'])/2.0))

				except:
					perturbation_indices.append(np.nan)


		elif 'cpf1' in nuclease.lower():
			for index, row in df.iterrows():
				try:
					if row['strand'] == '+':
						perturbation_indices.append(int(max(row['start'], row['stop']) - 1))
					elif row['strand'] == '-':
						perturbation_indices.append(int(min(row['start'], row['stop']) + 1))
					else:
						perturbation_indices.append(int((row['start'] + row['stop'])/2.0))

				except:
					perturbation_indices.append(np.nan)

		else:
			logger.info("Nuclease name not recognized (cas9 or cpf1). Continuing with perturbation index centered between sgRNA Start and Stop values ...")
			for index, row in df.iterrows():
				try:
					perturbation_indices.append(int((row['start'] + row['stop'])/2.0))

				except:
					perturbation_indices.append(np.nan)
	
	elif 'crispri' in perturbation.lower() or 'crispra' in perturbation.lower():
		for index, row in df.iterrows():
			try:
				perturbation_indices.append(int((row['start'] + row['stop'])/2.0))

			except:
				perturbation_indices.append(np.nan)
	elif 'be' in perturbation.lower():
		for index, row in df.iterrows():
			try:
				if row['strand'] == '+':
					perturbation_indices.append(int(min(row['start'], row['stop']) + 6))
				elif row['strand'] == '-':
					perturbation_indices.append(int(max(row['start'], row['stop']) - 6))
				else:
					perturbation_indices.append(int((row['start'] + row['stop'])/2.0))
			except:
				perturbation_indices.append(np.nan)
	else:
		logger.info("Perturbation type not recognized. Continuing with perturbation index centered between sgRNA Start and Stop values ...")		
		for index, row in df.iterrows():
			try:
				perturbation_indices.append(int((row['start'] + row['stop'])/2.0))

			except:
				perturbation_indices.append(np.nan)

	df['perturbation_index'] = perturbation_indices
	headers.append('perturbation_index')
	headers_out.append('Perturbation_Index')

	# Sort observation, negative, and positive control sgRNAs
	negative_controls_list = np.array(df.loc[df['sgrna_type'] == 'negative_control', 'sgrna_sequence']).flatten().tolist()
	positive_controls_list = np.array(df.loc[df['sgrna_type'] == 'positive_control', 'sgrna_sequence']).flatten().tolist()

	# Files to track problematic sgRNAs
	sgRNAs_missing_out = open(out_dir + '/sgRNAs_missing.csv', 'w')

	if not counted:

		if (control_fastqs == None) or (sample_fastqs == None):
			logger.error("FASTQ files for both control and sample are needed to perform counting process ...")
			sys.exit("FASTQ files for both control and sample are needed to perform counting process ...")

		# List control and sample FASTQs
		control_fastqs = list(control_fastqs)
		sample_fastqs = list(sample_fastqs)

		if len(control_fastqs) != len(sample_fastqs):
			logger.error("The number of control and sample FASTQs does not equal. If only a single control FASTQ is used, please just list it multiple times to match number of sample FASTQs ...")
			sys.exit("The number of control and sample FASTQs does not equal. If only a single control FASTQ is used, please just list it multiple times to match number of sample FASTQs ...")

		n_replicates = len(control_fastqs)

		# Dictionary for sgRNA counts
		sgRNA_counts = {}
		count_stats = {}

		for i in range(len(control_fastqs)):

			count_stats[control_fastqs[i]] = [0, 0]
			count_stats[sample_fastqs[i]] = [0, 0]

		# Count control FASTQs
		control_headers = []
		for n in range(len(control_fastqs)):

			headers.append('replicate%s_control_count' % str(n + 1))
			headers_out.append('Replicate%s_Control_Count' % str(n + 1))
			control_headers.append('replicate%s_control_count' % str(n + 1))
			logger.info('Starting to count sgRNAs for %s ...' % control_fastqs[n])

			# Check if zipped
			if control_fastqs[n].endswith('.gz'):
				infile = gzip.open(control_fastqs[n], 'r')

			else:
				infile = open(control_fastqs[n] , 'r')

			# Loop through control FASTQ
			total_reads = 0

			while infile.readline():

				read_sequence = infile.readline().strip()
				infile.readline()
				infile.readline()

				total_reads += 1

				if count_method.lower() == 'index':
					sgRNA_sequence = read_sequence[int(sgRNA_start):int(sgRNA_stop)].upper()

					if sgRNA_sequence not in sgRNA_counts:
						sgRNA_counts[sgRNA_sequence] = [0]*(len(control_fastqs) + len(sample_fastqs))

					sgRNA_counts[sgRNA_sequence][n] += 1

				elif count_method.lower() == 'tracrrna':
					tracrRNA_index = read_sequence.rfind(tracrRNA)

					if tracrRNA_index != -1:

						if tracrRNA_index - sgRNA_length >= 0:
							sgRNA_sequence = read_sequence[tracrRNA_index - sgRNA_length:tracrRNA_index].upper()

						# Rare case append a missing G (Fulco et al. 2016 data)
						else:
							sgRNA_sequence = 'G' + read_sequence[:tracrRNA_index].upper()

						if sgRNA_sequence not in sgRNA_counts:
							sgRNA_counts[sgRNA_sequence] = [0]*(len(control_fastqs) + len(sample_fastqs))

						sgRNA_counts[sgRNA_sequence][n] += 1

			count_stats[control_fastqs[n]][1] = total_reads

		# Count sample FASTQs
		sample_headers = []
		for n in range(len(sample_fastqs)):

			headers.append('replicate%s_sample_count' % str(n + 1))
			headers_out.append('Replicate%s_Sample_Count' % str(n + 1))
			sample_headers.append('replicate%s_sample_count' % str(n + 1))
			logger.info('Starting to count sgRNAs for %s ...' % sample_fastqs[n])

			if sample_fastqs[n].endswith('.gz'):
				infile = gzip.open(sample_fastqs[n], 'r')

			else:
				infile = open(sample_fastqs[n], 'r')

			# Loop through sample FASTQ
			total_reads = 0

			while infile.readline():

				read_sequence = infile.readline().strip()
				infile.readline()
				infile.readline()

				total_reads += 1

				if count_method.lower() == 'index':
					sgRNA_sequence = read_sequence[int(sgRNA_start):int(sgRNA_stop)].upper()

					if sgRNA_sequence in sgRNA_counts:
						sgRNA_counts[sgRNA_sequence][len(control_fastqs) + n] += 1

				elif count_method.lower() == 'tracrrna':
					tracrRNA_index = read_sequence.rfind(tracrRNA)

					if tracrRNA_index != -1:

						if tracrRNA_index - sgRNA_length >= 0:
							sgRNA_sequence = read_sequence[tracrRNA_index - sgRNA_length:tracrRNA_index].upper()

						# Rare case append a missing G (Fulco et al. 2016 data)
						else:
							sgRNA_sequence = 'G' + read_sequence[:tracrRNA_index].upper()

						if sgRNA_sequence in sgRNA_counts:
							sgRNA_counts[sgRNA_sequence][len(control_fastqs) + n] += 1

			count_stats[sample_fastqs[n]][1] = total_reads

		logger.info('Starting to filter sgRNAs based on provided sgRNA library file ...')

		# Filter sgRNAs from provided library
		sgRNAs_missing_out.write('sgRNA_Sequence' + ',' + ','.join(map(str, control_headers + sample_headers)) + ',' + 'Reason' + '\n')

		sgRNA_library_list = set(df['sgrna_sequence'].tolist())
		sgRNA_count_list = set(sgRNA_counts.keys())
		sgRNAs_missing_list = [x for x in sgRNA_library_list if x not in sgRNA_count_list]
		sgRNAs_delete_list = [x for x in sgRNA_count_list if x not in sgRNA_library_list]

		 # Output missing sgRNAs from provided list
		for sgRNA_sequence in sgRNAs_missing_list:
			df = df[df['sgrna_sequence'] != sgRNA_sequence]
			sgRNAs_missing_out.write(sgRNA_sequence + ',' + ','.join(map(str, [0]*(len(control_fastqs)*2))) + ',' + 'The sgRNA is not present in the FASTQ control files.' + '\n')

		# Delete sgRNAs that aren't in provided list
		for sgRNA_sequence in sgRNAs_delete_list:
			del sgRNA_counts[sgRNA_sequence]

	else:
		sgRNAs_missing_out.write('sgRNA_Sequence' + ',' + ','.join(map(str, count_headers)) + ',' + 'Reason' + '\n')

	# Apply minimum sgRNA count for control sample
	replicate_replacement = {}
	usable_replicates = {}
	top_replicate = {}

	# Eliminate sgRNAs with homopolymer T greater than length 3
	if str2bool(TTTT_penalty):
		for sgRNA_sequence in sgRNA_counts.keys():
				if sum([1 if 'T'*x in sgRNA_sequence else 0 for x in range(4, len(sgRNA_sequence))]) > 0:
					sgRNAs_missing_out.write(sgRNA_sequence + ',' + ','.join(map(str, sgRNA_counts[sgRNA_sequence])) + ',' + 'This sgRNA contains a homopolymer T sequence greater than length 3.' + '\n')
					del sgRNA_counts[sgRNA_sequence]
					df = df[df['sgrna_sequence'] != sgRNA_sequence]

	# Impose minimum sgRNA count and dropout penalty
	for sgRNA_sequence in sgRNA_counts.keys():

		sgRNA_control_count = [x for x in sgRNA_counts[sgRNA_sequence][:n_replicates]]

		if str2bool(dropout_penalty):
			dropout_combined = [1 if x.count(0) > 0 else 0 for x in zip(sgRNA_counts[sgRNA_sequence][:n_replicates], sgRNA_counts[sgRNA_sequence][n_replicates:])]
		else:
			dropout_combined = [0]*len(sgRNA_control_count)

		replicate_replacement[sgRNA_sequence] = [index for (index, count, dropout) in zip(range(len(sgRNA_control_count)), sgRNA_control_count, dropout_combined) if ((count < count_minimum) or (dropout == 1))]

		if n_replicates <= 3:

			if len(replicate_replacement[sgRNA_sequence]) > 0:
				sgRNAs_missing_out.write(sgRNA_sequence + ',' + ','.join(map(str, sgRNA_counts[sgRNA_sequence])) + ',' + 'More than half of the replicates failed minimum control read count requirement and/or had sgRNA dropout.' + '\n')
				del sgRNA_counts[sgRNA_sequence]
				df = df[df['sgrna_sequence'] != sgRNA_sequence]
				del replicate_replacement[sgRNA_sequence]

			else:
				usable_replicates[sgRNA_sequence] = [index for (index, count, dropout) in zip(range(len(sgRNA_control_count)), sgRNA_control_count, dropout_combined) if ((count >= count_minimum) and (dropout == 0))]

		else:

			if len(replicate_replacement[sgRNA_sequence]) > int(math.ceil(n_replicates/2.0)) - 1:
				sgRNAs_missing_out.write(sgRNA_sequence + ',' + ','.join(map(str, sgRNA_counts[sgRNA_sequence])) + ',' + 'More than half of the replicates failed minimum control read count requirement and/or had sgRNA dropout.' + '\n')
				del sgRNA_counts[sgRNA_sequence]
				df = df[df['sgrna_sequence'] != sgRNA_sequence]
				del replicate_replacement[sgRNA_sequence]

			else:
				usable_replicates[sgRNA_sequence] = [index for (index, count, dropout) in zip(range(len(sgRNA_control_count)), sgRNA_control_count, dropout_combined) if ((count >= count_minimum) and (dropout == 0))]

	sgRNAs_missing_out.close()

	logger.info('Normalizing and computing sgRNA L2FC enrichment scores ...')

	# Calculate counted reads in FASTQs
	if not counted:
		for sgRNA_sequence in sgRNA_counts.keys():
			for i in range(len(control_fastqs)):
				count_stats[control_fastqs[i]][0] += sgRNA_counts[sgRNA_sequence][i]
				count_stats[sample_fastqs[i]][0] += sgRNA_counts[sgRNA_sequence][len(control_fastqs) + i]

	# Normalize sgRNA counts

	if normalization == 'median':
		sgRNA_norm_counts = normalize_sgRNA_counts(sgRNA_counts, method='median')

	elif normalization == 'total':
		sgRNA_norm_counts = normalize_sgRNA_counts(sgRNA_counts, method='total')

	elif normalization == 'none':
		sgRNA_norm_counts = normalize_sgRNA_counts(sgRNA_counts, method='none')

	# Compute log2FC score for FASTQ pairs
	sgRNA_lfc = {}
	for sgRNA_sequence in sgRNA_norm_counts.keys():
		sgRNA_lfc[sgRNA_sequence] = []

		for i in range(n_replicates):
			control_value = sgRNA_norm_counts[sgRNA_sequence][i]
			sample_value = sgRNA_norm_counts[sgRNA_sequence][n_replicates + i]

			ratio = float(sample_value + 1.0)/float(control_value + 1.0)

			if str2bool(reverse_score):

				try:
					log2fc = -math.log(ratio, 2)

				except:
					log2fc = np.nan

			else:

				try:
					log2fc = math.log(ratio, 2)

				except:
					log2fc = np.nan

			sgRNA_lfc[sgRNA_sequence].append(log2fc)

	# Replace replicate lfc scores below count minimum and dropout using replicate with most conservative score
	for sgRNA_sequence in replicate_replacement:
		for replicate in replicate_replacement[sgRNA_sequence]:
			sgRNA_lfc[sgRNA_sequence][replicate] = min([sgRNA_lfc[sgRNA_sequence][x] for x in usable_replicates[sgRNA_sequence]], key = abs)

			# Using the mean of usable replicates
			# sgRNA_lfc[sgRNA_sequence][replicate] = np.mean([sgRNA_lfc[sgRNA_sequence][x] for x in usable_replicates[sgRNA_sequence]])

	logger.info('Outputting the sgRNAs summary table ...')

	#Update sgRNA library file with raw/normalized counts and L2FC enrichment scores
	sgRNA_order = df['sgrna_sequence'].tolist()
	if not counted:
		# Add raw, normalized and enrichment scores
		for n in range(len(control_headers)):

			# Raw
			df[control_headers[n]] = [sgRNA_counts[sgRNA_sequence][n] for sgRNA_sequence in sgRNA_order]
			df[control_headers[n].replace('count', 'normcount')] = [sgRNA_norm_counts[sgRNA_sequence][n] for sgRNA_sequence in sgRNA_order]

			# Norm
			df[sample_headers[n]] = [sgRNA_counts[sgRNA_sequence][len(control_headers) + n] for sgRNA_sequence in sgRNA_order]
			df[sample_headers[n].replace('count', 'normcount')] = [sgRNA_norm_counts[sgRNA_sequence][len(control_headers) + n] for sgRNA_sequence in sgRNA_order]

			# L2FC
			df['Log2FC_' + sample_headers[n].split('_')[0]] = [sgRNA_lfc[sgRNA_sequence][n] for sgRNA_sequence in sgRNA_order]

		# Update headers and headers_out
		for entry in control_headers:
			headers.append(entry.replace('count', 'normcount'))
			headers_out.append(entry.replace('replicate', 'Replicate').replace('control', 'Control').replace('count', 'NormCount'))

		for entry in sample_headers:
			headers.append(entry.replace('count', 'normcount'))
			headers_out.append(entry.replace('replicate', 'Replicate').replace('sample', 'Sample').replace('count', 'NormCount'))

		for entry in sample_headers:
			headers.append('Log2FC_' + entry.split('_')[0])
			headers_out.append('Log2FC_' + entry.split('_')[0].replace('replicate', 'Replicate'))

	else:
		# Add normalized and enrichment scores
		for n in range(len(count_headers)):

			# Norm
			df[count_headers[n].replace('count', 'normcount')] = [sgRNA_norm_counts[sgRNA_sequence][n] for sgRNA_sequence in sgRNA_order]

		for n in range(int(len(count_headers)/2)):
			# L2FC
			df['Log2FC_' + count_headers[n].split('_')[0]] = [sgRNA_lfc[sgRNA_sequence][n] for sgRNA_sequence in sgRNA_order]

		# Update headers and headers_out
		for entry in count_headers:
			headers.append(entry.replace('count', 'normcount'))
			headers_out.append(entry.replace('replicate', 'Replicate').replace('control', 'Control').replace('sample', 'Sample').replace('count', 'NormCount'))

		for entry in count_headers[:int(len(count_headers)/2)]:
			headers.append('Log2FC_' + entry.split('_')[0])
			headers_out.append('Log2FC_' + entry.split('_')[0].replace('replicate', 'Replicate'))

	# Re-order dataframe and add new columns
	df = df[headers]
	df.columns = headers_out

	df.to_csv(out_dir + '/sgRNAs_summary_table.csv', index = False, header = False)

	# Sort output sgRNAs summary table
	sb.call('sort -k1,1 -k2,2n %s > %s/tmp' % (out_dir + '/sgRNAs_summary_table.csv', out_dir), shell = True)
	sb.call('mv %s/tmp %s' % (out_dir, out_dir + '/sgRNAs_summary_table.csv'), shell = True)
	sb.call('echo %s | cat - %s > temp && mv temp %s' % (','.join(map(str, headers_out)), out_dir + '/sgRNAs_summary_table.csv', out_dir + '/sgRNAs_summary_table.csv'), shell = True)

	# Output FASTQ counting statistics
	if not counted:
		with open(out_dir + '/FASTQ_count_statistics.csv', 'w') as f:
			f.write(','.join(map(str, ['FASTQ', 'Reads_Counted', 'Reads_Total', 'Percentage_Counted'])) + '\n')

			for fastq in count_stats:
				f.write(str(fastq.split('/')[-1]) + ',' + str(count_stats[fastq][0]) + ',' + str(count_stats[fastq][1]) + ',' + str(float(count_stats[fastq][0])/float(count_stats[fastq][1])) + '\n')