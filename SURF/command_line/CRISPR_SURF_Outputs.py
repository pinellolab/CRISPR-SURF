#!/usr/bin/env python
import pandas as pd
import numpy as np
import math

def crispr_surf_sgRNA_summary_table_update(sgRNA_summary_table, gammas2betas, averaging_method, scale, guideindices2bin, simulation_n, padj_cutoffs, out_dir):

	"""
	Function to update sgRNA summary table with deconvolved signal and padj.-values
	If deconvolved signal doesn't have overlapping indices with sgRNA indices (scale > 1), nearest betas and padj.-values will be assigned.
	"""

	# Import sgRNA summary table
	df_summary_table = pd.read_csv(sgRNA_summary_table)
	replicates = len([x for x in df_summary_table.columns.tolist() if 'Log2FC_Replicate' in x])
	replicate_headers = [x for x in df_summary_table.columns.tolist() if 'Log2FC_Replicate' in x]
	gamma_chosen = gammas2betas['gamma_chosen']

	# Loop through sgRNA summary table and create deconvolved signal and padj.-value columns
	average_raw_scores = []
	replicate_deconvolutions = {}
	deconvolved_signal_combined = []
	p_values = []
	padj_values = []

	# Invert guideindices2bin
	if scale > 1:
		guideindices2bin_inverted = dict((v,k) for k in guideindices2bin for v in guideindices2bin[k])

	for index, row in df_summary_table.iterrows():

		if not math.isnan(row['Perturbation_Index']) and str(row['sgRNA_Type']) != 'negative_control':

			if scale > 1:
				closest_index = gammas2betas[1]['indices'].index(int(guideindices2bin_inverted[int(row['Perturbation_Index'])]))

			else:
				closest_index = gammas2betas[1]['indices'].index(int(row['Perturbation_Index']))

			for replicate in gammas2betas:
				if replicate != 'padj' and replicate != 'gamma_chosen' and replicate != 'combined' and replicate != 'indices'  and replicate != 'chr' and replicate != 'power' and replicate != 'p':

					if replicate not in replicate_deconvolutions:
						replicate_deconvolutions[replicate] = []

					replicate_deconvolutions[replicate].append(gammas2betas[replicate][gamma_chosen][closest_index])

			deconvolved_signal_combined.append(gammas2betas['combined'][closest_index])
			if gammas2betas['p'][closest_index] == 0:
				p_values.append('<%s' % (2.0/float(simulation_n)))
				padj_values.append('<%s' % min(padj_cutoffs))
			else:
				p_values.append(gammas2betas['p'][closest_index])
				padj_values.append(gammas2betas['padj'][closest_index])

			if averaging_method == 'mean':
				average_raw_scores.append(np.mean([row[x] for x in replicate_headers]))

			elif averaging_method == 'median':
				average_raw_scores.append(np.median([row[x] for x in replicate_headers]))

		else:

			for replicate in gammas2betas:
				if replicate != 'padj' and replicate != 'gamma_chosen' and replicate != 'combined' and replicate != 'indices'  and replicate != 'chr' and replicate != 'power' and replicate != 'p':

					if replicate not in replicate_deconvolutions:
						replicate_deconvolutions[replicate] = []

					replicate_deconvolutions[replicate].append('NaN')

			deconvolved_signal_combined.append('NaN')
			p_values.append('NaN')
			padj_values.append('NaN')

			if averaging_method == 'mean':
				average_raw_scores.append(np.mean([row[x] for x in replicate_headers]))

			elif averaging_method == 'median':
				average_raw_scores.append(np.median([row[x] for x in replicate_headers]))

	df_summary_table['Raw_Signal_Combined'] = average_raw_scores

	for replicate in replicate_deconvolutions:
		df_summary_table['Deconvolved_Signal_Replicate%s' % (str(replicate))] = replicate_deconvolutions[replicate]

	df_summary_table['Deconvolved_Signal_Combined'] = deconvolved_signal_combined
	df_summary_table['Pval.'] = p_values
	df_summary_table['Pval_adj.'] = padj_values

	# Output updated sgRNA summary table
	with open(out_dir + '/' + sgRNA_summary_table.replace('.csv', '_updated.csv'), 'w') as f:

		# Write columns
		f.write(','.join(map(str, df_summary_table.columns)) + '\n')

		# Write individual entries
		for index, row in df_summary_table.iterrows():

			f.write(','.join(map(str, row)) + '\n')

def complete_beta_profile(gammas2betas, simulation_n, padj_cutoffs, out_dir):
	"""
	Function to output total beta profile.
	"""

	chrom = gammas2betas['chr']
	indices = gammas2betas['indices']
	betas = gammas2betas['combined']
	pvals = gammas2betas['p']
	pvals_new = [x if x != 0 else '<%s' % (2.0/float(simulation_n)) for x in pvals]
	pvals_adj = gammas2betas['padj']
	pvals_adj_new = [x if x != 0 else '<%s' % (min(padj_cutoffs)) for x in pvals_adj]
	power = gammas2betas['power']

	df = pd.DataFrame({
		'Chr': chrom,
		'Index': indices,
		'Beta': betas,
		'Pval.': pvals_new,
		'Pval_adj.': pvals_adj_new,
		'Statistical_Power': power
		})

	df.to_csv(path_or_buf = (out_dir + '/beta_profile.csv'), index = False, header = True, columns = ['Chr','Index','Beta','Pval.','Pval_adj.','Statistical_Power'])

def crispr_surf_significant_regions(sgRNA_summary_table, gammas2betas, padj_cutoffs, scale, guideindices2bin, out_dir):

	"""
	Function to output significant regions text file.
	"""

	# Import sgRNA summary table
	df_summary_table = pd.read_csv(out_dir + '/' + sgRNA_summary_table)

	# Boundaries of inference
	diff_vec = [0] + list(np.diff(gammas2betas['indices']))

	boundaries = [0]
	for index, value in enumerate(diff_vec):
		if int(value) > int(scale):
			boundaries.append(index)

	boundaries.append(len(gammas2betas['indices']))

	with open(out_dir + '/significant_regions.csv', 'w') as f:

		# Header
		f.write(','.join(map(str, ['FDR', 'Chr', 'Start', 'Stop', 'Direction', 'Deconvolved_Signal_Area', 'Deconvolved_Signal_Mean', 'Padj_Mean', 'Supporting_sgRNAs', 'Supporting_sgRNA_Sequences'])) + '\n')

		for padj_cutoff in padj_cutoffs:

			# Find indices with significant padj-values
			significant_boundary_indices = []
			for i in range(len(boundaries) - 1):
				start_index, stop_index = boundaries[i], boundaries[i + 1]
				significant_indices = [1 if x < padj_cutoff else 0 for x in gammas2betas['padj'][start_index:stop_index]]
				significant_boundary_indices_tmp = [j + start_index for j, k in enumerate(np.diff([0] + significant_indices)) if k != 0]

				if len(significant_boundary_indices_tmp)%2 == 1:
					significant_boundary_indices_tmp.append(stop_index - 1)

				significant_boundary_indices += significant_boundary_indices_tmp

			# Invert guideindices2bin
			if scale > 1:
				guideindices2bin_inverted = dict((v,k) for k in guideindices2bin for v in guideindices2bin[k])

			for i, j in zip(significant_boundary_indices[0::2], significant_boundary_indices[1::2]):

				boundary_start = int(i)
				boundary_stop = int(j)

				genomic_boundary_start = gammas2betas['indices'][boundary_start]
				genomic_boundary_stop = gammas2betas['indices'][boundary_stop]

				# Find sgRNAs, padj, signal, and direction from significant region boundary indices
				# dff = df_summary_table.dropna().reset_index()

				chrom_index = gammas2betas['indices'].index(int(genomic_boundary_start))
				chrom = gammas2betas['chr'][chrom_index]

				associated_sgRNAs = np.array(df_summary_table.loc[(df_summary_table['Perturbation_Index'] >= genomic_boundary_start) & (df_summary_table['Perturbation_Index'] <= genomic_boundary_stop), 'sgRNA_Sequence']).flatten().tolist()
				padj_mean = np.mean(gammas2betas['padj'][boundary_start:boundary_stop])
				signal_mean = np.mean(gammas2betas['combined'][boundary_start:boundary_stop])
				signal_area = float(sum(gammas2betas['combined'][boundary_start:boundary_stop]))
				
				if signal_mean > np.median(gammas2betas['combined']):
					significance_direction = 'positive'

				else:
					significance_direction = 'negative'

				if len(associated_sgRNAs) > 0:
					f.write(','.join(map(str, [padj_cutoff, chrom, genomic_boundary_start, genomic_boundary_stop, significance_direction, signal_area, signal_mean, padj_mean, len(associated_sgRNAs), ','.join(map(str, associated_sgRNAs))])) + '\n')

def crispr_surf_IGV(sgRNA_summary_table, gammas2betas, padj_cutoffs, genome, scale, guideindices2bin, out_dir):

	"""
	Function to output tracks to load on IGV.
	Significant regions will be in .BED format.
	Raw and deconvolved scores will be in .BEDGRAPH format.
	"""

	# Import sgRNA summary table
	df_summary_table = pd.read_csv(out_dir + '/' + sgRNA_summary_table)
	replicates = len([x for x in df_summary_table.columns.tolist() if 'Log2FC_Replicate' in x])

	# Boundaries of inference
	diff_vec = [0] + list(np.diff(gammas2betas['indices']))

	boundaries = [0]
	for index, value in enumerate(diff_vec):
		if int(value) > int(scale):
			boundaries.append(index)

	boundaries.append(len(gammas2betas['indices']))

	with open(out_dir + '/positive_significant_regions.bed', 'w') as pos_bed, open(out_dir + '/negative_significant_regions.bed', 'w') as neg_bed:

		for padj_cutoff in padj_cutoffs:

			# Find indices with significant padj-values
			significant_boundary_indices = []
			for i in range(len(boundaries) - 1):
				start_index, stop_index = boundaries[i], boundaries[i + 1]
				significant_indices = [1 if x < padj_cutoff else 0 for x in gammas2betas['padj'][start_index:stop_index]]
				significant_boundary_indices_tmp = [j + start_index for j, k in enumerate(np.diff([0] + significant_indices)) if k != 0]

				if len(significant_boundary_indices_tmp)%2 == 1:
					significant_boundary_indices_tmp.append(stop_index - 1)

				significant_boundary_indices += significant_boundary_indices_tmp

			for i, j in zip(significant_boundary_indices[0::2], significant_boundary_indices[1::2]):

				boundary_start = int(i)
				boundary_stop = int(j)

				genomic_boundary_start = gammas2betas['indices'][boundary_start]
				genomic_boundary_stop = gammas2betas['indices'][boundary_stop]

				# Find sgRNAs, padj, signal, and direction from significant region boundary indices
				chrom_index = gammas2betas['indices'].index(int(genomic_boundary_start))
				chrom = gammas2betas['chr'][chrom_index]

				signal_mean = np.mean(gammas2betas['combined'][boundary_start:(boundary_stop + 1)])
				associated_sgRNAs = np.array(df_summary_table.loc[(df_summary_table['Perturbation_Index'] >= genomic_boundary_start) & (df_summary_table['Perturbation_Index'] <= genomic_boundary_stop), 'sgRNA_Sequence']).flatten().tolist()
				
				if len(associated_sgRNAs) > 0:
					if signal_mean > np.median(gammas2betas['combined']):
						pos_bed.write('\t'.join(map(str, [chrom, genomic_boundary_start, genomic_boundary_stop, padj_cutoff])) + '\n')

					else:
						neg_bed.write('\t'.join(map(str, [chrom, genomic_boundary_start, genomic_boundary_stop, padj_cutoff])) + '\n')

	# Output raw and deconvolved scores IGV track
	dff = df_summary_table[pd.notnull(df_summary_table['Chr']) & pd.notnull(df_summary_table['Perturbation_Index']) & pd.notnull(df_summary_table['Raw_Signal_Combined']) & pd.notnull(df_summary_table['Deconvolved_Signal_Combined'])]

	with open(out_dir + '/raw_scores.bedgraph', 'w') as raw_scores, open(out_dir + '/deconvolved_scores.bedgraph', 'w') as deconvolved_scores, open(out_dir + '/statistical_power.bedgraph', 'w') as statistical_power:

		for index, row in dff.iterrows():

			for r in range(1, replicates + 1):
				raw_scores.write('\t'.join(map(str, [row['Chr'], int(row['Perturbation_Index']), int(row['Perturbation_Index']), float(row['Log2FC_Replicate%s' % r]), row['sgRNA_Sequence']])) + '\n')

		for index in range(len(gammas2betas['indices'])):

			deconvolved_scores.write('\t'.join(map(str, [gammas2betas['chr'][index], int(gammas2betas['indices'][index]), int(gammas2betas['indices'][index]), float(gammas2betas['combined'][index])])) + '\n')
			statistical_power.write('\t'.join(map(str, [gammas2betas['chr'][index], int(gammas2betas['indices'][index]), int(gammas2betas['indices'][index]), float(gammas2betas['power'][index])])) + '\n')

	# Create IGV session
	with open('igv_session_template.xml', 'r') as f:
		igv_template = f.read()

	igv_template = igv_template.replace('#genome#', str(genome))

	with open(out_dir + '/igv_session.xml', 'w') as f:
		f.write(igv_template)