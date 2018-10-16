#!/usr/bin/env python
pymin = min
pyabs = abs
# from cvxpy import *
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import laplace
import random
import logging
import multiprocessing as mp
from statsmodels.sandbox.stats.multicomp import multipletests

from CRISPR_SURF_Deconvolution import crispr_surf_deconvolution_simulations, crispr_surf_statistical_power

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def crispr_surf_deconvolved_signal(gammas2betas, gamma_chosen, averaging_method, out_dir):

	"""
	Function to construct final deconvolved signal based on chosen gamma.
	The gammas2betas dictionary has structure: Keys: 1) Replicate, 2) Gamma
	Averaging method can be mean or median to combine all biological replicates.
	"""

	# Create combined deconvolved signals from replicates
	deconvolved_signal = {}
	for i in [x for x in gammas2betas.keys() if ((x != 'combined') and (x != 'gamma_chosen') and (x != 'padj') and (x != 'indices') and (x != 'chr') and (x != 'p'))]:

		for j in range(len(gammas2betas[i][gamma_chosen])):

			if j not in deconvolved_signal:
				deconvolved_signal[j] = []

			deconvolved_signal[j].append(gammas2betas[i][gamma_chosen][j])

	# Create mean or median profile
	if averaging_method == 'mean':
		gammas2betas['combined'] = [np.mean(deconvolved_signal[x]) for x in deconvolved_signal]

	elif averaging_method == 'median':
		gammas2betas['combined'] = [np.median(deconvolved_signal[x]) for x in deconvolved_signal]

	gammas2betas['gamma_chosen'] = gamma_chosen
	gammas2betas['indices'] = gammas2betas[1]['indices']
	gammas2betas['chr'] = gammas2betas[1]['chr']

	chrom = gammas2betas['chr']
	indices = gammas2betas['indices']
	betas = gammas2betas['combined']

	# df = pd.DataFrame({
	# 	'Chr': chrom,
	# 	'Index': indices,
	# 	'Beta': betas
	# 	})

	# df.to_csv(path_or_buf = (out_dir + '/beta_profile.bedgraph'), index = False, header = False, columns = ['Chr','Index','Index','Beta'], sep = '\t')

	return gammas2betas

def crispr_surf_statistical_significance(sgRNA_summary_table, sgRNA_indices, perturbation_profile, gammas2betas, null_distribution, simulation_n, test_type, guideindices2bin, averaging_method, padj_cutoffs, effect_size, limit, scale):

	"""
	Function to assess the statistical significance of deconvolved genomic signal.
	Calculates empirical p-values for each beta, then performs FDR correction through the Benjamini-Hochberg procedure for p.adj.-values.
	"""

	# Load sgRNA summary table
	df_summary_table = pd.read_csv(sgRNA_summary_table)
	replicates = len([x for x in df_summary_table.columns.tolist() if 'Log2FC_Replicate' in x])

	# Gamma chosen for downstream analysis
	gamma_chosen = gammas2betas['gamma_chosen']

	# Load estimated betas into dictionary
	beta_distributions = {}
	for i in range(len(gammas2betas['combined'])):

		beta_distributions[i] = gammas2betas['combined'][i]

	# Decide how to draw from null distribution and perform deconvolution on simulated null arrays
	logger.info('Performing %s simulations to construct beta null distributions ...' % (simulation_n))

	if null_distribution == 'negative_control':
		if 'negative_control' not in df_summary_table['sgRNA_Type'].unique().tolist():
			null_distribution = 'gaussian'

	replicate_parameters = []
	if null_distribution == 'negative_control':

		replicate_parameters.append('NA')

		# Grab all negative control sgRNA lfc scores
		negative_control_guide_scores = []
		for i in range(1, int(replicates) + 1):
			negative_control_guide_scores.append(np.array(df_summary_table.loc[(df_summary_table['sgRNA_Type'] == 'negative_control'), ['Log2FC_Replicate' + str(i)]]).flatten().tolist())

		# Construct many simulated null arrays to perform deconvolution
		beta_distributions_null = crispr_surf_deconvolution_simulations(negative_control_scores = ['negative_control_guides', negative_control_guide_scores], sgRNA_indices = sgRNA_indices, perturbation_profile = perturbation_profile, gamma_list = [gamma_chosen], simulations_n = simulation_n, replicates = replicates, guideindices2bin = guideindices2bin, averaging_method = averaging_method, scale = scale)

	elif null_distribution == 'laplace':

		# Parameterize observed signal with laplace distribution (assume majority of observation sgRNAs are null)
		for i in range(1, int(replicates) + 1):

			# # Remove distribution skew
			# sorted_nc = sorted(np.array(df_summary_table.loc[(df_summary_table['sgRNA_Type'] != 'positive_control'), ['Log2FC_Replicate' + str(i)]]).flatten().tolist())
			# median_val = np.median(sorted_nc)
			# left_tail_median = np.median(sorted_nc[:int(0.1*len(sorted_nc))])
			# right_tail_median = np.median(sorted_nc[-int(0.1*len(sorted_nc)):])

			# # Left skewed
			# if (median_val - left_tail_median) > (right_tail_median - median_val):
			# 	half_dist = [x for x in sorted_nc if x >= median_val]

			# # Right skewed
			# else:
			# 	half_dist = [x for x in sorted_nc if x <= median_val]
			
			# half_dist_mirrored = [(2*median_val - x) for x in half_dist]
			# total_dist = half_dist + half_dist_mirrored
			# replicate_parameters.append(laplace.fit(total_dist))

			# Parameterize distribution directly
			observation_median = np.median(np.array(df_summary_table.loc[(df_summary_table['sgRNA_Type'] != 'positive_control'), ['Log2FC_Replicate' + str(i)]]).flatten().tolist())
			laplace_loc, laplace_scale = laplace.fit(np.array(df_summary_table.loc[(df_summary_table['sgRNA_Type'] != 'positive_control'), ['Log2FC_Replicate' + str(i)]]).flatten().tolist())
			replicate_parameters.append([observation_median, laplace_scale])

		# Construct many simulated null arrays to perform deconvolution
		beta_distributions_null = crispr_surf_deconvolution_simulations(negative_control_scores = ['laplace', replicate_parameters], sgRNA_indices = sgRNA_indices, perturbation_profile = perturbation_profile, gamma_list = [gamma_chosen], simulations_n = simulation_n, replicates = replicates, guideindices2bin = guideindices2bin, averaging_method = averaging_method, scale = scale)

	elif null_distribution == 'gaussian':

		# Parameterize observed signal with gaussian distribution (assume majority of observation sgRNAs are null)
		for i in range(1, int(replicates) + 1):

			# # Remove distribution skew
			# sorted_nc = sorted(np.array(df_summary_table.loc[(df_summary_table['sgRNA_Type'] != 'positive_control'), ['Log2FC_Replicate' + str(i)]]).flatten().tolist())
			# median_val = np.median(sorted_nc)
			# left_tail_median = np.median(sorted_nc[:int(0.1*len(sorted_nc))])
			# right_tail_median = np.median(sorted_nc[-int(0.1*len(sorted_nc)):])

			# # Left skewed
			# if (median_val - left_tail_median) > (right_tail_median - median_val):
			# 	half_dist = [x for x in sorted_nc if x >= median_val]

			# # Right skewed
			# else:
			# 	half_dist = [x for x in sorted_nc if x <= median_val]
			
			# half_dist_mirrored = [(2*median_val - x) for x in half_dist]
			# total_dist = half_dist + half_dist_mirrored
			# replicate_parameters.append(norm.fit(total_dist))

			# Parameterize distribution directly
			observation_median = np.median(np.array(df_summary_table.loc[(df_summary_table['sgRNA_Type'] != 'positive_control'), ['Log2FC_Replicate' + str(i)]]).flatten().tolist())
			gaussian_loc, gaussian_scale = norm.fit(np.array(df_summary_table.loc[(df_summary_table['sgRNA_Type'] != 'positive_control'), ['Log2FC_Replicate' + str(i)]]).flatten().tolist())
			replicate_parameters.append([observation_median, gaussian_scale])

		# Construct many simulated null arrays to perform deconvolution
		beta_distributions_null = crispr_surf_deconvolution_simulations(negative_control_scores = ['gaussian', replicate_parameters], sgRNA_indices = sgRNA_indices, perturbation_profile = perturbation_profile, gamma_list = [gamma_chosen], simulations_n = simulation_n, replicates = replicates, guideindices2bin = guideindices2bin, averaging_method = averaging_method, scale = scale)

	# Calculate p-values
	logger.info('Calculating p. values for %s betas ...' % (len(beta_distributions)))

	beta_pvals = []
	if test_type == 'nonparametric':

		for i in range(len(beta_distributions)):

			if (i + 1)%500 == 0:
				logger.info('Calculated p. values for %s out of %s betas ...' % ((i + 1), len(beta_distributions)))

			estimated_beta = beta_distributions[i]
			# null_betas = beta_distributions_null[i]
			# beta_pvals.append(2.0*float(max(0.0, min(sum(x >= estimated_beta for x in null_betas), sum(x <= estimated_beta for x in null_betas))))/float(len(null_betas)))

			null_betas = np.array(beta_distributions_null[i])
			beta_pvals.append(2.0 * min((null_betas >= estimated_beta).sum(), (null_betas <= estimated_beta).sum()) / float(len(null_betas)))

	elif test_type == 'parametric':

		for i in range(len(beta_distributions)):

			if (i + 1)%500 == 0:
				logger.info('Calculated p. values for %s out of %s betas ...' % ((i + 1), len(beta_distributions)))

			estimated_beta = beta_distributions[i]
			null_betas_loc, null_betas_scale = norm.fit(beta_distributions_null[i])
			
			beta_pvals.append(2.0*float(max(0.0, min([norm(loc = null_betas_loc, scale = null_betas_scale).sf(estimated_beta), 1.0 - norm(loc = null_betas_loc, scale = null_betas_scale).sf(estimated_beta)]))))

	logger.info('Calculated p. values for %s out of %s betas ...' % (len(beta_distributions), len(beta_distributions)))

	beta_pvals_adj = multipletests(pvals = beta_pvals, alpha = 0.05, method = 'fdr_bh')[1]
	gammas2betas['p'] = beta_pvals
	gammas2betas['padj'] = beta_pvals_adj

	new_p_cutoff = beta_pvals[pymin(range(len(beta_pvals_adj)), key=lambda i: pyabs(beta_pvals_adj[i] - float(padj_cutoffs[0])))]
	
	# Estimate statistical power
	beta_statistical_power = []
	if scale > 1:
		beta_corrected_effect_size = crispr_surf_statistical_power(sgRNA_indices = guideindices2bin.keys(), gammas2betas = gammas2betas, effect_size = effect_size, gamma_chosen = gamma_chosen, perturbation_profile = perturbation_profile, scale = scale)

	else:
		beta_corrected_effect_size = crispr_surf_statistical_power(sgRNA_indices = sgRNA_indices, gammas2betas = gammas2betas, effect_size = effect_size, gamma_chosen = gamma_chosen, perturbation_profile = perturbation_profile, scale = scale)

	for i in range(len(beta_corrected_effect_size)):

		# shifted_distribution = [x + beta_corrected_effect_size[i] for x in beta_distributions_null[i]]
		# percentile_cutoff = np.percentile(beta_distributions_null[i], (100.0 - float(new_p_cutoff)*100.0/2.0))

		beta_dist_null = np.array(beta_distributions_null[i])
		shifted_distribution = beta_dist_null + beta_corrected_effect_size[i]
		percentile_cutoff = np.percentile(beta_dist_null, (100.0 - float(new_p_cutoff)*100.0/2.0))

		if (i + 1)%500 == 0:
			logger.info('Calculated statistical power for %s out of %s betas ...' % ((i + 1), len(beta_distributions)))

		# beta_statistical_power.append(float(sum(x >= percentile_cutoff for x in shifted_distribution))/float(len(shifted_distribution)))

		beta_statistical_power.append((shifted_distribution > percentile_cutoff).sum() / float(len(shifted_distribution)))

	gammas2betas['power'] = beta_statistical_power

	return gammas2betas, replicate_parameters