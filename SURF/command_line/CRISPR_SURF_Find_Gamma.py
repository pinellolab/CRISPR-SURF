#!/usr/bin/env python
import itertools
import pandas as pd
import scipy as sp

def crispr_surf_find_gamma(gammas2betas, correlation_ratio_start, correlation_ratio_stop, correlation_ratio_opt, out_dir):

	"""
	Function to find optimal gamma range to be used for regularization based on empirical simulations.
	Input is gammas2betas dictionary which is the output of function crispr_surf_deconvolution. Dictionary structure: Keys: 1) Replicate, 2) Gamma
	The characteristic_perturbation_length and noise parameters inform the Pearson correlation ratio between R2_max and R2_gamma_opt.
	"""

	# Enumerate all pairs of biological replicates and calculate pearson correlation across gamma range
	replicate_pair_correlations = {}
	for pair in itertools.combinations([x for x in gammas2betas.keys() if ((x != 'combined') and (x != 'gamma_chosen') and (x != 'padj') and (x != 'indices') and (x != 'chr'))], 2):

		replicate1, replicate2 = pair[0], pair[1]
		replicate_pair_id = '-'.join(map(str, [replicate1, replicate2]))
		gamma_list = sorted([x for x in gammas2betas[replicate1].keys() if ((x != 'combined') and (x != 'gamma_chosen') and (x != 'padj') and (x != 'indices') and (x != 'chr'))])

		replicate_pair_correlations[replicate_pair_id] = []

		for gamma in gamma_list:

			replicate_pair_correlations[replicate_pair_id].append(sp.stats.pearsonr(gammas2betas[replicate1][gamma], gammas2betas[replicate2][gamma])[0])

	# Rescale all R2 curves and aggregate to find optimal gamma range
	correlation_curve = {}
	for pair in replicate_pair_correlations:

		for i in range(len(replicate_pair_correlations[pair])):

			if i not in correlation_curve:
				correlation_curve[i] = []

			correlation_curve[i].append(replicate_pair_correlations[pair][i])

	# Normalize each correlation curve and then average across curves
	correlation_curve_averaged = [float(sum(correlation_curve[x]))/float(len(correlation_curve[x])) for x in correlation_curve]
	correlation_curve_rescaled = [float(x)/float(max(correlation_curve_averaged)) for x in correlation_curve_averaged]

	max_index = correlation_curve_rescaled.index(max(correlation_curve_rescaled))
	correlation_curve_rescaled_capped = correlation_curve_rescaled[:max_index + 1]

	# Return optimal gamma range
	gamma_index_start = min(range(len(correlation_curve_rescaled_capped)), key=lambda i: abs(correlation_curve_rescaled_capped[i] - float(correlation_ratio_start)))
	gamma_index_stop = min(range(len(correlation_curve_rescaled_capped)), key=lambda i: abs(correlation_curve_rescaled_capped[i] - float(correlation_ratio_stop)))
	gamma_index_opt = min(range(len(correlation_curve_rescaled_capped)), key=lambda i: abs(correlation_curve_rescaled_capped[i] - float(correlation_ratio_opt)))

	gamma_range = [x for x in zip(gamma_list[gamma_index_start:], correlation_curve_rescaled_capped[gamma_index_start:])]

	gamma_start = gamma_list[gamma_index_start]
	gamma_stop = gamma_list[gamma_index_stop]
	gamma_opt = gamma_list[gamma_index_opt]

	with open(out_dir + '/correlation_curve_gamma.csv', 'w') as f:
		f.write(','.join(map(str, ['Gamma', 'Corr_Avg', 'Corr_Avg_Scaled'])) + '\n')

		for i in zip(gamma_list, correlation_curve_averaged, correlation_curve_rescaled):
			f.write(','.join(map(str, i)) + '\n')

		f.write('\n')
		f.write('Correlation Ratio:,%s' % correlation_ratio_opt + '\n')
		f.write('Gamma Chosen:,%s' % gamma_opt + '\n')

	return gamma_range, gamma_opt