### CRISPR-SURF Functions for Dash Web Application ###
pymin = min
pyabs = abs
import math
import numpy as np
import pandas as pd
from cvxpy import *
import logging
from scipy.stats import norm
from scipy.stats import laplace
import itertools
import scipy as sp
import random
import json
import time
from statsmodels.sandbox.stats.multicomp import multipletests

# Support Functions #
def str2bool(v):

	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		sys.exit(1)

nt_complement = dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})

def reverse_complement(seq):

	return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def total_count_normalization(gRNA_dict):

	samples_n = len(gRNA_dict[gRNA_dict.keys()[0]])
	guides_n = len(gRNA_dict)
	sample_sums = [0]*samples_n

	for gRNA in gRNA_dict:
		for index in range(samples_n):
			sample_sums[index] += gRNA_dict[gRNA][index]

	sample_average = float(sum(sample_sums))/float(samples_n)
	total_count_normalization_factor = [sample_average/x for x in sample_sums]

	return total_count_normalization_factor

def median_normalization(gRNA_dict):

	samples_n = len(gRNA_dict[list(gRNA_dict)[0]])
	guides_n = len(gRNA_dict)
	means = {k:math.exp((float(sum([math.log(v2 + 1.0) for v2 in v]))/samples_n)) for (k,v) in gRNA_dict.items() if sum(v) > 0}
	means_filter = {k:(lambda x: x if x > 0 else 1)(v) for (k,v) in means.items()}

	median_normalization_factor = [0]*samples_n
	for i in range(samples_n):
		mean_factor = [v[i]/means_filter[k] for (k,v) in gRNA_dict.items() if k in means_filter]
		corrected_factor = sorted(mean_factor)[len(mean_factor)//2]
		if corrected_factor > 0.0:
			median_normalization_factor[i] = 1.0/float(corrected_factor)

	return median_normalization_factor

def normalize_sgRNA_counts(gRNA_dict, method = 'total'):

	samples_n = len(gRNA_dict[list(gRNA_dict)[0]])

	if method == 'total':
		norm_factor = total_count_normalization(gRNA_dict)
	elif method == 'median':
		norm_factor = median_normalization(gRNA_dict)
	elif method == 'none':
		norm_factor = [1.0]*samples_n
	else:
		sys.exit(1)
    
	normalized_gRNA_dict = {k:[norm_factor[i]*v[i] for i in range(samples_n)] for (k,v) in gRNA_dict.items()}
	return normalized_gRNA_dict

# Deconvolution Functions #

def gaussian_pattern(characteristic_perturbation_range, scale, limit):

	"""
	Function to create gaussian function for screen-specific perturbation range.

	The characteristic_perturbation_range parameter refers to the average perturbation length for the screening strategy used.
	General guidelines: Cas9 ~ 10 bps, CRISPRa ~ 100 bps, CRISPRi ~ 200 bps

	The scale parameter scales the gaussian pattern depending on the resolution of your screen.
	Saturating mutagenesis screens using genetic perturbations (Cas9, Cpf1, base-editors) should have scale = 1.
	Enhancer discovery screens using epigenetic perturbations (CRISPRi, CRISPRa) should have scales of ~10 or ~20, respectively.

	The limit parameter determines how far the perturbation profile will reach.
	Any distances yielding a value below this limit (default: 0.05) will not be incorporated in the convolution.
	"""

	# Constant that incorporates characteristic perturbation length into gaussian pattern: gaussian(characteristic_perturbation_range) = 0.5
	constant = math.sqrt(-((float(characteristic_perturbation_range)**2.0)/(2.0*math.log(0.5))))

	# Build the gaussian perturbation profile
	gaussian_perturbation = [math.exp(-((float(distance)**2.0)/(2.0*float(constant)**2))) for distance in range(-limit, limit + 1)]

	# Rescale the perturbation profile
	gaussian_perturbation_rescaled = []
	if scale > 1:

		# New indices aggregating set of bps
		rescale_indices = np.arange(limit + 1, len(gaussian_perturbation), scale)
		for i in range(len(rescale_indices) - 1):
			gaussian_perturbation_rescaled.append(float(sum(gaussian_perturbation[rescale_indices[i]:rescale_indices[i + 1]]))/float(scale))

		# Mirror perturbation profile
		gaussian_perturbation_final = gaussian_perturbation_rescaled[::-1] + [1.0] + gaussian_perturbation_rescaled

	else:
		gaussian_perturbation_final = gaussian_perturbation

	return gaussian_perturbation_final

def crispr_surf_deconvolution(observations, chromosomes, sgRNA_indices, perturbation_profile, gamma_list, scale, uploads_dir, results_dir):

	"""
	Function to perform deconvolution for sgRNA LFC scores.
	The observations are a list of sgRNA LFC scores.
	The sgRNA_indices are a list of sgRNA perturbation indices.
	The perturbation_profile is the output of function gaussian_pattern.
	The gamma_list is a list of gammas to solve for betas.
	"""

	# time_per_gamma = {}
	# for g in gamma_list:
	# 	time_per_gamma[str(g)] = 0

	starttime = time.time()
	gammas2betas = {}
	# Loop through observations
	for rep_n in range(1, len(observations) + 1):

		print 'Replicate %s ' % rep_n
		# logger.info('Performing deconvolution on Replicate %s' % rep_n)

		gammas2betas[rep_n] = {}

		# Sort observations based on sgRNA indices
		sgRNA_indices, observations_rep, chromosomes = zip(*sorted(zip(sgRNA_indices, observations[rep_n - 1], chromosomes)))

		# Specify convolution shift and maximum distance for inference
		convolution_shift = len(perturbation_profile)/2
		maximum_distance = len(perturbation_profile)

		# Rescale observations and sgRNA indices
		guidescores2bin = {}
		guideindices2bin = {}

		# Make sure scale is integer
		scale = int(scale)

		if scale > 1:

			# Rescale genomic indices where sgRNAs are tiling
			bin_indices = np.arange(min(sgRNA_indices), max(sgRNA_indices), scale)
			rescaled_sgRNA_indices = np.array([(x + int(scale/2.0)) for x in bin_indices]) #  bin_indices[:-1]])
			rescaled_chromosomes = []

			# Assign sgRNAs to rescaled genomic indices
			current_rescaled_index = 0
			sgRNA_indices = sorted(sgRNA_indices)
			for sgRNA_index in range(len(sgRNA_indices)):

				for rescaled_index in range(current_rescaled_index, len(rescaled_sgRNA_indices)):

					if ((rescaled_sgRNA_indices[rescaled_index] - scale) <= sgRNA_indices[sgRNA_index] < (rescaled_sgRNA_indices[rescaled_index] + scale)):

						if rescaled_sgRNA_indices[rescaled_index] not in guidescores2bin:
							guidescores2bin[rescaled_sgRNA_indices[rescaled_index]] = []
							guideindices2bin[rescaled_sgRNA_indices[rescaled_index]] = []
							rescaled_chromosomes.append(chromosomes[sgRNA_index])

						guidescores2bin[rescaled_sgRNA_indices[rescaled_index]].append(observations_rep[sgRNA_index])
						guideindices2bin[rescaled_sgRNA_indices[rescaled_index]].append(sgRNA_indices[sgRNA_index])
						break

					else:
						current_rescaled_index += 1

			# Average across sgRNA scores within each rescale genomic index
			rescaled_sgRNA_indices_w_obs = sorted(guidescores2bin.keys())
			rescaled_observations = [np.mean(guidescores2bin[x]) for x in rescaled_sgRNA_indices_w_obs]

		else:
			rescaled_sgRNA_indices = np.array(sgRNA_indices)
			rescaled_sgRNA_indices_w_obs = np.array(sgRNA_indices)
			rescaled_observations = np.array(observations_rep)
			rescaled_chromosomes = np.array(chromosomes)

		# Identify groups within CRISPR tiling screen based on specified perturbation range
		group_boundaries = [0] + [(i + 1) for (i, j) in zip(range(len(rescaled_sgRNA_indices_w_obs) - 1), np.diff(rescaled_sgRNA_indices_w_obs)) if j > (maximum_distance*scale)] + [len(rescaled_sgRNA_indices_w_obs)]
		groups = []
		for i in range(1, len(group_boundaries)):
			groups += [i]*(group_boundaries[i] - group_boundaries[i - 1])

		# Set up regularized deconvolution optimization problem
		df = pd.DataFrame({'pos':rescaled_sgRNA_indices_w_obs, 'lfc':rescaled_observations, 'group':groups, 'chr': rescaled_chromosomes})

		# df.to_csv('rescaled_indices.csv')

		genomic_coordinates = []
		chromosomes_final = []
		delete_gammas = []

		# Iterate through groups and perform deconvolution
		for group in df.group.unique():

			# print 'Group %s' % group

			# Filtered dataframe to separate individual groups
			dff = df[df.group == group]

			# Make sure >1 sgRNA exists per group
			# if len(dff.index) > 1:

			# Assign relevant variables for optimization problem
			# y = dff.lfc.tolist()
			y = [0]*1 + dff.lfc.tolist() + [0]*1
			betas = Variable(len(np.arange(dff.pos.tolist()[0], dff.pos.tolist()[-1], scale).tolist()) + maximum_distance)

			# x_shift = [int(maximum_distance + (x - dff.pos.tolist()[0])/int(scale) - 1) for x in dff.pos.tolist()]
			x_shift_tmp = [int(maximum_distance + (x - dff.pos.tolist()[0])/int(scale) - 1) for x in dff.pos.tolist()]

			x_shift = list(range(x_shift_tmp[0] - 1, x_shift_tmp[0])) + x_shift_tmp + list(range(x_shift_tmp[-1] + 1, x_shift_tmp[-1] + 2))

			gamma = Parameter(sign = "positive")

			genomic_coordinates += np.arange(int(dff.pos.tolist()[0]), int(dff.pos.tolist()[-1]) + scale, scale).tolist()
			chromosomes_final += [dff['chr'].tolist()[0]]*len(np.arange(int(dff.pos.tolist()[0]), int(dff.pos.tolist()[-1]) + scale, scale).tolist())

			# Formulate optimization problem
			objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*sum_entries(abs(diff(betas))))
			p = Problem(objective)

			# Solve for varying lambdas
			for g in gamma_list:

				# timestart = time.time()

				g = float(g)

				# Make sure solver converges, otherwise delete gammas that fail
				try:

					if g not in gammas2betas[rep_n]:
						gammas2betas[rep_n][g] = []

					gamma.value = g
					result = p.solve()
					gammas2betas[rep_n][g] += np.array(betas.value).reshape(-1).tolist()[int(maximum_distance/2):-int(maximum_distance/2)]

				except:

					delete_gammas.append(g)
					continue

				# timestop = time.time()

				# time_per_gamma[str(g)] += timestop - timestart

		# Delete gammas that failed to converge
		for g in delete_gammas:
			print rep_n, g
			del gammas2betas[rep_n][g]

	endtime = time.time()

	time_per_gamma = (endtime - starttime)/float(len(gamma_list))/60.0

	print time_per_gamma

	gammas2betas['indices'] = genomic_coordinates
	gammas2betas['chr'] = chromosomes_final

	with open(uploads_dir + '/data.json', 'r') as f:
		json_string = f.readline().strip()
		data_dict = json.loads(json_string)

	data_dict['gammas2betas'] = gammas2betas
	data_dict['guideindices2bin'] = guideindices2bin
	data_dict['time_per_gamma'] = time_per_gamma

	with open(uploads_dir + '/data.json', 'w') as f:
		new_json_string = json.dumps(data_dict)
		f.write(new_json_string + '\n')

	return gammas2betas, guideindices2bin

def crispr_surf_deconvolution_simulations(negative_control_scores, sgRNA_indices, perturbation_profile, gamma_list, simulations_n, replicates, guideindices2bin, averaging_method, scale):

	"""
	Function to perform deconvolution for negative control sgRNA LFC scores.
	The negative controls are a list of negative control sgRNA LFC scores.
	The sgRNA_indices are a list of sgRNA perturbation indices.
	The perturbation_profile is the output of function gaussian_pattern.
	The gamma_list is a list of gammas to solve for betas.
	The simulations_n is the number of simulations to perform.
	"""

	# Specify convolution shift and maximum distance for inference
	convolution_shift = int((len(perturbation_profile) - 1)/2)
	maximum_distance = len(perturbation_profile)

	# Make sure scale is integer
	scale = int(scale)

	if scale > 1:
		tmp = [int(x) for x in guideindices2bin.keys()]
		rescaled_sgRNA_indices_w_obs = sorted(tmp)

	else:
		tmp = [int(x) for x in sgRNA_indices]
		rescaled_sgRNA_indices_w_obs = sorted(tmp)

	# Identify groups within CRISPR tiling screen based on specified perturbation range
	group_boundaries = [0] + [(i + 1) for (i, j) in zip(range(len(rescaled_sgRNA_indices_w_obs) - 1), np.diff(rescaled_sgRNA_indices_w_obs)) if j > (maximum_distance*scale)] + [len(rescaled_sgRNA_indices_w_obs)]
	groups = []
	for i in range(1, len(group_boundaries)):
		groups += [i]*(group_boundaries[i] - group_boundaries[i - 1])

	# Iterate through n simulations
	beta_distributions = {}
	for n in range(1, simulations_n + 1):

		if n%100 == 0:
			print 'Simulation %s out of %s ...' % (str(n), str(simulations_n))
			# logger.info('Simulation %s out of %s ...' % (str(n), str(simulations_n)))

		replicate_store = {}
		for r in range(replicates):

			if negative_control_scores[0] == 'gaussian':
				if scale > 1:
					rescaled_observations = []
					for scaled_index in rescaled_sgRNA_indices_w_obs:
						rescaled_observations.append(np.mean(norm.rvs(loc = negative_control_scores[1][r][0], scale = negative_control_scores[1][r][1], size = len(guideindices2bin[str(scaled_index)]))))

				else:
					rescaled_observations = norm.rvs(loc = negative_control_scores[1][r][0], scale = negative_control_scores[1][r][1], size = len(rescaled_sgRNA_indices_w_obs))

			elif negative_control_scores[0] == 'laplace':
				if scale > 1:
					rescaled_observations = []
					for scaled_index in rescaled_sgRNA_indices_w_obs:
						rescaled_observations.append(np.mean(laplace.rvs(loc = negative_control_scores[1][r][0], scale = negative_control_scores[1][r][1], size = len(guideindices2bin[str(scaled_index)]))))

				else:
					rescaled_observations = laplace.rvs(loc = negative_control_scores[1][r][0], scale = negative_control_scores[1][r][1], size = len(rescaled_sgRNA_indices_w_obs))

			elif negative_control_scores[0] == 'negative_control_guides':
				if scale > 1:
					rescaled_observations = []
					for scaled_index in rescaled_sgRNA_indices_w_obs:
						rescaled_observations.append(np.mean(np.random.choice(negative_control_scores[1][r], len(guideindices2bin[str(scaled_index)]), replace = True)))

				else:
					rescaled_observations = np.random.choice(negative_control_scores[1][r], len(rescaled_sgRNA_indices_w_obs), replace = True)

			# Set up regularized deconvolution optimization problem
			df = pd.DataFrame({'pos':rescaled_sgRNA_indices_w_obs, 'lfc':rescaled_observations, 'group':groups})

			genomic_coordinates = []
			gammas2betas = {}
			delete_gammas = []

			# Iterate through groups and perform deconvolution
			for group in df.group.unique():

				# Filtered dataframe to separate individual groups
				dff = df[df.group == group]

				# Make sure >1 sgRNA exists per group
				# if len(dff.index) > 1:

				# Assign relevant variables for optimization problem
				y = dff.lfc.tolist()
				betas = Variable(len(np.arange(dff.pos.tolist()[0], dff.pos.tolist()[-1], scale).tolist()) + maximum_distance)

				x_shift = [int(maximum_distance + (x - dff.pos.tolist()[0])/int(scale)) for x in dff.pos.tolist()]

				gamma = Parameter(sign = "positive")

				genomic_coordinates += np.arange(int(dff.pos.tolist()[0]), int(dff.pos.tolist()[-1]) + scale, scale).tolist()

				# Formulate optimization problem
				objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*sum_entries(abs(diff(betas))))
				p = Problem(objective)

				# Solve for varying lambdas
				for g in gamma_list:

					g = float(g)

					# Make sure solver converges, otherwise delete gammas that fail
					try:

						if g not in gammas2betas:
							gammas2betas[g] = []

						gamma.value = g
						result = p.solve()
						gammas2betas[g] += np.array(betas.value).reshape(-1).tolist()[int(maximum_distance/2):-int(maximum_distance/2)]

					except:

						delete_gammas.append(g)
						continue

			# Delete gammas that failed to converge
			for g in delete_gammas:
				del gammas2betas[g]

			gammas2betas['indices'] = genomic_coordinates

			# Add to replicate store
			replicate_store[r] = gammas2betas[float(gamma_list[0])]

		# Create combined deconvolved signals from replicates for simulation
		deconvolved_signal = {}
		for i in replicate_store.keys():

			for j in range(len(replicate_store[i])):

				if j not in deconvolved_signal:
					deconvolved_signal[j] = []

				deconvolved_signal[j].append(replicate_store[i][j])

		# Create mean or median profile
		if averaging_method == 'mean':
			combine_simulations = [np.mean(deconvolved_signal[x]) for x in deconvolved_signal]

		elif averaging_method == 'median':
			combine_simulations = [np.median(deconvolved_signal[x]) for x in deconvolved_signal]

		for i in range(len(combine_simulations)):
			try:
				beta_distributions[i].append(combine_simulations[i])
			except:
				beta_distributions[i] = [combine_simulations[i]]

	return beta_distributions

def crispr_surf_statistical_power(sgRNA_indices, gammas2betas, effect_size, gamma_chosen, perturbation_profile, scale):
	"""
	Function to estimate statistical power for beta profile.
	"""

	# Create perturbation profiles
	maximum_distance = len(perturbation_profile)

	# Loop through beta indices
	beta_indices = [int(x) for x in gammas2betas['indices']]
	sgRNA_indices = [int(x) for x in sgRNA_indices]
	beta_statistical_power = []
	total_betas = len(beta_indices)
	for beta_index in range(len(beta_indices)):

		if (beta_index + 1)%100 == 0:
			print 'Estimating statistical power for %s out of %s units ...' % (str(beta_index + 1), str(total_betas))
			# logger.info('Estimating statistical power for %s out of %s units ...' % (str(beta_index + 1), str(total_betas)))

		# Construct underlying truth
		bp_range = np.arange(beta_indices[beta_index] - (maximum_distance)*scale, beta_indices[beta_index] + (maximum_distance + 1)*scale, scale)
		landscape = [0]*len(bp_range)
		landscape[maximum_distance] += effect_size

		# Convolve with CRISPR perturbation
		crispr_convolved_landscape = np.convolve(perturbation_profile, landscape, mode = 'same')
		sgRNA_indices_truth = [x[0] for x in zip(bp_range, crispr_convolved_landscape) if int(x[0]) in sgRNA_indices]
		signal_truth = [x[1] for x in zip(bp_range, crispr_convolved_landscape) if int(x[0]) in sgRNA_indices]

		rescaled_sgRNA_indices_w_obs = np.array(sgRNA_indices_truth)
		rescaled_observations = np.array(signal_truth)

		# Set up regularized deconvolution optimization problem
		df = pd.DataFrame({'pos':rescaled_sgRNA_indices_w_obs, 'lfc':rescaled_observations})

		# Assign relevant variables for optimization problem
		y = df.lfc.tolist()
		betas = Variable(len(np.arange(df.pos.tolist()[0], df.pos.tolist()[-1], scale).tolist()) + maximum_distance)

		x_shift = [int(maximum_distance + (x - df.pos.tolist()[0])/int(scale) - 1) for x in df.pos.tolist()]

		gamma = Parameter(sign = "positive")

		if maximum_distance%2 == 1:
			genomic_coordinates = np.arange(int(df.pos.tolist()[0]) - int(maximum_distance/2)*scale, int(df.pos.tolist()[-1]) + int(maximum_distance/2)*scale + scale, scale).tolist()
		else:
			genomic_coordinates = np.arange(int(dff.pos.tolist()[0]) - int(maximum_distance/2)*scale, int(dff.pos.tolist()[-1]) + int(maximum_distance/2)*scale, scale).tolist()

		# Formulate optimization problem
		objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*sum_entries(abs(diff(betas))))
		p = Problem(objective)

		gamma.value = float(gamma_chosen)
		result = p.solve()
		betas_inferred = np.array(betas.value).reshape(-1).tolist()

		try:
			index_truth = genomic_coordinates.index(beta_indices[beta_index])
			beta_statistical_power.append(betas_inferred[index_truth])
		except:
			index_truth = min(range(len(genomic_coordinates)), key=lambda i: abs(genomic_coordinates[i] - beta_indices[beta_index]))
			beta_statistical_power.append(betas_inferred[index_truth])

	print 'Estimating statistical power for %s out of %s units ...' % (str(total_betas), str(total_betas))
	# logger.info('Estimating statistical power for %s out of %s units ...' % (str(total_betas), str(total_betas)))

	return beta_statistical_power

# Gamma Functions #
def crispr_surf_find_gamma(gammas2betas, correlation_ratio_start, correlation_ratio_stop, correlation_ratio_opt, uploads_dir, results_dir):

	"""
	Function to find optimal gamma range to be used for regularization based on empirical simulations.
	Input is gammas2betas dictionary which is the output of function crispr_surf_deconvolution. Dictionary structure: Keys: 1) Replicate, 2) Gamma
	The characteristic_perturbation_length and noise parameters inform the Pearson correlation ratio between R2_max and R2_gamma_opt.
	"""
	# print gammas2betas
	# Enumerate all pairs of biological replicates and calculate pearson correlation across gamma range
	replicate_pair_correlations = {}
	for pair in itertools.combinations([x for x in gammas2betas.keys() if ((x != 'combined') and (x != 'gamma_chosen') and (x != 'padj') and (x != 'indices') and (x != 'chr'))], 2):

		replicate1, replicate2 = pair[0], pair[1]
		replicate_pair_id = '-'.join(map(str, [replicate1, replicate2]))
		gamma_list = sorted([float(x) for x in gammas2betas[replicate1].keys() if ((x != 'combined') and (x != 'gamma_chosen') and (x != 'padj') and (x != 'indices') and (x != 'chr'))])

		replicate_pair_correlations[replicate_pair_id] = []
		for gamma in gamma_list:

			replicate_pair_correlations[replicate_pair_id].append(sp.stats.pearsonr(gammas2betas[replicate1][str(gamma)], gammas2betas[replicate2][str(gamma)])[0])

	# Rescale all R2 curves and aggregate to find optimal gamma range
	correlation_curve = {}
	for pair in replicate_pair_correlations:

		for i in range(len(replicate_pair_correlations[pair])):

			if i not in correlation_curve:
				correlation_curve[i] = []

			correlation_curve[i].append(replicate_pair_correlations[pair][i])

	# print correlation_curve

	# Normalize each correlation curve and then average across curves
	correlation_curve_averaged = [float(sum(correlation_curve[x]))/float(len(correlation_curve[x])) for x in correlation_curve]
	correlation_curve_rescaled = [float(x)/float(max(correlation_curve_averaged)) for x in correlation_curve_averaged]

	# print correlation_curve_rescaled

	max_index = correlation_curve_rescaled.index(max(correlation_curve_rescaled))
	correlation_curve_rescaled_capped = correlation_curve_rescaled[:max_index + 1]

	# print correlation_curve_rescaled_capped

	# Return optimal gamma range
	gamma_index_start = min(range(len(correlation_curve_rescaled_capped)), key=lambda i: (correlation_curve_rescaled_capped[i] - float(correlation_ratio_start)).__abs__())
	gamma_index_stop = min(range(len(correlation_curve_rescaled_capped)), key=lambda i: (correlation_curve_rescaled_capped[i] - float(correlation_ratio_stop)).__abs__())
	gamma_index_opt = min(range(len(correlation_curve_rescaled_capped)), key=lambda i: (correlation_curve_rescaled_capped[i] - float(correlation_ratio_opt)).__abs__())

	gamma_range = [x for x in zip(gamma_list[gamma_index_start:], correlation_curve_rescaled_capped[gamma_index_start:])]

	gamma_start = gamma_list[gamma_index_start]
	gamma_stop = gamma_list[gamma_index_stop]
	gamma_opt = gamma_list[gamma_index_opt]

	with open(results_dir + '/correlation_curve_gamma.csv', 'w') as f:
		f.write(','.join(map(str, ['Gamma', 'Corr_Avg', 'Corr_Avg_Scaled'])) + '\n')

		for i in zip(gamma_list, correlation_curve_averaged, correlation_curve_rescaled):
			f.write(','.join(map(str, i)) + '\n')

		f.write('\n')
		f.write('Correlation Ratio:,%s' % correlation_ratio_opt + '\n')
		f.write('Gamma Chosen:,%s' % gamma_opt + '\n')

	with open(uploads_dir + '/data.json', 'r') as f:
		json_string = f.readline().strip()
		data_dict = json.loads(json_string)

	data_dict['gamma_range'] = gamma_range

	with open(uploads_dir + '/data.json', 'w') as f:
		new_json_string = json.dumps(data_dict)
		f.write(new_json_string + '\n')

	with open(results_dir + '/flag1.txt', 'w') as f:
		f.write('Completed deconvolution and replicate correlation steps!')

	return gamma_range, gamma_opt

# Statistical Significance Functions #
def crispr_surf_deconvolved_signal(gammas2betas, gamma_chosen, averaging_method, uploads_dir, results_dir):

	"""
	Function to construct final deconvolved signal based on chosen gamma.
	The gammas2betas dictionary has structure: Keys: 1) Replicate, 2) Gamma
	Averaging method can be mean or median to combine all biological replicates.
	"""

	# Create combined deconvolved signals from replicates
	gamma_chosen = str(gamma_chosen)
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

	# hack
	json_good = False
	while not json_good:
	    with open(uploads_dir + '/data2.json', 'r') as f:
	        json_string = f.readline().strip()
	        try:
	            data_dict2 = json.loads(json_string)
	            json_good = True
	        except:
	            pass

	data_dict2['gammas2betas'] = gammas2betas

	with open(uploads_dir + '/data2.json', 'w') as f:
		new_json_string = json.dumps(data_dict2)
		f.write(new_json_string + '\n')

	return gammas2betas

def crispr_surf_statistical_significance(sgRNA_summary_table, sgRNA_indices, perturbation_profile, gammas2betas, simulation_type, simulation_n, guideindices2bin, averaging_method, padj_cutoffs, effect_size, limit, scale, rapid_mode, uploads_dir, results_dir):

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
	print 'Performing %s simulations to construct beta null distributions ...' % (simulation_n)
	# logger.info('Performing %s simulations to construct beta null distributions ...' % (simulation_n))

	if simulation_type == 'negative_control':
		if 'negative_control' not in df_summary_table['sgRNA_Type'].unique().tolist():
			simulation_type = 'gaussian'

	replicate_parameters = []
	if simulation_type == 'negative_control':

		replicate_parameters.append('NA')

		# Grab all negative control sgRNA lfc scores
		negative_control_guide_scores = []
		for i in range(1, int(replicates) + 1):
			negative_control_guide_scores.append(np.array(df_summary_table.loc[(df_summary_table['sgRNA_Type'] == 'negative_control'), ['Log2FC_Replicate' + str(i)]]).flatten().tolist())

		# Construct many simulated null arrays to perform deconvolution
		beta_distributions_null = crispr_surf_deconvolution_simulations(negative_control_scores = ['negative_control_guides', negative_control_guide_scores], sgRNA_indices = sgRNA_indices, perturbation_profile = perturbation_profile, gamma_list = [gamma_chosen], simulations_n = simulation_n, replicates = replicates, guideindices2bin = guideindices2bin, averaging_method = averaging_method, scale = scale)

	elif simulation_type == 'laplace':

		# Parameterize observed signal with laplace distribution (assume majority of observation sgRNAs are null)
		observation_median = np.median()
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

	elif simulation_type == 'gaussian':

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
	print 'Calculating p. values for %s betas ...' % (len(beta_distributions))
	# logger.info('Calculating p. values for %s betas ...' % (len(beta_distributions)))

	beta_pvals = []
	if not rapid_mode:
		for i in range(len(beta_distributions)):

			if (i + 1)%100 == 0:
				print 'Calculated p. values for %s out of %s betas ...' % ((i + 1), len(beta_distributions))
				# logger.info('Calculated p. values for %s out of %s betas ...' % ((i + 1), len(beta_distributions)))

			estimated_beta = beta_distributions[i]
			null_betas = beta_distributions_null[i]

			beta_pvals.append(2.0*float(max(0.0, min(sum(x >= estimated_beta for x in null_betas), sum(x <= estimated_beta for x in null_betas))))/float(len(null_betas)))
	else:
		null_betas = []

		for i in range(len(beta_distributions_null)):
			null_betas += beta_distributions_null[i]

		null_betas = sorted(null_betas)
		print 'Aggregated beta null distribution size: %s ...' % (len(null_betas))
		# logger.info('Aggregated beta null distribution size: %s ...' % (len(null_betas)))

		for i in range(len(beta_distributions)):

			if (i + 1)%100 == 0:
				print 'Calculated p. values for %s out of %s betas ...' % ((i + 1), len(beta_distributions))
				# logger.info('Calculated p. values for %s out of %s betas ...' % ((i + 1), len(beta_distributions)))

			estimated_beta = beta_distributions[i]

			# Decide which tail to calculate
			if estimated_beta < np.median(null_betas):

				for j in range(len(null_betas)):
					if null_betas[j] > estimated_beta:
						pval = 2.0*max(0.0, float(j))/float(len(null_betas))
						beta_pvals.append(pval)
						break
			else:

				for j in range(len(null_betas))[::-1]:
					if null_betas[j] < estimated_beta:
						pval = 2.0*max(0.0, float(len(null_betas) - j - 1))/float(len(null_betas))
						beta_pvals.append(pval)
						break

	print 'Calculated p. values for %s out of %s betas ...' % (len(beta_distributions), len(beta_distributions))
	# logger.info('Calculated p. values for %s out of %s betas ...' % (len(beta_distributions), len(beta_distributions)))

	beta_pvals_adj = multipletests(pvals = beta_pvals, alpha = 0.05, method = 'fdr_bh')[1]
	gammas2betas['p'] = list(beta_pvals)
	gammas2betas['padj'] = list(beta_pvals_adj)

	new_p_cutoff = beta_pvals[pymin(range(len(beta_pvals_adj)), key=lambda i: pyabs(beta_pvals_adj[i] - float(padj_cutoffs[0])))]

	# Estimate statistical power
	beta_statistical_power = []
	if scale > 1:
		beta_corrected_effect_size = crispr_surf_statistical_power(sgRNA_indices = guideindices2bin.keys(), gammas2betas = gammas2betas, effect_size = effect_size, gamma_chosen = gamma_chosen, perturbation_profile = perturbation_profile, scale = scale)

	else:
		beta_corrected_effect_size = crispr_surf_statistical_power(sgRNA_indices = sgRNA_indices, gammas2betas = gammas2betas, effect_size = effect_size, gamma_chosen = gamma_chosen, perturbation_profile = perturbation_profile, scale = scale)

	if not rapid_mode:

		for i in range(len(beta_corrected_effect_size)):

			shifted_distribution = [x + beta_corrected_effect_size[i] for x in beta_distributions_null[i]]
			percentile_cutoff = np.percentile(beta_distributions_null[i], (100.0 - float(new_p_cutoff)*100.0/2.0))
			if (i + 1)%100 == 0:
				print 'Calculated statistical power for %s out of %s betas ...' % ((i + 1), len(beta_distributions))
				# logger.info('Calculated statistical power for %s out of %s betas ...' % ((i + 1), len(beta_distributions)))

			beta_statistical_power.append(float(sum(x >= percentile_cutoff for x in shifted_distribution))/float(len(shifted_distribution)))

	else:

		percentile_cutoff = np.percentile(null_betas, (100.0 - float(new_p_cutoff)*100.0/2.0))
		for i in range(len(beta_corrected_effect_size)):

			for j in range(len(null_betas)):
				if (null_betas[j] + beta_corrected_effect_size[i]) > percentile_cutoff:
					beta_statistical_power.append(float(len(null_betas) - j)/float(len(null_betas)))
					break

			if (i + 1)%100 == 0:
				print 'Calculated statistical power for %s out of %s betas ...' % ((i + 1), len(beta_distributions))
				# logger.info('Calculated statistical power for %s out of %s betas ...' % ((i + 1), len(beta_distributions)))

	gammas2betas['power'] = list(beta_statistical_power)

	# hack
	json_good = False
	while not json_good:
	    with open(uploads_dir + '/data2.json', 'r') as f:
	        json_string = f.readline().strip()
	        try:
	            data_dict2 = json.loads(json_string)
	            json_good = True
	        except:
	            pass

	data_dict2['gammas2betas'] = gammas2betas

	with open(uploads_dir + '/data2.json', 'w') as f:
		new_json_string = json.dumps(data_dict2)
		f.write(new_json_string + '\n')

	with open(results_dir + '/flag2.txt', 'w') as f:
		f.write('Completed simulation and statistical significance steps!')

	return gammas2betas, replicate_parameters

# Outputs Functions #
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
				closest_index = gammas2betas['indices'].index(int(guideindices2bin_inverted[int(row['Perturbation_Index'])]))

			else:
				closest_index = gammas2betas['indices'].index(int(row['Perturbation_Index']))

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
	with open(out_dir + '/' + 'sgRNAs_summary_table_updated.csv', 'w') as f:

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
