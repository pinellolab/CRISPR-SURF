#!/usr/bin/env python
import math
import numpy as np
import pandas as pd
from cvxpy import *
import logging
from scipy.stats import norm
from scipy.stats import laplace

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

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

def crispr_surf_deconvolution(observations, chromosomes, sgRNA_indices, perturbation_profile, gamma_list, scale):

	"""
	Function to perform deconvolution for sgRNA LFC scores.
	The observations are a list of sgRNA LFC scores.
	The sgRNA_indices are a list of sgRNA perturbation indices.
	The perturbation_profile is the output of function gaussian_pattern.
	The gamma_list is a list of gammas to solve for betas.
	"""

	# Sort observations based on sgRNA indices
	sgRNA_indices, observations, chromosomes = zip(*sorted(zip(sgRNA_indices, observations, chromosomes)))

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

					guidescores2bin[rescaled_sgRNA_indices[rescaled_index]].append(observations[sgRNA_index])
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
		rescaled_observations = np.array(observations)
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
	gammas2betas = {}
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
		# y = np.array(y).reshape(len(y), 1)
		betas = Variable(len(np.arange(dff.pos.tolist()[0], dff.pos.tolist()[-1], scale).tolist()) + maximum_distance)

		# x_shift = [int(maximum_distance + (x - dff.pos.tolist()[0])/int(scale) - 1) for x in dff.pos.tolist()]
		x_shift_tmp = [int(maximum_distance + (x - dff.pos.tolist()[0])/int(scale) - 1) for x in dff.pos.tolist()]

		x_shift = list(range(x_shift_tmp[0] - 1, x_shift_tmp[0])) + x_shift_tmp + list(range(x_shift_tmp[-1] + 1, x_shift_tmp[-1] + 2))

		gamma = Parameter(sign = "positive")
		# gamma = Parameter(nonneg = True)

		genomic_coordinates += np.arange(int(dff.pos.tolist()[0]), int(dff.pos.tolist()[-1]) + scale, scale).tolist()
		chromosomes_final += [dff['chr'].tolist()[0]]*len(np.arange(int(dff.pos.tolist()[0]), int(dff.pos.tolist()[-1]) + scale, scale).tolist())

		# Formulate optimization problem
		objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*sum_entries(abs(diff(betas))))
		# objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*tv(betas))
		p = Problem(objective)

		# Solve for varying lambdas
		for g in gamma_list:

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
	gammas2betas['chr'] = chromosomes_final

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
		rescaled_sgRNA_indices_w_obs = sorted(guideindices2bin.keys())

	else:
		rescaled_sgRNA_indices_w_obs = np.array(sorted(sgRNA_indices))

	# Identify groups within CRISPR tiling screen based on specified perturbation range
	group_boundaries = [0] + [(i + 1) for (i, j) in zip(range(len(rescaled_sgRNA_indices_w_obs) - 1), np.diff(rescaled_sgRNA_indices_w_obs)) if j > (maximum_distance*scale)] + [len(rescaled_sgRNA_indices_w_obs)]
	groups = []
	for i in range(1, len(group_boundaries)):
		groups += [i]*(group_boundaries[i] - group_boundaries[i - 1])

	# Iterate through n simulations
	beta_distributions = {}
	for n in range(1, simulations_n + 1):

		if n%100 == 0:
			logger.info('Simulation %s out of %s ...' % (str(n), str(simulations_n)))

		replicate_store = {}
		for r in range(replicates):

			if negative_control_scores[0] == 'gaussian':
				if scale > 1:
					rescaled_observations = []
					for scaled_index in rescaled_sgRNA_indices_w_obs:
						rescaled_observations.append(np.mean(norm.rvs(loc = negative_control_scores[1][r][0], scale = negative_control_scores[1][r][1], size = len(guideindices2bin[scaled_index]))))

				else:
					rescaled_observations = norm.rvs(loc = negative_control_scores[1][r][0], scale = negative_control_scores[1][r][1], size = len(rescaled_sgRNA_indices_w_obs))

			elif negative_control_scores[0] == 'laplace':
				if scale > 1:
					rescaled_observations = []
					for scaled_index in rescaled_sgRNA_indices_w_obs:
						rescaled_observations.append(np.mean(laplace.rvs(loc = negative_control_scores[1][r][0], scale = negative_control_scores[1][r][1], size = len(guideindices2bin[scaled_index]))))

				else:
					rescaled_observations = laplace.rvs(loc = negative_control_scores[1][r][0], scale = negative_control_scores[1][r][1], size = len(rescaled_sgRNA_indices_w_obs))

			elif negative_control_scores[0] == 'negative_control_guides':
				if scale > 1:
					rescaled_observations = []
					for scaled_index in rescaled_sgRNA_indices_w_obs:
						rescaled_observations.append(np.mean(np.random.choice(negative_control_scores[1][r], len(guideindices2bin[scaled_index]), replace = True)))

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
				# y = np.array(y).reshape(len(y), 1)
				betas = Variable(len(np.arange(dff.pos.tolist()[0], dff.pos.tolist()[-1], scale).tolist()) + maximum_distance)

				x_shift = [int(maximum_distance + (x - dff.pos.tolist()[0])/int(scale)) for x in dff.pos.tolist()]

				gamma = Parameter(sign = "positive")
				# gamma = Parameter(nonneg = True)

				genomic_coordinates += np.arange(int(dff.pos.tolist()[0]), int(dff.pos.tolist()[-1]) + scale, scale).tolist()

				# Formulate optimization problem
				objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*sum_entries(abs(diff(betas))))
				# objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*tv(betas))
				p = Problem(objective)

				# Solve for varying lambdas
				for g in gamma_list:

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
			replicate_store[r] = gammas2betas[gamma_list[0]]

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
	beta_indices = gammas2betas['indices']
	beta_statistical_power = []
	total_betas = len(beta_indices)
	for beta_index in range(len(beta_indices)):

		if (beta_index + 1)%100 == 0:
			logger.info('Estimating statistical power for %s out of %s units ...' % (str(beta_index + 1), str(total_betas)))

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
		# y = np.array(y).reshape(len(y), 1)
		betas = Variable(len(np.arange(df.pos.tolist()[0], df.pos.tolist()[-1], scale).tolist()) + maximum_distance)

		x_shift = [int(maximum_distance + (x - df.pos.tolist()[0])/int(scale) - 1) for x in df.pos.tolist()]

		gamma = Parameter(sign = "positive")
		# gamma = Parameter(nonneg = True)

		if maximum_distance%2 == 1:
			genomic_coordinates = np.arange(int(df.pos.tolist()[0]) - int(maximum_distance/2)*scale, int(df.pos.tolist()[-1]) + int(maximum_distance/2)*scale + scale, scale).tolist()
		else:
			genomic_coordinates = np.arange(int(dff.pos.tolist()[0]) - int(maximum_distance/2)*scale, int(dff.pos.tolist()[-1]) + int(maximum_distance/2)*scale, scale).tolist()

		# Formulate optimization problem
		objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*sum_entries(abs(diff(betas))))
		# objective = Minimize(0.5*sum_squares(y - conv(perturbation_profile, betas)[x_shift]) + gamma*tv(betas))
		p = Problem(objective)

		gamma.value = gamma_chosen
		result = p.solve()
		betas_inferred = np.array(betas.value).reshape(-1).tolist()

		try:
			index_truth = genomic_coordinates.index(beta_indices[beta_index])
			beta_statistical_power.append(betas_inferred[index_truth])
		except:
			index_truth = min(range(len(genomic_coordinates)), key=lambda i: abs(genomic_coordinates[i] - beta_indices[beta_index]))
			beta_statistical_power.append(betas_inferred[index_truth])

	logger.info('Estimating statistical power for %s out of %s units ...' % (str(total_betas), str(total_betas)))

	return beta_statistical_power