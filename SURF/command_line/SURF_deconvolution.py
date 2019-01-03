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
from CRISPR_SURF_Deconvolution import gaussian_pattern, crispr_surf_deconvolution
from CRISPR_SURF_Find_Lambda import crispr_surf_find_lambda
from CRISPR_SURF_Statistical_Significance import crispr_surf_deconvolved_signal, crispr_surf_statistical_significance
from CRISPR_SURF_Outputs import crispr_surf_sgRNA_summary_table_update, complete_beta_profile, crispr_surf_significant_regions, crispr_surf_IGV

# Time stamp analysis start
timestr = time.strftime("%Y-%m-%d_%H.%M.%S")

##### Argument handeling
parser = argparse.ArgumentParser(description = 'CRISPR-SURF command line tool to analyze CRISPR tiling screen data across various screening modalities.')
parser.add_argument('-f', '--sgRNAs_summary_table', type = str, required = True, help = 'Input sgRNA summary file with the following required columns: Chr, Start, Stop, Perturbation_Index, sgRNA_Sequence, Strand, sgRNA_Type, Log2FC_Replicate1. Use crisprsurf_counts_CL.py to generate this sgRNA summary file from FASTQs.')
parser.add_argument('-pert', '--perturbation_type', type = str, choices = ['be', 'cas9', 'cpf1', 'crispri', 'crispra', 'cas9_avg'], required = True, help = 'Perturbation type used (cas9, cpf1, crispri, crispra, be, cas9_avg).')
parser.add_argument('-range', '--characteristic_perturbation_range', type = int, default = 0, help = 'Characteristic perturbation length. If 0 (default), the --perturbation_type argument will be used to set an appropriate perturbation range.')
parser.add_argument('-scale', '--scale', type = int, default = 0, help = 'Scaling factor to efficiently perform deconvolution with negligible consequences. If 0 (default), the --characteristic_perturbation_range argument will be used to set an appropriate scaling factor.')
parser.add_argument('-limit', '--limit', type = int, default = 0, help = 'Maximum distance between two sgRNAs to perform inference on bp in-between. Sets the boundaries of the gaussian profile to perform efficient deconvolution. If 0 (default), the --perturbation_type argument will be used to set an appropriate limit.')
parser.add_argument('-avg', '--averaging_method', type = str, default = 'median', choices = ['mean', 'median'], help = 'The averaging method to be performed to combine biological replicates (mean, median).')
parser.add_argument('-null_dist', '--null_distribution', type = str, default = 'gaussian', choices = ['negative_control', 'gaussian', 'laplace'], help = 'The method of building a null distribution for each smoothed beta score (negative_control, gaussian, laplace).')
parser.add_argument('-sim_n', '--simulation_n', type = int, default = 1000, help = 'The number of simulations to perform for construction of the null distribution.')
parser.add_argument('-test_type', '--test_type', type = str, choices = ['parametric', 'nonparametric'], default = 'parametric', help = 'Parametric or non-parametric test for betas.')
parser.add_argument('-lambda_list', '--lambda_list', default = 0, nargs = '+', help = 'List of lambdas (regularization parameter) to use during deconvolution step. If 0 (default), the --perturbation_type argument will be used to set a reasonable lambda list. Example: 1,2,3,4,5,6,7,8,9,10')
parser.add_argument('-lambda_val', '--lambda_val', type = float, default = 0, help = 'The lambda to use to use during deconvolution step. If 0 (default), the --lambda_list argument will be used.')
parser.add_argument('-corr', '--correlation', type = float, default = 0, help = 'The correlation between biological replicates to determine a reasonable lambda for the deconvolution operation. if 0 (default), the --characteristic_perturbation_range argument will be used to set an appropriate correlation.')
parser.add_argument('-genome', '--genome', type = str, default = 'hg19', help = 'The genome to be used to create the IGV session file (hg19, hg38, mm9, mm10, etc.).')
parser.add_argument('-effect_size', '--effect_size', type = float, default = 1, help = 'Effect size to estimate statistical power.')
parser.add_argument('-estimate_power', '--estimate_statistical_power', type = str, choices = ['yes', 'no'], default = 'no', help = 'Whether or not to compute a track estimating statistical power for the CRISPR tiling screen data.')
parser.add_argument('-padjs', '--padj_cutoffs', default = 0, nargs = '+', help = 'List of p-adj. (Benjamini-Hochberg) cut-offs for determining significance of regulatory regions in the CRISPR tiling screen.')
parser.add_argument('-out_dir', '--out_directory', type = str, default = 'CRISPR_SURF_Analysis_%s' % (timestr), help = 'The name of the output directory to place CRISPR-SURF analysis files.')
args = parser.parse_args()

##### Initialize arguments
sgRNAs_summary_table = args.sgRNAs_summary_table
perturbation_type = args.perturbation_type
characteristic_perturbation_range = args.characteristic_perturbation_range
scale = args.scale
limit = args.limit
averaging_method = args.averaging_method
null_distribution = args.null_distribution
simulation_n = args.simulation_n
test_type = args.test_type
gamma_list = args.lambda_list
gamma = args.lambda_val
correlation = args.correlation
genome = args.genome
padj_cutoffs = args.padj_cutoffs
effect_size = args.effect_size
estimate_statistical_power = args.estimate_statistical_power
out_dir = args.out_directory

##### Create output directory
if not os.path.exists(out_dir):
	os.makedirs(out_dir)

##### Establish logging mechanism
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

fh = logging.FileHandler(out_dir + '/crispr-surf.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

##### Argument handling and assigning dependent defaults when necessary
logger.info('Performing argument handling and setting defaults when necessary ...')

# Setting default for characteristic perturbation range parameter
if characteristic_perturbation_range == 0:
	if (perturbation_type == 'cas9') or (perturbation_type == 'cpf1') or (perturbation_type == 'cas9_avg'):
		characteristic_perturbation_range = 7
		logger.info('The characteristic perturbation range parameter is set to default 7 bp based on Cas9/Cpf1 perturbation ...')

	elif (perturbation_type == 'crispri') or (perturbation_type == 'crispra'):
		characteristic_perturbation_range = 200
		logger.info('The characteristic perturbation range parameter is set to default 200 bp based on CRISPRi/a perturbation ...')

	elif (perturbation_type == 'be'):
		characteristic_perturbation_range = 5
		logger.info('The characteristic perturbation range parameter is set to default 5 bp based on base-editing perturbation ...')

else:
	logger.info('The characteristic perturbation range parameter is set to %s bp ...' % characteristic_perturbation_range)	

# Setting default for limit parameter
if limit == 0:
	if (perturbation_type == 'cas9') or (perturbation_type == 'cpf1') or (perturbation_type == 'cas9_avg'):
		limit = 30
		logger.info('The limit parameter is set to default 30 bp based on Cas9/Cpf1 perturbation ...')

	elif (perturbation_type == 'crispri') or (perturbation_type == 'crispra'):
		limit = 300
		logger.info('The limit parameter is set to default 300 bp based on CRISPRi/a perturbation ...')

	elif (perturbation_type == 'be'):
		limit = 10
		logger.info('The limit parameter is set to default 10 bp based on base-editing perturbation ...')

else:
	logger.info('The limit parameter is set to %s bp ...' % limit)


# Setting default for scale parameter
if scale == 0:
	if characteristic_perturbation_range >= 50:
		scale = int(int(characteristic_perturbation_range)/10)

	else:
		scale = 1

logger.info('The scale parameter is set to %s ...' % (scale))

if padj_cutoffs == 0:
	padj_cutoffs = [0.05, 0.01, 0.001, 0.0001]
	
else:
	try:
		padj_cutoffs = sorted([float(x) for x in padj_cutoffs])
		logger.info('The following p adj. values will be used %s ...' % padj_cutoffs)
	except:
		logger.error('The p adj. argument input could not be converted into a list. Example: 0.01 0.05 0.1 0.25 ...')
		sys.exit('The p adj. argument input could not be converted into a list of floats. Example: 0.01 0.05 0.1 0.25 ...')

logger.info('The following p adj. values will be used %s ...' % padj_cutoffs)

# Setting default for correlation parameter
if correlation == 0:
	if characteristic_perturbation_range <= 50:
		correlation = 0.8

	else:
		correlation = 0.9

logger.info('The correlation parameter is set to %s ...' % correlation)

# Setting default gamma list parameter
if gamma == 0:
	if gamma_list == 0:
		if characteristic_perturbation_range <= 50:
			gamma_list = [0.02,0.04,0.06,0.08,0.1,0.2,0.4,0.6,0.8,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0]
			logger.info('The lambda list parameter is set to %s based on the characteristic perturbation range parameter of %s bp ...' % (gamma_list, characteristic_perturbation_range))
		else:
			gamma_list = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0]
			logger.info('The lambda list parameter is set to %s based on the characteristic perturbation range parameter of %s bp ...' % (gamma_list, characteristic_perturbation_range))

	else:
		try:
			gamma_list = sorted([float(x) for x in gamma_list])
			logger.info('The lambda list parameter is set to %s ...' % (gamma_list))
		except:
			logger.error('The lambda list argument input could not be converted into a list. Example: 1 2 3 4 5 ...')
			sys.exit('The lambda list argument input could not be converted into a list of floats. Example: 1 2 3 4 5 ...')

else:
	try:
		gamma_list = [float(gamma)]
		logger.info('The lambda parameter is set to %s ...' % (gamma_list[0]))
	except:
		logger.error('The lambda argument input could not be converted into a float.')
		sys.exit('The lambda argument input could not be converted into a float.')

logger.info('Finished argument handling and setting defaults when necessary ...')

##### Check if sgRNA summary csv exists, load dataframe and make sure formatting is correct
try:
	df = pd.read_csv(sgRNAs_summary_table)
	required_columns = ['Chr', 'Start', 'Stop', 'Perturbation_Index', 'sgRNA_Sequence', 'Strand', 'sgRNA_Type', 'Log2FC_Replicate1']
	df_columns = df.columns.tolist()

except:
	logger.error('The input sgRNA summary file named %s does not exist and cannot be opened. Please make sure the specified location and file name is correct ...' % sgRNAs_summary_table)
	sys.exit('The input sgRNA summary file named %s does not exist and cannot be opened. Please make sure the specified location and file name is correct ...' % sgRNAs_summary_table)

if '.csv' not in sgRNAs_summary_table:
	logger.error('Please make sure input sgRNA summary file is in .CSV format ...')
	sys.exit('Please make sure input sgRNA summary file is in .CSV format ...')

if len(set(required_columns + df_columns)) == len(df_columns):
	replicates = len([x for x in df.columns.tolist() if 'Log2FC_Replicate' in x])
	logger.info('Input sgRNA summary dataframe successfully loaded ...')

else:
	logger.error('Required columns not detected in input sgRNA summary file (Chr, Start, Stop, Perturbation_Index, sgRNA_Sequence, Strand, sgRNA_Type, Log2FC_Replicate1, (add more replicate scores here)) ...')
	sys.exit('Required columns not detected in input sgRNA summary file (Chr, Start, Stop, Perturbation_Index, sgRNA_Sequence, Strand, sgRNA_Type, Log2FC_Replicate1, (add more replicate scores here)) ...')

##### Construct gaussian perturbation pattern
if perturbation_type == 'cas9_avg':
	perturbation_profile = [0.005291005,0.015873016,0.021164021,0.031746032,0.031746032,0.031746032,0.047619048,0.042328042,0.047619048,0.068783069,0.063492063,0.068783069,0.095238095,0.08994709,0.095238095,0.126984127,0.121693122,0.126984127,0.158730159,0.153439153,0.158730159,0.206349206,0.19047619,0.201058201,0.26984127,0.28042328,0.264550265,0.28042328,0.44973545,0.544973545,1,0.544973545,0.44973545,0.28042328,0.264550265,0.28042328,0.26984127,0.201058201,0.19047619,0.206349206,0.158730159,0.153439153,0.158730159,0.126984127,0.121693122,0.126984127,0.095238095,0.08994709,0.095238095,0.068783069,0.063492063,0.068783069,0.047619048,0.042328042,0.047619048,0.031746032,0.031746032,0.031746032,0.021164021,0.015873016,0.005291005]
else:
	perturbation_profile = gaussian_pattern(characteristic_perturbation_range = characteristic_perturbation_range, scale = scale, limit = limit)

##### Deconvolve underlying functional signal per biological replicate
logger.info('Deconvolving CRISPR tiling screen signal for biological replicates ...')
gammas2betas = {}

try:
	sgRNA_indices = [int(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Perturbation_Index']).flatten().tolist()]
	chromosomes = [str(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Chr']).flatten().tolist()]

except:
	logger.error('Please make sure all values in column Perturbation_Index are numbers for sgRNAs that are NOT classified as negative controls ...')
	sys.exit('Please make sure all values in column Perturbation_Index are numbers for sgRNAs that are NOT classified as negative controls ...')

for i in range(1, replicates + 1):
	logger.info('Deconvolving signal for Replicate %s' % str(i))

	try:
		observations = [float(x) for x in np.array(df.loc[df['sgRNA_Type'] != 'negative_control', 'Log2FC_Replicate' + str(i)]).flatten().tolist()]

	except:
		logger.error('Please make sure all values in column Log2FC_Replicate%s are numbers ...' % str(i))
		sys.exit('Please make sure all values in column Log2FC_Replicate%s are numbers ...' % str(i))

	try:
		gammas2betas[i], guideindices2bin = crispr_surf_deconvolution(observations = observations, chromosomes = chromosomes, sgRNA_indices = sgRNA_indices, perturbation_profile = perturbation_profile, gamma_list = gamma_list, scale = scale)
		
	except:
		logger.error('Deconvolution of Replicate %s was not successful. The scale parameter may need to be adjusted ...' % str(i))
		sys.exit('Deconvolution of Replicate %s was not successful. The scale parameter may need to be adjusted ...' % str(i))

logger.info('Successfully deconvolved all biological replicates ...')

##### Find gamma from correlation ratio
if len(gamma_list) == 1:
	gamma_use = gamma_list[0]

else:
	try:
		(gamma_range, gamma_use) = crispr_surf_find_lambda(gammas2betas = gammas2betas, correlation_start = correlation, correlation_stop = correlation, correlation_opt = correlation, out_dir = out_dir)
		logger.info('Identified lambda range to be used for downstream deconvolution statistics')

	except:
		logger.error('Failed to identify lambda range ...')
		sys.exit('Failed to identify lambda range ...')

##### Combine biological replicate deconvolved signals
try:
	gammas2betas_updated = crispr_surf_deconvolved_signal(gammas2betas = gammas2betas, gamma_chosen = gamma_use, averaging_method = averaging_method, out_dir = out_dir)
	logger.info('Combined all biological replicates into single deconvolution profile for downstream statistics ...')

except:
	logger.error('Combined all biological replicates into single deconvolution profile for downstream statistics ...')
	sys.exit('Combined all biological replicates into single deconvolution profile for downstream statistics ...')

##### Bootstrap deconvolution analysis to assign statistical significance
try:
	gammas2betas_updated, replicate_parameters = crispr_surf_statistical_significance(sgRNA_summary_table = sgRNAs_summary_table, sgRNA_indices = sgRNA_indices, perturbation_profile = perturbation_profile, gammas2betas = gammas2betas_updated, null_distribution = null_distribution, simulation_n = simulation_n, test_type = test_type, guideindices2bin = guideindices2bin, averaging_method = averaging_method, padj_cutoffs = padj_cutoffs, effect_size = effect_size, limit = limit, scale = scale, estimate_statistical_power = estimate_statistical_power)
	logger.info('Finished simulations to assess statistical significance of deconvolution profile ...')

except:
	logger.error('Simulations for statistical significance failed ...')
	sys.exit('Simulations for statistical significance failed ...')

##### Update sgRNA summary file
try:
	crispr_surf_sgRNA_summary_table_update(sgRNA_summary_table = sgRNAs_summary_table, gammas2betas = gammas2betas_updated, averaging_method = averaging_method, scale = scale, guideindices2bin = guideindices2bin, simulation_n = simulation_n, padj_cutoffs = padj_cutoffs, out_dir = out_dir)
	logger.info('Successfully updated sgRNA summary table ...')

except:
	logger.error('Failed to update sgRNA summary table ...')
	sys.exit('Failed to update sgRNA summary table ...')

##### Output beta profile
try:
	complete_beta_profile(gammas2betas = gammas2betas_updated, simulation_n = simulation_n, padj_cutoffs = padj_cutoffs, estimate_statistical_power = estimate_statistical_power, out_dir = out_dir)
	logger.info('Successfully created beta profile file ...')

except:
	logger.error('Failed to create beta profile file ...')
	sys.exit('Failed to create beta profile file ...')

##### Output significant regions file
try:
	crispr_surf_significant_regions(sgRNA_summary_table = sgRNAs_summary_table.split('/')[-1].replace('.csv', '_updated.csv'), gammas2betas = gammas2betas_updated, padj_cutoffs = padj_cutoffs, scale = scale, guideindices2bin = guideindices2bin, out_dir = out_dir)
	logger.info('Successfully created significant regions file ...')

except:
	logger.error('Failed to create significant regions file ...')
	sys.exit('Failed to create significant regions file ...')

##### Output IGV tracks
try:
	crispr_surf_IGV(sgRNA_summary_table = sgRNAs_summary_table.split('/')[-1].replace('.csv', '_updated.csv'), gammas2betas = gammas2betas_updated, padj_cutoffs = padj_cutoffs, genome = genome, scale = scale, guideindices2bin = guideindices2bin, estimate_statistical_power = estimate_statistical_power, out_dir = out_dir)
	logger.info('Successfully created IGV tracks ...')

except:
	logger.error('Failed to create IGV tracks ...')
	sys.exit('Failed to create IGV tracks ...')

##### Write parameters file
parameters = {
'sgRNAs_summary_table': sgRNAs_summary_table,
'perturbation_type': perturbation_type,
'characteristic_perturbation_range': characteristic_perturbation_range,
'scale': scale,
'limit': limit,
'averaging_method': averaging_method,
'null_distribution': null_distribution,
'replicate_parameters [loc/scale]': ' '.join(map(str, replicate_parameters)).replace(', ', '/').replace(' ', ','),
'simulation_n': simulation_n,
'test_type': test_type,
'lambda_list': ' '.join(map(str, gamma_list)),
'lambda': gamma,
'lambda_used': gamma_use,
'correlation': correlation,
'genome': genome,
'effect_size': effect_size,
'estimate_statistical_power': estimate_statistical_power,
'padj_cutoffs': ' '.join(map(str, padj_cutoffs)),
'out_directory': out_dir
}

with open(out_dir + '/crispr-surf_parameters.csv', 'w') as f:
	f.write(','.join(map(str, ['Argument','Value'])) + '\n')
	for parameter in sorted(parameters):
		f.write(','.join(map(str, [parameter, parameters[parameter]])) + '\n')