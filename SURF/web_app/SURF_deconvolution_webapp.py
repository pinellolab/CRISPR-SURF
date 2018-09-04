### Call crispr_surf_deconvolution function through command line ###

import argparse
import ast
import json

from CRISPR_SURF_Core_WebApp import crispr_surf_deconvolution, crispr_surf_find_lambda

parser = argparse.ArgumentParser(description = 'CRISPR-SURF deconvolution command line tool for webapp.')
parser.add_argument('-uploads_dir', '--uploads_directory', type = str)
parser.add_argument('-results_dir', '--results_directory', type = str)
args = parser.parse_args()

##### Initialize arguments
UPLOADS_FOLDER = str(args.uploads_directory)
RESULTS_FOLDER = str(args.results_directory)

# hack
json_good = False
while not json_good:
    with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
        json_string = f.readline().strip()
        try:
            data_dict = json.loads(json_string)
            json_good = True
        except:
            pass

# Run CRISPR-SURF deconvolution
print '1 STARTING TO DECONVOLVE'
gammas2betas, guideindices2bin = crispr_surf_deconvolution(observations = data_dict['observations'], chromosomes = data_dict['chromosomes'], sgRNA_indices = data_dict['sgRNA_indices'], perturbation_profile = data_dict['perturbation_profile'], gamma_list = data_dict['gamma_list'], scale = data_dict['scale'], uploads_dir = UPLOADS_FOLDER, results_dir = RESULTS_FOLDER)
print '2 FINISHED DECONVOLVING'

# hack
json_good = False
while not json_good:
    with open(UPLOADS_FOLDER + '/data.json', 'r') as f:
        json_string = f.readline().strip()
        try:
            data_dict = json.loads(json_string)
            json_good = True
        except:
            pass

print '3 STARTING TO FIND GAMMA RANGE'
gamma_range, gamma_use = crispr_surf_find_lambda(gammas2betas = data_dict['gammas2betas'], correlation_ratio_start = 0.5, correlation_ratio_stop = 1.0, correlation_ratio_opt = 0.8, uploads_dir = UPLOADS_FOLDER, results_dir = RESULTS_FOLDER)
print '4 FINISHED FINDING GAMMA RANGE'