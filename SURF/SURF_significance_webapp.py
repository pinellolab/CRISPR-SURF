### Call crispr_surf_deconvolution function through command line ###

import argparse
import ast
import json

from CRISPR_SURF_Core_WebApp import crispr_surf_deconvolution_simulations, crispr_surf_statistical_power, crispr_surf_deconvolved_signal, crispr_surf_statistical_significance

parser = argparse.ArgumentParser(description = 'CRISPR-SURF statistical significance command line tool for webapp.')
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
print '1 STARTING TO COMBINE SIGNAL'
gammas2betas_updated = crispr_surf_deconvolved_signal(gammas2betas = data_dict['gammas2betas'], gamma_chosen = data_dict['gamma_use'], averaging_method = data_dict['avg'], uploads_dir = UPLOADS_FOLDER, results_dir = RESULTS_FOLDER)
print '2 FINISHED COMBINING SIGNAL'

# hack
json_good = False
while not json_good:
    with open(UPLOADS_FOLDER + '/data2.json', 'r') as f:
        json_string = f.readline().strip()
        try:
            data_dict2 = json.loads(json_string)
            json_good = True
        except:
            pass

print '3 STARTING TO FIND REGIONS'
gammas2betas_updated, replicate_parameters = crispr_surf_statistical_significance(sgRNA_summary_table = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv', sgRNA_indices = data_dict['sgRNA_indices'], perturbation_profile = data_dict['perturbation_profile'], gammas2betas = data_dict2['gammas2betas'], simulation_type = data_dict['sim_type'], simulation_n = data_dict['sim_n'], guideindices2bin = data_dict['guideindices2bin'], averaging_method = data_dict['avg'], padj_cutoffs = [0.05], effect_size = data_dict['effect_size'], limit = data_dict['limit'], scale = data_dict['scale'], rapid_mode = data_dict['rapid'], uploads_dir = UPLOADS_FOLDER, results_dir = RESULTS_FOLDER)
print '4 FINISHED FINDING REGIONS'