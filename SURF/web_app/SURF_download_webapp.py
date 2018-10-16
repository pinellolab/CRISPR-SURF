### Call crispr_surf_deconvolution function through command line ###

import argparse
import ast
import json

from CRISPR_SURF_Core_WebApp import crispr_surf_sgRNA_summary_table_update, complete_beta_profile, crispr_surf_significant_regions, crispr_surf_IGV

parser = argparse.ArgumentParser(description = 'CRISPR-SURF download command line tool for webapp.')
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

# Update sgRNA summary table
print '1 UPDATING SUMMARY TABLE'
crispr_surf_sgRNA_summary_table_update(sgRNA_summary_table = UPLOADS_FOLDER + '/sgRNAs_summary_table.csv', gammas2betas = data_dict2['gammas2betas'], averaging_method = data_dict['avg'], scale = data_dict['scale'], guideindices2bin = data_dict['guideindices2bin'], simulation_n = data_dict['sim_n'], padj_cutoffs = [data_dict2['fdr']], out_dir = RESULTS_FOLDER)
print '2 FINISHED UPDATING SUMMARY TABLE'

# Output complete beta profile
print '1 OUTPUTTING BETAS'
complete_beta_profile(gammas2betas = data_dict2['gammas2betas'], simulation_n = data_dict['sim_n'], padj_cutoffs = [data_dict2['fdr']], out_dir = RESULTS_FOLDER)
print '2 FINISHED OUTPUTTING BETAS'

# Output significant regions
print '1 OUTPUTTING SIGNIFICANT REGIONS'
crispr_surf_significant_regions(sgRNA_summary_table = 'sgRNAs_summary_table_updated.csv', gammas2betas = data_dict2['gammas2betas'], padj_cutoffs = [data_dict2['fdr']], scale = data_dict['scale'], guideindices2bin = data_dict['guideindices2bin'], out_dir = RESULTS_FOLDER)
print '2 FINISHED OUTPUTTING SIGNIFICANT REGIONS'

# Output IGV tracks
print '1 OUTPUTTING IGV TRACKS'
crispr_surf_IGV(sgRNA_summary_table = 'sgRNAs_summary_table_updated.csv', gammas2betas = data_dict2['gammas2betas'], padj_cutoffs = [data_dict2['fdr']], genome = data_dict['genome'], scale = data_dict['scale'], guideindices2bin = data_dict['guideindices2bin'], out_dir = RESULTS_FOLDER)
print '2 FINISHED OUTPUTTING IGV TRACKS'

##### Write parameters file
parameters = {
'sgRNAs_summary_table': UPLOADS_FOLDER + '/sgRNAs_summary_table.csv',
'perturbation_type': data_dict['pert'],
'characteristic_perturbation_range': data_dict['range'],
'scale': data_dict['scale'],
'limit': data_dict['limit'],
'averaging_method': data_dict['avg'],
'null_distribution': data_dict['null_distribution'],
'simulation_n': data_dict['sim_n'],
'test_type': data_dict['test_type'],
'gamma_list': ' '.join(map(str, data_dict['gamma_list'])),
'gamma_used': data_dict['gamma_use'],
'correlation': data_dict['correlation'],
'genome': data_dict['genome'],
'effect_size': data_dict['effect_size'],
'padj_cutoffs': ' '.join(map(str, [data_dict2['fdr']])),
'rapid_mode': data_dict['rapid'],
'out_directory': RESULTS_FOLDER
}

with open(RESULTS_FOLDER + '/crispr-surf_parameters.csv', 'w') as f:
    f.write(','.join(map(str, ['Argument','Value'])) + '\n')
    for parameter in sorted(parameters):
        f.write(','.join(map(str, [parameter, parameters[parameter]])) + '\n')