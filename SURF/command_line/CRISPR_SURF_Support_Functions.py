#!/usr/bin/env python
import math

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