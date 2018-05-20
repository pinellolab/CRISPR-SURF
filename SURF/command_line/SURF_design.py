# Find sgRNAs within specified genomic ranges

##### Import libraries
import os
import sys
import re
import argparse
import logging
from bioutilities import Genome_2bit, Coordinate

##### Argument handeling
parser = argparse.ArgumentParser(description = 'Find sgRNAs across genomic regions for CRISPR tiling screen')
parser.add_argument('-f', '--f', required = True, type = str, help = 'Input bed file to make DNA excisions.')
parser.add_argument('-genome', '--genome', required = True, type = str, help = 'Input genome 2bit file.')
parser.add_argument('-pams', '--pams', required = True, type = str, nargs = '+', default = [], help = 'Specification of different CRISPR PAMs ([ATCG]GG, TTT[ACG, etc.]).')
parser.add_argument('-orient', '--orientations', required = True, choices = ['left', 'right'], type = str, nargs = '+', default = [], help = 'Orientation of spacer relative to PAM (Cas9 -> left, Cpf1 -> right).')
parser.add_argument('-guide_l', '--guide_length', default = 20, type = int, help = 'Length of sgRNA')
parser.add_argument('-g_constraint', '--g_constraint', default = 'false', choices = ['true', 'false', 'True', 'False'], help = "Constraint forcing the 5' sgRNA bp to be G base.")
parser.add_argument('-out', '--out_dir', default = '.', type = str, help = 'Output directory.')

args = parser.parse_args()

##### Initialize arguments
input_f = args.f
genome = args.genome
pams = args.pams
spacer_orientations = args.orientations
gRNA_length = args.guide_length
g_constraint = args.g_constraint
out_dir = args.out_dir

if not os.path.exists(out_dir):
	os.makedirs(out_dir)

# Initialize logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

fh = logging.FileHandler(out_dir + '/surf_design.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

for pam in pams:
	if not all(x in ['A', 'T', 'C', 'G', '[', ']'] for x in set(list(pam))):
		logger.error('The PAM %s includes an unidentifiable character. Please only use A, T, C, G bases.\nBrackets [] can be used to specify multiple bases at a single position.\nFor example: NGG = [ATCG]GG ...' % pam)
		sys.exit(1)

##### Input genome
logger.info('Importing %s genome ...' % genome)
genome = Genome_2bit(genome)

##### Input regions file
logger.info('Reading in %s file ...' % input_f)

total_regions = 0
bed_dict = {}
with open(input_f, 'r') as f:
	for line in f:

		total_regions += 1

		if 'csv' in input_f:
			line = line.strip('\n').split(',')
		else:
			line = line.strip('\n').split()

		if len(line) < 3:
			logger.error('The input file %s does not have at least three columns (chr, start, stop).')
			sys.exit(1)

		chrom, start, stop = line[:3]
		if len(line[3:]) > 0:
			info = '_'.join(map(str, line[3:]))
			id_name = chrom + '_' + start + '_' + stop + '_' + info
		else:
			id_name = chrom + '_' + start + '_' + stop

		if chrom not in bed_dict:
			bed_dict[chrom] = []

		bed_dict[chrom].append([int(start), int(stop), id_name])

##### Reverse complement function
def reverse_complement(sequence):
	sequence = sequence.upper()
	new_sequence = ''
	for base in sequence:
		if base == 'A':
			new_sequence += 'T'
		elif base == 'T':
			new_sequence += 'A'
		elif base == 'C':
			new_sequence += 'G'
		elif base == 'G':
			new_sequence += 'C'
		elif base == '[':
			new_sequence += ']'
		elif base == ']':
			new_sequence += '['
	return(new_sequence[::-1])

def pam_length_calc(pam):
	pam_l = 0
	pam_split_1 = pam.split('[')
	for i in pam_split_1:
		if ']' in i:
			pam_l += 1
			pam_split_2 = i.split(']')
			pam_l += len(pam_split_2[-1])
		else:
			pam_l += len(i)
	return pam_l

##### PAM site initialization and position weights dependent on nuclease selection
pam_list = []

logger.info('Searching for the following PAMs:')
pam_headers = []
for i in range(len(pams)):

	logger.info('---------- PAM Group %s ---------' % (i + 1))

	pam_groups = pams[i].split(',')
	pam_headers.append(pams[i].replace(',', '-'))
	pam_entry = []

	for pam in pam_groups:
		pam_entry.append([pam, reverse_complement(pam), pam_length_calc(pam), spacer_orientations[i]])

		if spacer_orientations[i] == 'right':
			logger.info('(+)strand: %s   (-)strand: %s   Length: %s   sgRNA: %s' % (pam, reverse_complement(pam), pam_length_calc(pam), pam + 'N'*int(gRNA_length)))
		elif spacer_orientations[i] == 'left':
			logger.info('(+)strand: %s   (-)strand: %s   Length: %s   sgRNA: %s' % (pam, reverse_complement(pam), pam_length_calc(pam), 'N'*int(gRNA_length) + pam))

	pam_list.append(pam_entry)

##### Initialize data storage for output
total_sgRNAs_dict = {}
pams_found = {}

##### Perform gRNA search
logger.info('Searching for PAMs for each genomic site ...')

counter = 0
for chrom in bed_dict:
	for site in bed_dict[chrom]:

		counter += 1

		if counter%1000 == 0:
			logger.info('Completed PAM search for %s out of %s sites ...' % (counter, total_regions))

		pams_found_list = []

		start, stop, sample = site
		total_sgRNAs_dict[sample] = []

		pams_found[sample] = [0]*len(pam_list)

		genomic_region = genome.extract_sequence(Coordinate(chrom, start, stop)).upper()

		for i in range(len(pam_list)):

			for pam_entry in pam_list[i]:

				find_guides_top = [[m.start()] for m in re.finditer('(?=%s)' % pam_entry[0], genomic_region, re.IGNORECASE)]
				find_guides_bottom = [[m.start()] for m in re.finditer('(?=%s)' % pam_entry[1], genomic_region, re.IGNORECASE)]

				if find_guides_top:

					for match in find_guides_top:

						start_g = match[0]
						stop_g = start_g + pam_entry[2]

						if pam_entry[3] == 'right':
							pam_seq = genomic_region[start_g:stop_g]
							gRNA_seq = genomic_region[stop_g:stop_g + gRNA_length]

							start_g_adj = start + start_g - 1
							stop_g_adj = start + stop_g + gRNA_length - 1

							if (len(gRNA_seq) + len(pam_seq)) == (gRNA_length + pam_entry[2]):

								if g_constraint.lower() == 'true' and gRNA_seq[0] != 'G':
									continue

								pams_found[sample][i] += 1
								entry = [chrom, start_g_adj, stop_g_adj, gRNA_seq, pam_seq, '+', sample]
								total_sgRNAs_dict[sample].append(entry)


						elif pam_entry[3] == 'left':
							pam_seq = genomic_region[start_g:stop_g]
							gRNA_seq = genomic_region[start_g - gRNA_length:start_g].upper()

							start_g_adj = start + start_g - gRNA_length - 1
							stop_g_adj = start + stop_g - 1

							if (len(gRNA_seq) + len(pam_seq)) == (gRNA_length + pam_entry[2]):

								if g_constraint.lower() == 'true' and gRNA_seq[0] != 'G':
									continue

								pams_found[sample][i] += 1
								entry = [chrom, start_g_adj, stop_g_adj, gRNA_seq, pam_seq, '+', sample]
								total_sgRNAs_dict[sample].append(entry)

				if find_guides_bottom:

					for match in find_guides_bottom:

						start_g = match[0]
						stop_g = start_g + pam_entry[2]

						if pam_entry[3] == 'right':
							pam_seq = reverse_complement(genomic_region[start_g:stop_g])
							gRNA_seq = reverse_complement(genomic_region[start_g - gRNA_length:start_g])

							start_g_adj = start + start_g - gRNA_length - 1
							stop_g_adj = start + stop_g - 1

							if (len(gRNA_seq) + len(pam_seq)) == (gRNA_length + pam_entry[2]):

								if g_constraint.lower() == 'true' and gRNA_seq[0] != 'G':
									continue

								pams_found[sample][i] += 1
								entry = [chrom, start_g_adj, stop_g_adj, gRNA_seq, pam_seq, '-', sample]
								total_sgRNAs_dict[sample].append(entry)


						elif pam_entry[3] == 'left':
							pam_seq = reverse_complement(genomic_region[start_g:stop_g])
							gRNA_seq = reverse_complement(genomic_region[stop_g:stop_g + gRNA_length])

							start_g_adj = start + start_g - 1
							stop_g_adj = start + stop_g + gRNA_length - 1

							if (len(gRNA_seq) + len(pam_seq)) == (gRNA_length + pam_entry[2]):

								if g_constraint.lower() == 'true' and gRNA_seq[0] != 'G':
									continue

								pams_found[sample][i] += 1
								entry = [chrom, start_g_adj, stop_g_adj, gRNA_seq, pam_seq, '-', sample]
								total_sgRNAs_dict[sample].append(entry)

logger.info('Completed PAM search for %s out of %s sites ...' % (counter, total_regions))

##### Output single gRNAs stats
logger.info('Outputting identified sgRNAs for each genomic site ...')

sgRNAs_out = 'sgRNAs_from_PAMs_%s_on_%s.csv' % ('_'.join(map(str, pam_headers)), input_f.split('.')[0])
sgRNAs_bed_out = 'sgRNAs_from_PAMs_%s_on_%s.bed' % ('_'.join(map(str, pam_headers)), input_f.split('.')[0])

with open(out_dir + '/%s' % sgRNAs_out, 'w') as f1, open(out_dir + '/%s' % sgRNAs_bed_out, 'w') as f2:
	f1.write(','.join(map(str, ['chr', 'start', 'stop', 'sgRNA_sequence', 'pam', 'strand', 'target'])) + '\n')
	for sample in total_sgRNAs_dict:
		for sgRNA_entry in total_sgRNAs_dict[sample]:
			f1.write(','.join(map(str, sgRNA_entry)) + '\n')
			f2.write('\t'.join(map(str, sgRNA_entry)) + '\n')

##### Output counts for each PAM
logger.info('Outputting PAM counts across entire BED file ...')

pam_comparison = 'pam_comparison_%s_on_%s.csv' % ('_'.join(map(str, pam_headers)), input_f.split('.')[0])

with open(out_dir + '/%s' % pam_comparison, 'w') as f:
	f.write('Site' + ',' + ','.join(map(str, pam_headers)) + '\n')
	for site in sorted(pams_found):
		f.write(str(site) + ',' + ','.join(map(str, pams_found[site])) + '\n')