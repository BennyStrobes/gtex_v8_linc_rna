import numpy as np 
import os
import sys
import pdb


def get_dictionary_lists_of_linc_rna_and_protein_coding_genes(gtex_v8_gene_list, linc_rna_gene_list):
	linc_rna = {}
	pc = {}
	# First extract linc rnas
	f = open(linc_rna_gene_list)
	for line in f:
		line = line.rstrip()
		linc_rna[line] = 0
	f.close()
	# Next protein coding
	f = open(gtex_v8_gene_list)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if data[1] == 'protein_coding':
			pc[data[0]] = 1
	f.close()
	return linc_rna, pc

def get_dictionary_lists_of_linc_rna_and_protein_coding_clusters(linc_rna_genes, pc_genes, cluster_info_file):
	linc_rna_clusters = {}
	pc_clusters = {}
	ambiguous_clusters = {}
	f = open(cluster_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# skip header line
		if head_count == 0:
			head_count = head_count + 1
			continue
		cluster_id = data[0]
		ensamble_ids = data[2].split(',')

		linc_rna_binary = 0
		pc_binary = 0
		for ensamble_id in ensamble_ids:
			if ensamble_id in linc_rna_genes:
				linc_rna_binary = 1
			elif ensamble_id in pc_genes:
				pc_binary = 1
			else:
				print('assumption eroror!')
		if linc_rna_binary == 0 and pc_binary == 0:
			print('assumption error!')
			pdb.set_trace()
		elif linc_rna_binary == 1 and pc_binary == 0:
			linc_rna_clusters[cluster_id] = 1
		elif linc_rna_binary == 0 and pc_binary == 1:
			pc_clusters[cluster_id] = 1
		elif linc_rna_binary == 1 and pc_binary == 1:
			ambiguous_clusters[cluster_id] = 1
		else:
			print('assumption eororor!')
			pdb.set_trace()
	f.close()
	return linc_rna_clusters, pc_clusters, ambiguous_clusters

def get_list_of_tissue_names(tissue_names_file):
	f = open(tissue_names_file)
	tissues = []
	for line in f:
		line = line.rstrip()
		tissues.append(line)
	return tissues

def filter_outlier_file_with_gene_list(input_outlier_file, output_outlier_file, gene_list):
	f = open(input_outlier_file)
	t = open(output_outlier_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		gene_id = data[0]
		if gene_id not in gene_list:
			continue
		t.write(line + '\n')
	f.close()
	t.close()

tissue_names_file = sys.argv[1]
gtex_v8_gene_list = sys.argv[2]
linc_rna_gene_list = sys.argv[3]
gtex_v8_outlier_calls_dir = sys.argv[4]
cluster_info_file = sys.argv[5]
output_dir = sys.argv[6]

# Extact list of tissue names
tissue_names = get_list_of_tissue_names(tissue_names_file)

# Extract dictionary lists of linc RNA genes and protein coding genes
linc_rna_genes, pc_genes = get_dictionary_lists_of_linc_rna_and_protein_coding_genes(gtex_v8_gene_list, linc_rna_gene_list)

# Extract dictionary lists of LeafCutter clusters that map exclusively to linc RNA genes, and leafcutter clusters that map exclusively to protein coding genes
# If cluster maps to both a lincRNA and a protein coding then ignore it (ambiguous clusters)
linc_rna_clusters, pc_clusters, ambiguous_clusters = get_dictionary_lists_of_linc_rna_and_protein_coding_clusters(linc_rna_genes, pc_genes, cluster_info_file)


# Split up cluster level outlier files into lincRNA and protein coding files
for tissue in tissue_names:
	input_outlier_file = gtex_v8_outlier_calls_dir + tissue + '_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue.txt'
	# linc RNA
	linc_outlier_file = output_dir + tissue + '_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_lincrna.txt'
	filter_outlier_file_with_gene_list(input_outlier_file, linc_outlier_file, linc_rna_clusters)
	# PC
	pc_outlier_file = output_dir + tissue + '_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_proteincoding.txt'
	filter_outlier_file_with_gene_list(input_outlier_file, pc_outlier_file, pc_clusters)

# Split up gene level outlier files into lincRNA and protein coding files
for tissue in tissue_names:
	input_outlier_file = gtex_v8_outlier_calls_dir + tissue + '_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_gene_level.txt'
	# linc RNA
	linc_outlier_file = output_dir + tissue + '_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_gene_level_lincrna.txt'
	filter_outlier_file_with_gene_list(input_outlier_file, linc_outlier_file, linc_rna_genes)
	# PC
	pc_outlier_file = output_dir + tissue + '_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_gene_level_proteincoding.txt'
	filter_outlier_file_with_gene_list(input_outlier_file, pc_outlier_file, pc_genes)

# Split up global outliers for Leafcutter Cluster
input_outlier_file = gtex_v8_outlier_calls_dir + 'cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue.txt'
# linc RNA
linc_outlier_file = output_dir + 'cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_lincrna.txt'
filter_outlier_file_with_gene_list(input_outlier_file, linc_outlier_file, linc_rna_clusters)
# PC
pc_outlier_file = output_dir + 'cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_proteincoding.txt'
filter_outlier_file_with_gene_list(input_outlier_file, pc_outlier_file, pc_clusters)


# Split up global outliers for genes
input_outlier_file = gtex_v8_outlier_calls_dir + 'cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'
# linc RNA
linc_outlier_file = output_dir + 'cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level_lincrna.txt'
filter_outlier_file_with_gene_list(input_outlier_file, linc_outlier_file, linc_rna_genes)
# PC
pc_outlier_file = output_dir + 'cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level_proteincoding.txt'
filter_outlier_file_with_gene_list(input_outlier_file, pc_outlier_file, pc_genes)