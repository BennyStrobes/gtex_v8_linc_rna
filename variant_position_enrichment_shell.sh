#!/bin/bash -l

#SBATCH
#SBATCH --time=9:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


rare_variant_dir="$1"
variant_enrichment_dir="$2"
variant_position_enrichment_dir="$3"
outlier_calls_dir="$4"
european_ancestry_individual_list="$5"
gencode_gene_annotation_file="$6"
cluster_info_file="$7"
exon_file="$8"
tissue_names_file="$9"



# Whether to take 'top_outlier' per cluster of 'all' variants per cluster
enrichment_version="all"


##################################
# Tissue-by Tissue variant outlier enrichment
##################################
if false; then
# Range of Distances
distances=( "2" "4" "6" "8" "10" "100" "1000")
# Range of pvalue thresholds
pvalue_thresholds=( ".000001" ".00001" ".0001")
# Loop through distances
for distance_window in "${distances[@]}"; do
	# Loop through range of pvalue thresholds
	for pvalue_threshold in "${pvalue_thresholds[@]}"; do
		echo "Distance window="$distance_window" pvalue_threshold="$pvalue_threshold
		# Run Tissue by tissue enrichment of rare variants within outlier calls (after filtering individuals that were "global outliers")
		variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"

		# For protein coding
		gene_type="proteincoding"
		outlier_suffix="_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_"$gene_type".txt"
		output_root=$variant_enrichment_dir"tbt_variant_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version"_gene_type_"$gene_type
		python variant_tbt_enrichment_quantification.py $variant_bed_file $output_root $outlier_calls_dir $outlier_suffix $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold

		# For lincrna
		gene_type="lincrna"
		outlier_suffix="_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_"$gene_type".txt"
		output_root=$variant_enrichment_dir"tbt_variant_outlier_enrichment_pvalue_"$pvalue_threshold"_distance_"$distance_window"_version_"$enrichment_version"_gene_type_"$gene_type
		python variant_tbt_enrichment_quantification.py $variant_bed_file $output_root $outlier_calls_dir $outlier_suffix $european_ancestry_individual_list $tissue_names_file $enrichment_version $pvalue_threshold

	done
done
fi



##################################
# Cross Tissue variant outlier enrichment
##################################
previous_distance="2"
distances=( "2" "4" "6" "8" "10" "100" "1000")
# Loop through distances
for distance_window in "${distances[@]}"; do
	echo "Multi-tissue distance="$distance_window
	if false; then
	# Run cross tissue enrichment using this sized window for LINCRNA
	splicing_outlier_file=$outlier_calls_dir"cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_lincrna.txt"
	variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
	previous_variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$previous_distance".txt"
	output_root=$variant_enrichment_dir"cross_tissue_variant_outlier_enrichment_mutually_exclusive_distance_"$distance_window"_version_"$enrichment_version"_gene_type_lincrna"
	python variant_cross_tissue_mutually_exclusive_distance_enrichment_quantification.py $variant_bed_file $previous_variant_bed_file $output_root $splicing_outlier_file $european_ancestry_individual_list $enrichment_version
	
	# Run cross tissue enrichment using this sized window for PROTEINCODING
	splicing_outlier_file=$outlier_calls_dir"cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_proteincoding.txt"
	variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
	previous_variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$previous_distance".txt"
	output_root=$variant_enrichment_dir"cross_tissue_variant_outlier_enrichment_mutually_exclusive_distance_"$distance_window"_version_"$enrichment_version"_gene_type_proteincoding"
	python variant_cross_tissue_mutually_exclusive_distance_enrichment_quantification.py $variant_bed_file $previous_variant_bed_file $output_root $splicing_outlier_file $european_ancestry_individual_list $enrichment_version
	fi

	# Run cross tissue enrichment using this sized window for LINCRNA
	splicing_outlier_file=$outlier_calls_dir"cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_lincrna.txt"
	variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
	output_root=$variant_enrichment_dir"cross_tissue_variant_outlier_enrichment_distance_"$distance_window"_version_"$enrichment_version"_gene_type_lincrna"
	python variant_cross_tissue_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_file $european_ancestry_individual_list $enrichment_version

	# Run cross tissue enrichment using this sized window for PROTEINCODING
	splicing_outlier_file=$outlier_calls_dir"cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_proteincoding.txt"
	variant_bed_file=$rare_variant_dir"variant_cluster_only_bed_"$distance_window".txt"
	output_root=$variant_enrichment_dir"cross_tissue_variant_outlier_enrichment_distance_"$distance_window"_version_"$enrichment_version"_gene_type_proteincoding"
	python variant_cross_tissue_enrichment_quantification.py $variant_bed_file $output_root $splicing_outlier_file $european_ancestry_individual_list $enrichment_version



	previous_distance=$distance_window
done

