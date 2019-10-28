######## Ben Strober
######## 11/27/18




#######################################
# Input Data
#######################################

# Ordered list of GTEx v8 tissue names
tissue_names_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtex_v8_tissues.txt"

# File containing list of genes (ensamble IDs) to be used in this analysis
# File provided by Nicoloe Ferraro
gtex_v8_gene_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gencode.v26.GRCh38.genes_genetypes_autosomal_PCandlinc_only.txt"

# File containing list of lincRNAs to be used in this analysis
linc_rna_gene_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8_linc_rna/input_data/gtexV8_lincRNA_genes.txt"

# Gencode Gene annotation file
# Downloaded from scp.nygenome.org:/data/delivery/gtex-rare/data/gencode.v26.GRCh38.genes.gtf on 1/10/19
# Xin put the file there. He said he got in from the GTEx Exchange
gencode_gene_annotation_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gencode.v26.GRCh38.genes.gtf"

# List of individuals used in each tissue generated by Nicole Ferraro
# Downloaded from scp.nygenome.org:/data/delivery/gtex-rare/data/gtexV8_inds_by_tissue.txt on 11/30/18
individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/gtexV8_inds_by_tissue.txt"

# List of european ancestry individuals (extracted from Nicole Ferraro's variant list: "all_rare_variants_SNPs_10kb_genebody.txt")
european_ancestry_individual_list="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/input_data/variant_calls/european_ancestry_individuals.txt"

# Directory containing outlier calls for gtex v8
gtex_v8_outlier_calls_dir="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/splicing_outlier_calls/"

# File containing all cluster_ids and their corresponding junction positions
# Generated with "outlier_calling" (https://github.com/BennyStrobes/gtex_v8_rare_splice/tree/master/outlier_calling)
cluster_info_file="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8/splicing/outlier_calling/clusters_filtered/cluster_info.txt"

#############################################################
#Used Directories (directories need to be created and empty before starting)
#############################################################

# Root directories that all other directories will be placed in 
output_root="/work-zfs/abattle4/bstrober/rare_variant/gtex_v8_linc_rna/"

# Directory containing lincRNA specific outlier calls
outlier_calls_dir=$output_root"outlier_calls/"










#############################################################
# Scripts: Run the following X parts in Series
#############################################################

sh generate_linc_rna_only_outlier_files.sh $tissue_names_file $gtex_v8_gene_list $linc_rna_gene_list $gtex_v8_outlier_calls_dir $cluster_info_file $outlier_calls_dir

