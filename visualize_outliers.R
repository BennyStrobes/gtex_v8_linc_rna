args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)




gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}



bar_plot_showing_number_of_clusters_per_tissue <- function(tissues, outlier_calls_dir, rna_type) {
	tissue_name_arr <- c()
	num_cluster_arr <- c()
	for (tissue_iter in 1:length(tissues)) {
		tissue_name <- tissues[tissue_iter]
		outlier_file <- paste0(outlier_calls_dir, tissue_name, "_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_", rna_type, ".txt")
		outlier_data <- read.table(outlier_file, header=TRUE)
		num_clusters <- dim(outlier_data)[1]
		tissue_name_arr <- c(tissue_name_arr, tissue_name)
		num_cluster_arr <- c(num_cluster_arr, num_clusters)
	}
	df <- data.frame(tissue=tissue_name_arr, number_of_clusters=num_cluster_arr)
	
	plotter <- ggplot(df,aes(x=tissue, y=number_of_clusters)) + 
    	geom_bar(stat="identity") +
    	labs(x="",y="Number of clusters", title=rna_type) + 
    	gtex_v8_figure_theme() +
    	theme(axis.text.x = element_text(angle=90,hjust=1)) +
    	theme(plot.title = element_text(hjust = 0.5))
    return(plotter)
}

box_plot_showing_number_of_outlier_clusters_per_tissue <- function(tissues, outlier_calls_dir, rna_type) {
	tissue_name_arr <- c()
	num_outlier_arr <- c()
	pvalue_arr <- c()

	thresholds <- c(.0001, .00001, .000001)

	for (tissue_iter in 1:length(tissues)) {
		tissue_name <- tissues[tissue_iter]
		outlier_file <- paste0(outlier_calls_dir, tissue_name, "_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_", rna_type, ".txt")
		outlier_data <- read.table(outlier_file, header=TRUE)
		outlier_data <- outlier_data[,2:(dim(outlier_data)[2])]

		for (threshold_iter in 1:length(thresholds)) {
			threshold <- thresholds[threshold_iter]
			num_outliers <- sum(rowSums(outlier_data < threshold) != 0.0)

			tissue_name_arr <- c(tissue_name_arr, tissue_name)
			num_outlier_arr <- c(num_outlier_arr, num_outliers)
			pvalue_arr <- c(pvalue_arr, threshold)
		}
	}
	df <- data.frame(tissue=tissue_name_arr, number_of_outliers=num_outlier_arr, pvalue_threshold=factor(pvalue_arr, levels=c("1e-04", "1e-05", "1e-06")))
	
	plotter <- ggplot(df,aes(x=pvalue_threshold, y=number_of_outliers)) + 
    	geom_boxplot() + 
    	labs(x="Outlier p-value threshold",y="Number of outlier clusters per tissue", title=rna_type) + 
    	gtex_v8_figure_theme() +
    	theme(plot.title = element_text(hjust = 0.5))

    return(plotter)

}


box_plot_showing_fraction_of_clusters_with_outlier_per_tissue <- function(tissues, outlier_calls_dir) {
	tissue_name_arr <- c()
	frac_outlier_arr <- c()
	pvalue_arr <- c()
	rna_type_arr <- c()

	thresholds <- c(.0001, .00001, .000001)

	rna_type <- "lincrna"
	for (tissue_iter in 1:length(tissues)) {
		tissue_name <- tissues[tissue_iter]
		outlier_file <- paste0(outlier_calls_dir, tissue_name, "_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_", rna_type, ".txt")
		outlier_data <- read.table(outlier_file, header=TRUE)
		outlier_data <- outlier_data[,2:(dim(outlier_data)[2])]

		for (threshold_iter in 1:length(thresholds)) {
			threshold <- thresholds[threshold_iter]
			num_outliers <- sum(rowSums(outlier_data < threshold) != 0.0)

			tissue_name_arr <- c(tissue_name_arr, tissue_name)
			#num_outlier_arr <- c(num_outlier_arr, num_outliers)
			frac_outlier_arr <- c(frac_outlier_arr, num_outliers/dim(outlier_data)[1])
			pvalue_arr <- c(pvalue_arr, threshold)
			rna_type_arr <- c(rna_type_arr, rna_type)
		}
	}
	rna_type <- "proteincoding"
	for (tissue_iter in 1:length(tissues)) {
		tissue_name <- tissues[tissue_iter]
		outlier_file <- paste0(outlier_calls_dir, tissue_name, "_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_", rna_type, ".txt")
		outlier_data <- read.table(outlier_file, header=TRUE)
		outlier_data <- outlier_data[,2:(dim(outlier_data)[2])]

		for (threshold_iter in 1:length(thresholds)) {
			threshold <- thresholds[threshold_iter]
			num_outliers <- sum(rowSums(outlier_data < threshold) != 0.0)

			tissue_name_arr <- c(tissue_name_arr, tissue_name)
			#num_outlier_arr <- c(num_outlier_arr, num_outliers)
			frac_outlier_arr <- c(frac_outlier_arr, num_outliers/dim(outlier_data)[1])
			pvalue_arr <- c(pvalue_arr, threshold)
			rna_type_arr <- c(rna_type_arr, rna_type)
		}
	}


	df <- data.frame(rna_type=factor(rna_type_arr), tissue=tissue_name_arr, fraction_outliers=frac_outlier_arr, pvalue_threshold=factor(pvalue_arr, levels=c("1e-04", "1e-05", "1e-06")))
	
	plotter <- ggplot(df,aes(x=pvalue_threshold, y=fraction_outliers, fill=rna_type)) + 
    	geom_boxplot() + 
    	labs(x="Outlier p-value threshold",y="Fraction of clusters with an outlier", fill="") + 
    	gtex_v8_figure_theme() +
    	theme(plot.title = element_text(hjust = 0.5))

    return(plotter)

}

get_color_vector <- function(tissue_colors, tissue_names) {
	colors <- c()
	for (iter in 1:length(tissue_names)) {
		tissue_name <- as.character(tissue_names[iter])
		if (tissue_name == "Brain_Spinal_cord_cervical_c-1") {
			tissue_name = "Brain_Spinal_cord_cervical_c1"
		}
		if (tissue_name == "Cells_EBV-transformed_lymphocytes") {
			tissue_name = "Cells_EBVtransformed_lymphocytes"
		}

		indices = tissue_name==tissue_colors$tissue_id

		color <- tissue_colors$tissue_color_hex[indices]
		
		colors <- c(colors, paste0('#',color))
	}
	return(colors)
}


tbt_variant_outlier_enrichment_errorbar_plot <- function(tissue_by_tissue_enrichment_file, color_vector, gene_type) {
	enrichments <- read.table(tissue_by_tissue_enrichment_file, header=FALSE)
	num_tissues <- dim(enrichments)[1]

	# Initialize vectors
	tissue_names <- c()
	odds_ratios <- c()
	lower_bounds <- c()
	upper_bounds <- c()

	# Loop through tissues
	for (tissue_number in 1:num_tissues) {
		# Extract tissue name for this line
		tissue_name <- as.character(enrichments[tissue_number, 1])
		# Compute odds ratios for this line
		a <- enrichments[tissue_number, 2]
		b <- enrichments[tissue_number, 3]
		c <- enrichments[tissue_number, 4]
		d <- enrichments[tissue_number, 5]
		if (a>0.0) {
			if (c==0.0) {
				c = c+1
			}
			orat <- (a/b)/(c/d)
			# Compute error bars for this orat
			log_orat <- log(orat)
			log_bounds <- 1.96*sqrt((1.0/a) - (1.0/b) + (1.0/c) - (1.0/d))
			#upper_bound <- exp(log_orat + log_bounds)
			#lower_bound <- exp(log_orat - log_bounds) 
			upper_bound <- orat*exp(log_bounds)
			lower_bound <- orat*exp(-log_bounds)

			# Add information to vectors
			tissue_names <- c(tissue_names, tissue_name)
			odds_ratios <- c(odds_ratios, orat)
			lower_bounds <- c(lower_bounds, lower_bound)
			if (upper_bound > 400) {
				upper_bound = 400
			}
			upper_bounds <- c(upper_bounds, upper_bound)
		} else {
			tissue_names <- c(tissue_names, tissue_name)
			odds_ratios <- c(odds_ratios, NA)
			lower_bounds <- c(lower_bounds, NA)
			upper_bounds <- c(upper_bounds, NA)
		}
	}

	tissue_names <- gsub("_", " ", tissue_names)
	# Add information to data frame
	df <- data.frame(tissue_names=factor(tissue_names), odds_ratios=odds_ratios, lower_bounds=lower_bounds, upper_bounds=upper_bounds)
	print(df)
	error_bar_plot <- ggplot() + geom_errorbar(data=df, mapping=aes(x=tissue_names,ymin=lower_bounds, ymax=upper_bounds), width=0, colour=color_vector) +
					geom_point(data=df, mapping=aes(x=tissue_names, y=odds_ratios), colour=color_vector) +
					labs(x = "Tissue", y = "Relative risk", title=gene_type) +
					geom_hline(yintercept=1) + 
					gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5)) +
					ylim(0,400) +
					theme(plot.title = element_text(hjust = 0.5))

	return(error_bar_plot)
}

tissue_names_file <- args[1]
outlier_calls_dir <- args[2]
variant_enrichment_dir <- args[3]
visualization_dir <- args[4]
tissue_colors_file <- args[5]


########################
# Extract tissue color information
########################
# Extract vector tissue names
tissue_names <- as.character(unlist(read.table(tissue_names_file,header=FALSE), use.names=FALSE))
tissue_colors = read.table(tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")
# Get vector of hex colors in correct order
color_vector <- get_color_vector(tissue_colors, tissue_names)

tissues <- tissue_names


###################################################
# Make plot showing number of clusters per tissue
###################################################
output_file <- paste0(visualization_dir, "number_of_lincrna_clusters_per_tissue.pdf")
lincrna_bar_plot_num_clusters_per_tissue <- bar_plot_showing_number_of_clusters_per_tissue(tissues, outlier_calls_dir, "lincrna")
ggsave(lincrna_bar_plot_num_clusters_per_tissue, file=output_file,width=7.2, height=6,units="in")

output_file <- paste0(visualization_dir, "number_of_protein_coding_clusters_per_tissue.pdf")
lincrna_bar_plot_num_clusters_per_tissue <- bar_plot_showing_number_of_clusters_per_tissue(tissues, outlier_calls_dir, "proteincoding")
ggsave(lincrna_bar_plot_num_clusters_per_tissue, file=output_file,width=7.2, height=6,units="in")




###################################################
# Make boxplot showing number outlier clusters per tissue as a function of p-value threshold
###################################################
output_file <- paste0(visualization_dir, "number_of_lincrna_outlier_clusters_per_tissue.pdf")
lincrna_box_plot_num_outlier_clusters_per_tissue <- box_plot_showing_number_of_outlier_clusters_per_tissue(tissues, outlier_calls_dir, "lincrna")
ggsave(lincrna_box_plot_num_outlier_clusters_per_tissue, file=output_file,width=7.2, height=5,units="in")

output_file <- paste0(visualization_dir, "number_of_proteincoding_outlier_clusters_per_tissue.pdf")
proteincoding_box_plot_num_outlier_clusters_per_tissue <- box_plot_showing_number_of_outlier_clusters_per_tissue(tissues, outlier_calls_dir, "proteincoding")
ggsave(proteincoding_box_plot_num_outlier_clusters_per_tissue, file=output_file,width=7.2, height=5,units="in")



###################################################
# Make boxplot showing fraction of tested clusters with an outlier per tissue
###################################################
output_file <- paste0(visualization_dir, "fraction_of_clusters_that_have_an_outlier_per_tissue.pdf")
fraction_of_clusters_with_outlier_boxplot <- box_plot_showing_fraction_of_clusters_with_outlier_per_tissue(tissues, outlier_calls_dir)
ggsave(fraction_of_clusters_with_outlier_boxplot, file=output_file, width=7.2, height=5, units="in")



#######################
# Make errorbar plot showing enrichment of rare variants nearby splice sites within splicing outliers
#########################
pvalue_threshold = '.00001'
distance = '6'
gene_type <- "lincrna"
tissue_by_tissue_enrichment_file <- paste0(variant_enrichment_dir, "tbt_variant_outlier_enrichment_pvalue_", pvalue_threshold, "_distance_", distance, "_version_all_gene_type_", gene_type,".txt")
output_file <- paste0(visualization_dir, "tbt_variant_outlier_enrichment_pvalue_", pvalue_threshold, "_distance_", distance, "_errorbar_gene_type_", gene_type, ".pdf")
tbt_variant_outlier_errorbar_plot <- tbt_variant_outlier_enrichment_errorbar_plot(tissue_by_tissue_enrichment_file, color_vector, gene_type)
ggsave(tbt_variant_outlier_errorbar_plot, file=output_file, width=7.2, height=5, units="in")

#######################
# Make errorbar plot showing enrichment of rare variants nearby splice sites within splicing outliers
#########################
pvalue_threshold = '.00001'
distance = '6'
gene_type <- "proteincoding"
tissue_by_tissue_enrichment_file <- paste0(variant_enrichment_dir, "tbt_variant_outlier_enrichment_pvalue_", pvalue_threshold, "_distance_", distance, "_version_all_gene_type_", gene_type,".txt")
output_file <- paste0(visualization_dir, "tbt_variant_outlier_enrichment_pvalue_", pvalue_threshold, "_distance_", distance, "_errorbar_gene_type_", gene_type, ".pdf")
tbt_variant_outlier_errorbar_plot <- tbt_variant_outlier_enrichment_errorbar_plot(tissue_by_tissue_enrichment_file, color_vector, gene_type)
ggsave(tbt_variant_outlier_errorbar_plot, file=output_file, width=7.2, height=5, units="in")

