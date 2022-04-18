library(cluster)
library(pheatmap)
library(EnhancedVolcano)
library(tidyverse)
library(Biostrings)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# DE Functions
write_files <- function(results, numerator, denominator){
  # these are all the genes that are differentially expressed between the two conditions, not just the significant ones
  no_outliers <- results[!is.na(results$padj),]
  write.csv(no_outliers, paste0(output_dir,numerator,"_",denominator,"_all.csv"), row.names = TRUE, col.names = TRUE)
  
  # these are the genes that are significantly differentially expressed by FDR 10% and abs(log2fc) > log2(1.5)
  sig_padj_genes <- no_outliers[no_outliers$padj < 0.1,]
  sig_padj_fc_genes <- sig_padj_genes[abs(sig_padj_genes$log2FoldChange) > lfc.threshold,]
  write.csv(sig_padj_fc_genes, paste0(output_dir,numerator,"_",denominator,"_significant.csv"), row.names = TRUE, col.names = TRUE)  
}

generate_de_section <- function(dds_obj, numerator_condition, denominator_condition){
  res1 <- results(dds_obj, contrast = c( "condition", numerator_condition, denominator_condition), test="Wald")
  res <- lfcShrink(dds = dds_obj, res = res1, type = "normal", contrast = c( "condition", numerator_condition, denominator_condition))
  res$gene_biotype <- ensembl.genes$gene_biotype[match(row.names(res1), ensembl.genes$gene_id)]
  res$external_gene_name <- ensembl.genes$external_gene_name[match(row.names(res1), ensembl.genes$gene_id)]
  print(paste("Number of significant genes (padj < 0.1 & log2FoldChange < log2(1.5)):", sum(res$padj < 0.1 & abs(res$log2FoldChange) > lfc.threshold, na.rm = T)))
  hist(res$pvalue)
  ### Writing out .csv files
  write_files(res, numerator_condition, denominator_condition)
  EnhancedVolcano(res, lab=res$external_gene_name, pCutoff = 0.05, FCcutoff = lfc.threshold,
                  x = "log2FoldChange", y = "padj",
                  title = NULL,
                  subtitle = NULL,
                  legendPosition = 'right',
                  legendLabSize = 8)
}

## This function generates the zscore table ordered by cluster number
## It takes in the full zscore matrix and the number of clusters desired

generate_data <- function(zscores, n_clust, fun = c("kmeans","pam"), seed = 2){
  set.seed(seed)
  if(fun == "kmeans") {
    km = kmeans(zscores, n_clust) ## Running kmeans will return 9 components and var[1] will give back the first component which is "cluster", var[2] will return "centers" etc.
    clust <- cbind(zscores, km$cluster)
  } else {
    dd = as.dist((1 - cor(t(zscores)))/2)
    pam_clust = pam(dd, n_clust, diss = TRUE)
    clust <- cbind(zscores, pam_clust$clustering)
  }
  n_col = dim(clust)[2]
  colnames(clust)[n_col] = "Cluster"
  ordered <- order(clust[,n_col])
  clust <- clust[ordered, ]
  return(clust)
}


## This function generates the heatmap using the results from the above function
## It takes in the data frame, annotation data, and annotation colours

generate_heatmap <- function(clust_data, annotation, annotation_col){
  return(pheatmap(clust_data[,1:(ncol(clust_data)-1)],
                  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                  fontsize_row = 5.5,
                  annotation_col = annotation,
                  annotation_colors = anno_colours,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE))
}


## This function plots the barplots using results from EnrichR 
## It takes in an enrichR result, name of the plot, and the number of terms to show (default being 20)

plot_enrichr <- function(data.frame, name, showCategory = 20){
  plot = data.frame %>%
    mutate(Term = gsub("\\([^()]*\\)", "", Term),
           Term = factor(Term, levels = rev(Term))) %>%
    mutate(Annotated = as.numeric(str_extract(as.character(Overlap), "\\d+$")),
           Significant = as.numeric(str_extract(as.character(Overlap), "^\\d+")),
           Ratio = Significant/Annotated) %>%
    arrange(Adjusted.P.value) %>%
    head(showCategory) %>%
    ggplot(aes(x = reorder(Term,dplyr::desc(Adjusted.P.value)), y = Ratio, fill = Adjusted.P.value)) +
    geom_bar(stat = "identity") +
    ggpubr::rotate() +
    xlab(NULL) + 
    ylab("Gene Ratio") +
    scale_fill_continuous(low = "red", high = "blue") + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    labs(title = name,
         fill = "Adjusted p-value") +
    guides(fill = guide_colorbar(title = "Adjusted p-value", reverse = TRUE)) +
    theme(text = element_text(size=15), legend.position="bottom", legend.direction = "vertical")
  return(plot)
}


get_wanted_terms <- function(data.frame, output_dir, name, wanted_terms){
  df <- data.frame %>%
    mutate(Annotated = as.numeric(str_extract(as.character(Overlap), "\\d+$")),
           Significant = as.numeric(str_extract(as.character(Overlap), "^\\d+")),
           Ratio = Significant/Annotated) %>%
    dplyr::filter(Term %in% wanted_terms) %>%
    arrange(Adjusted.P.value)
    if(nrow(df) != 0) {
     saveRDS(df, paste0(output_dir,name,".rds"))
    }
}

## This function takes the cluster of interest and returns the fasta file in the stated output directory
generate_fasta <- function(output_dir, sig_de_granges, data_frame, cluster_name){
  granges <- sig_de_granges[sig_de_granges$gene_id %in% rownames(data_frame),]
  fasta <- getSeq(Dmelanogaster, granges)
  writeXStringSet(fasta, paste(output_dir,cluster_name,".fa", sep = ""))
}


## This function takes the cluster of interest and returns the fasta file of EVERYTHING BUT THE CLUSTER OF INTEREST in the stated output directory (differentially expressed genes)
extract_control_deg <- function(output_dir, sig_de_granges, cluster_of_interest, cluster_name){
  deg_control <- sig_de_granges[-which(sig_de_granges$gene_id %in% rownames(cluster_of_interest)),] #note the order
  deg_fasta <- getSeq(Dmelanogaster, deg_control)
  writeXStringSet(deg_fasta, paste(output_dir,cluster_name,".fa", sep = ""))
}


## This function takes the cluster of interest and returns the fasta file of EVERYTHING BUT THE CLUSTER OF INTEREST in the stated output directory (all expressed genes)
extract_control_all <- function(output_dir, all_expressed_granges, cluster_of_interest, cluster_name){
  all_control <- all_expressed_granges[-which(all_expressed_granges$gene_id %in% rownames(cluster_of_interest)),] #note the order
  all_fasta <- getSeq(Dmelanogaster, all_control)
  writeXStringSet(all_fasta, paste(output_dir,cluster_name,".fa", sep = ""))
}

## plots box plots for each cluster

cluster_boxplot <- function(cluster_df, name, control, sen_oksm, oksm){
  return(cluster_df %>%
           dplyr::select(-gene_name) %>%
           gather() %>%
           mutate(condition = ifelse(grepl(control, key), "Control", 
                                     ifelse(grepl(sen_oksm, key), "Senolytic_OKSM",
                                            ifelse(grepl(oksm, key, fixed = TRUE), "OKSM", "Senolytic")))) %>%
           mutate(condition = factor(condition, levels = c("Control", "OKSM", "Senolytic", "Senolytic_OKSM"))) %>%
           ggplot(aes(x = condition, y = value, fill = condition)) + 
           geom_boxplot() + 
           theme_classic() + 
           scale_fill_manual(values = c("#BF3100", "#E9B44C", "#1B998B", "#5D576B")) +
           theme(axis.text.x = element_text(angle = 90, colour="black", hjust = 1)) + 
           theme(axis.text.x = element_blank(), legend.position="bottom") +
           xlab(NULL) +
           ylab(NULL) +
           labs(title = name) #+
         #scale_y_continuous(limits = c(lower_y_limit, upper_y_limit))
  )
}
