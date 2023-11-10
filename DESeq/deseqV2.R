#!/usr/bin/env Rscript

#BiocManager::install("DESeq2")
#BiocManager::install("tximport")
#install.packages("argparse")

library("DESeq2")
library("tximport")
library("ggplot2")
#library("patchwork")
library("argparse")
library("ggrepel") #for showing labels

parser = ArgumentParser()

parser$add_argument("-st", "--star_matrix", type="character", default="./raw_counts_matrix.txt", help="Path to STAR raw counts matrix file")
parser$add_argument("-sl", "--salmon_folder", type="character", default="./salmon", help="Path to salmon count folders")
parser$add_argument("-t2g", "--tx2gene", type="character", default="./tx2gene.gencode.v29.csv", help="path to mock theoretical composition (%) in csv (genus,composition)")
parser$add_argument("-m", "--metadata", type="character", default="./meta_sorted.tsv", help="Path to metadata file for DESeq2")
parser$add_argument("-t", "--data_type", type="character", default="STAR", help="Run the script on STAR or salmon output [STAR,salmon]")
parser$add_argument("-g", "--gene_counts", type="character", help="Output gene counts for specific gene")

# Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

star_counts = args$star_matrix
metadata = args$metadata
tx2gene = args$tx2gene
sal_counts = args$salmon_folder
data_type = args$data_type
gene = args$gene_counts


#star_counts = "raw_counts_matrix.txt" #star raw counts location
#metadata = "meta_sorted.tsv" #metadata for DESeq2
#tx2gene = "tx2gene.gencode.v29.csv" #tx2 gene file for salmon
#sal_counts = "salmon" #base salmon counts folder

# Gene name should be accroding to gene name in tx2gene file, not accession!! Review tx2gene file if necessary

# FUNCTIONS

# Read tables

read_table = function(x, y=NULL, z=TRUE) {
  read.table(x, sep = "\t", header = z, row.names = y, check.names=FALSE)
}

# Plot gene counts

plot_gene = function(fitted, gene, group, title) { # group not as string
  gene = unique(tx2gene[tx2gene$V3==gene,]$V2) # get gene id from tx2gene file
  plot_data = plotCounts(fitted, gene = gene, 
                         intgroup = group, 
                         returnData=TRUE) # return dataframe for plotting
  plot_gene = ggplot(plot_data, 
                     aes(.data[[group]], .data[["count"]], label=rownames(plot_data))) + # .data[[]] is very useful for function
    geom_point(position=position_jitter(w=0.1,h=0)) +
    geom_label_repel() +
    labs(title = title)
  plot_gene
}

# Plot PCA

plot_PCA = function(fitted, title, group) {
  #group_str = deparse(substitute(group))
  vsd = vst(fitted)
  plot_data = plotPCA(object = vsd,
                      intgroup = group, returnData = TRUE)
  percentVar = round(100 * attr(plot_data, "percentVar"))
  plot_pca = ggplot(plot_data, aes(PC1, PC2, color=.data[[group]], label=rownames(plot_data))) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    geom_label_repel() + 
    labs(title = title) +
    coord_fixed()
  plot_pca
}

# ANALISYS

# Read tables

metadata = read_table(metadata)

tx2gene = read_table(tx2gene, NULL, FALSE) #no header

# DESeq2 requires this

rownames(metadata) <- metadata$SampleName

### STAR

if (data_type == "STAR") {
  
  print(data_type)

  star_counts = read_table(star_counts, 1) # 1 for rownames
  
  se_star_matrix <- DESeqDataSetFromMatrix(countData = star_counts,
                                         colData = metadata,
                                         design = ~ Genotype)

  # Check number of rows

  print(nrow(se_star_matrix))

  # Filter counts (10 is standard values)

  se_star_matrix = se_star_matrix[rowSums(counts(se_star_matrix)) > 10, ]

  # Check number of rows again: less genes

  print(nrow(se_star_matrix))

  # Differential expression analysis

  se_star_matrix2 = DESeq(se_star_matrix)

  # ask about + 1

  # Normalize counts

  norm_counts <- log2(counts(se_star_matrix2, normalized = TRUE)+1)

  # add the gene symbols
  norm_counts_symbols = merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(norm_counts), norm_counts), by=1, all=F)

  # write normalized counts to text file
  write.table(norm_counts_symbols, "normalized_counts_star.txt", quote=F, col.names=T, row.names=F, sep="\t")

  print(resultsNames(se_star_matrix2))

  head(se_star_matrix2)

  #Select gene from interest according to tx2gene table and plot

  if (!is.null(gene) ) {
  
    png("STAR_gene.png")
    print(plot_gene(se_star_matrix2, gene, "Genotype", "STAR - Bmal1"))
    dev.off()
  
  }

  png("STAR_PCA.png", width = 1000, height = 1000)
  print(plot_PCA(se_star_matrix2, "STAR", "Genotype"))
  dev.off()

  ### Salmon
  
} else if (data_type == "salmon") {

  files <- dir(sal_counts, recursive=TRUE, pattern="quant.sf", full.names=TRUE)

  # Add file names from folder names

  names(files) <- dir(sal_counts)
  
  # Import counts

  txi <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene)

  se_salmon_matrix <- DESeqDataSetFromTximport(txi,
                                             colData = metadata,
                                             design = ~ Genotype)

  se_salmon_matrix = se_salmon_matrix[rowSums(counts(se_salmon_matrix)) > 10, ]

  se_salmon_matrix2 = DESeq(se_salmon_matrix)

  #ask about + 1

  norm_counts = log2(counts(se_salmon_matrix2, normalized = TRUE)+1)

  # add the gene symbols
  norm_counts_symbols = merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(norm_counts), norm_counts), by=1, all=F)

  # write normalized counts to text file
  write.table(norm_counts_symbols, "normalized_counts_salmon.txt", quote=F, col.names=T, row.names=F, sep="\t")

  #View results
  #res_salmon = results(se_salmon_matrix2)

  #head(res_salmon)

  if (!is.null(gene) ) {
  
    png("salmon_gene.png")
    print(plot_gene(se_salmon_matrix2, gene, "Genotype", "Salmon - Bmal1"))
    dev.off()
  
  }

  png("salmon_PCA.png", width = 1000, height = 1000)
  print(plot_PCA(se_salmon_matrix2, "STAR", "Genotype"))
  dev.off()

}