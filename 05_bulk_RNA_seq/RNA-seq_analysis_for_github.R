load libraries
```{r}
library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(edgeR)
library(tidyverse)
```


functions for plotting
```{r}
plot_PCA <- function(vsd.obj, ntop = 500) {
  pcaData <- plotPCA(vsd.obj, intgroup = c("condition"), returnData = TRUE, ntop = ntop)
  percentVar <- round(100 * attr(pcaData, "percentVar"))


  # Calculate the limits of the gray area (plot panel)
  x_limit <- max(abs(pcaData$PC1))+1
  y_limit <- max(abs(pcaData$PC2))+1
  

  ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_rect(xmin = -x_limit, xmax = x_limit, ymin = -y_limit, ymax = y_limit, fill = "transparent", color = "black") +
    geom_point(size = 3) +
    labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance")
    ) +
    ggrepel::geom_text_repel(aes(label = name), color = "black") +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'))
}

generate_DESeq_object <- function (my_data, groups) {
  data_subset1 <- my_data[,grep(str_c("^", groups[1]), colnames(my_data))]
  data_subset2 <- my_data[,grep(str_c("^", groups[2]), colnames(my_data))]
  my_countData <- cbind(data_subset1, data_subset2)
  condition <- c(rep(groups[1],ncol(data_subset1)), rep(groups[2],ncol(data_subset2)))
  my_colData <- as.data.frame(condition)
  rownames(my_colData) <- colnames(my_countData)
  print(my_colData)
  dds <- DESeqDataSetFromMatrix(countData = my_countData,
                              colData = my_colData,
                              design = ~ condition)
  dds <- DESeq(dds, quiet = T)
  return(dds)
  

}

```
```{r}
generate_DE_results <- function (dds, comparisons, padjcutoff = 0.05, log2cutoff = 0.5, cpmcutoff = 2) {
  # generate average counts per million metric from raw count data 
  raw_counts <- counts(dds, normalized = F)
  cpms <- enframe(rowMeans(edgeR::cpm(raw_counts)))
  colnames(cpms) <- c("ensembl_id", "avg_cpm")
  
  # extract DESeq results between the comparisons indicated
  res <- results(dds, contrast = c("condition", comparisons[1], comparisons[2]))[,-c(3,4)]
  
  # annotate the data with gene name and average counts per million value
  res <- as_tibble(res, rownames = "ensembl_id")
  # read in the annotation and append it to the data
  my_annotation <- read.csv("C:/projects/genome_sequence/GRCm38.vM10.annotation.csv", header = T, stringsAsFactors = F,sep='\t')
  res <- left_join(res, my_annotation, by = c("ensembl_id" = "gene_id"))
  # append the average cpm value to the results data
  res <- left_join(res, cpms, by = c("ensembl_id" = "ensembl_id"))
  
  # combine normalized counts with entire DE list
  normalized_counts <- round(counts(dds, normalized = TRUE),3)
  pattern <- str_c(comparisons[1], "|", comparisons[2])
  combined_data <- as_tibble(cbind(res, normalized_counts[,grep(pattern, colnames(normalized_counts))] ))
  combined_data <- combined_data[order(combined_data$log2FoldChange, decreasing = T),]
  
  # make ordered rank file for GSEA, selecting only protein coding genes
  res_prot <- res[which(res$Gene.type == "protein_coding"),]
  res_prot_ranked <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("Gene.name", "log2FoldChange")]
  res_prot_ranked <- na.omit(res_prot_ranked)
  res_prot_ranked$Gene.name <- str_to_upper(res_prot_ranked$Gene.name)
  
  # generate sorted lists with the indicated cutoff values
  res <- res[order(res$log2FoldChange, decreasing=TRUE ),]
  de_genes_padj <- res[which(res$padj < padjcutoff),]
  de_genes_log2f <- res[which(abs(res$log2FoldChange) > log2cutoff & res$padj < padjcutoff),]
  de_genes_cpm <- res[which(res$avg_cpm > cpmcutoff & res$padj < padjcutoff),]
  
  # write output to files
  write.csv (de_genes_padj, file = paste0(comparisons[1], "_vs_", comparisons[2], "_padj_cutoff.csv"), row.names =F)
  write.csv (de_genes_log2f, file = paste0(comparisons[1], "_vs_", comparisons[2], "_log2f_cutoff.csv"), row.names =F)
  write.csv (de_genes_cpm, file = paste0(comparisons[1], "_vs_", comparisons[2], "_cpm_cutoff.csv"), row.names =F)
  write.csv (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_allgenes.csv"), row.names =F)
  write.table (res_prot_ranked, file = paste0(comparisons[1], "_vs_", comparisons[2], "_rank.rnk"), sep = "\t", row.names = F, quote = F)
  
  writeLines( paste0("For the comparison: ", comparisons[1], "_vs_", comparisons[2], ", out of ", nrow(combined_data), " genes, there were: \n", 
               nrow(de_genes_padj), " genes below padj ", padjcutoff, "\n",
               nrow(de_genes_log2f), " genes below padj ", padjcutoff, " and above a log2FoldChange of ", log2cutoff, "\n",
               nrow(de_genes_cpm), " genes below padj ", padjcutoff, " and above an avg cpm of ", cpmcutoff, "\n",
               "Gene lists ordered by log2fchange with the cutoffs above have been generated.") )
  gene_count <- tibble (cutoff_parameter = c("padj", "log2fc", "avg_cpm" ), 
                        cutoff_value = c(padjcutoff, log2cutoff, cpmcutoff), 
                        signif_genes = c(nrow(de_genes_padj), nrow(de_genes_log2f), nrow(de_genes_cpm)))
  invisible(gene_count)
}
```
```{r}
### when you want to highlight genes from another list in the volcano plot,here the list is from TBXT binding sites
plot_volcano3 <- function(res2, padj_cutoff, nlabel = 10, label.by = "padj") {
  # assign significance to results based on padj
  res2 <- mutate(res2, significance = ifelse(res2$padj < padj_cutoff, paste0("padj < ", padj_cutoff), paste0("padj > ", padj_cutoff)))
  res2 <- res2[!is.na(res2$significance), ]
  significant_genes <- res2 %>% filter(significance == paste0("padj < ", padj_cutoff))
  
  # get labels for the highest or lowest genes according to either padj or log2FoldChange

  res2$InGeneList <- ifelse(trimws(res2$GeneName) %in% trimws(pnas$gene), "In List", "Not In List")
  
  if (nrow(significant_genes) == 0) {
    stop("No significant genes found for the given padj_cutoff.")
  }
  
  if (label.by == "padj") { 
    top_genes <- res2[trimws(res2$GeneName) %in% trimws(pnas$gene), ] %>% arrange(padj) %>% head(nlabel)
    bottom_genes <- significant_genes %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% head(nlabel)
  } else if (label.by == "log2FoldChange") {
    top_genes <- head(arrange(significant_genes, desc(log2FoldChange)), nlabel)
    bottom_genes <- head(arrange(significant_genes, log2FoldChange), nlabel)
  } else {
    stop("Invalid label.by argument. Choose either padj or log2FoldChange.")
  }

  ggplot(res2, aes(x = log2FoldChange, y = -log10(padj), color = InGeneList)) + 
    geom_point(alpha = ifelse(res2$InGeneList == "In List", 4, 0.5), size = 1,shape=16) +
    theme_bw(base_size = 12) + 
    xlab("Log2(Fold change)") + 
    ylab("-Log10(P.adj)") + 
    theme(plot.title = element_text(size = 15, hjust = 0.5)) + 
    scale_colour_manual(values = c('In List' = 'red', 'Not In List' = 'gray')) + 
    geom_hline(yintercept = -log10(0.05), lty = 4,color = "steelblue") +
    geom_vline(xintercept = c(-0.5, 0.5), lty = 4,color = "steelblue") +
    ggrepel::geom_text_repel(data = top_genes,
                             aes(label = GeneName, color = Change), size = 5) +
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'))
}

```

```{r}

generate_DE_results <- function (dds, comparisons, padjcutoff = 0.05, log2cutoff = 0.5, cpmcutoff = 2) {
  # generate average counts per million metric from raw count data 
  raw_counts <- counts(dds, normalized = F)
  cpms <- enframe(rowMeans(edgeR::cpm(raw_counts)))
  colnames(cpms) <- c("ensembl_id", "avg_cpm")
  
  # extract DESeq results between the comparisons indicated
  res <- results(dds, contrast = c("condition", comparisons[1], comparisons[2]))[,-c(3,4)]
  
  # annotate the data with gene name and average counts per million value
  res <- as_tibble(res, rownames = "ensembl_id")
  # read in the annotation and append it to the data
  my_annotation <- read.csv("C:/projects/genome_sequence/GRCm38.vM10.annotation.csv", header = T, stringsAsFactors = F,sep='\t')
  res <- left_join(res, my_annotation, by = c("ensembl_id" = "gene_id"))
  # append the average cpm value to the results data
  res <- left_join(res, cpms, by = c("ensembl_id" = "ensembl_id"))
  
  # combine normalized counts with entire DE list
  normalized_counts <- round(counts(dds, normalized = TRUE),3)
  pattern <- str_c(comparisons[1], "|", comparisons[2])
  combined_data <- as_tibble(cbind(res, normalized_counts[,grep(pattern, colnames(normalized_counts))] ))
  combined_data <- combined_data[order(combined_data$log2FoldChange, decreasing = T),]
  
  # make ordered rank file for GSEA, selecting only protein coding genes
  res_prot <- res[which(res$X == "protein_coding"),]
  res_prot_ranked <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("GeneName", "log2FoldChange")]
  res_prot_ranked <- na.omit(res_prot_ranked)
  res_prot_ranked$GeneName <- str_to_upper(res_prot_ranked$GeneName)
  
  # generate sorted lists with the indicated cutoff values
  res <- res[order(res$log2FoldChange, decreasing=TRUE ),]
  de_genes_padj <- res[which(res$padj < padjcutoff),]
  de_genes_log2f <- res[which(abs(res$log2FoldChange) > log2cutoff & res$padj < padjcutoff),]
  de_genes_cpm <- res[which(res$avg_cpm > cpmcutoff & res$padj < padjcutoff),]
  
  # write output to files
  write.csv (de_genes_padj, file = paste0(comparisons[1], "_vs_", comparisons[2], "_padj_cutoff.csv"), row.names =F)
  write.csv (de_genes_log2f, file = paste0(comparisons[1], "_vs_", comparisons[2], "_log2f_cutoff.csv"), row.names =F)
  write.csv (de_genes_cpm, file = paste0(comparisons[1], "_vs_", comparisons[2], "_cpm_cutoff.csv"), row.names =F)
  write.csv (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_allgenes.csv"), row.names =F)
  write.table (res_prot_ranked, file = paste0(comparisons[1], "_vs_", comparisons[2], "_rank.rnk"), sep = "\t", row.names = F, quote = F)
  
  writeLines( paste0("For the comparison: ", comparisons[1], "_vs_", comparisons[2], ", out of ", nrow(combined_data), " genes, there were: \n", 
               nrow(de_genes_padj), " genes below padj ", padjcutoff, "\n",
               nrow(de_genes_log2f), " genes below padj ", padjcutoff, " and above a log2FoldChange of ", log2cutoff, "\n",
               nrow(de_genes_cpm), " genes below padj ", padjcutoff, " and above an avg cpm of ", cpmcutoff, "\n",
               "Gene lists ordered by log2fchange with the cutoffs above have been generated.") )
  gene_count <- tibble (cutoff_parameter = c("padj", "log2fc", "avg_cpm" ), 
                        cutoff_value = c(padjcutoff, log2cutoff, cpmcutoff), 
                        signif_genes = c(nrow(de_genes_padj), nrow(de_genes_log2f), nrow(de_genes_cpm)))
  invisible(gene_count)
}
```




```{r}
variable_gene_heatmap <- function (vsd.obj, num_genes = 500, annotation, title = "") {
  brewer_palette <- "RdBu"
  # Ramp the color in order to get the scale.
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  # get the stabilized counts from the vsd object
  stabilized_counts <- assay(vsd.obj)
  # calculate the variances by row(gene) to find out which genes are the most variable across the samples.
  row_variances <- rowVars(stabilized_counts)
  # get the top most variable genes
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=T)[1:num_genes],]
  # subtract out the means from each row, leaving the variances for each gene
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=T)
  # replace the ensembl ids with the gene names
  gene_names <- annotation$Gene.name[match(rownames(top_variable_genes), annotation$Gene.stable.ID)]
  rownames(top_variable_genes) <- gene_names
  # reconstruct colData without sizeFactors for heatmap labeling
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  # draw heatmap using pheatmap
  pheatmap::pheatmap(top_variable_genes, color = mr, annotation_col = coldata, fontsize_col = 8, fontsize_row = 250/num_genes, border_color = NA, main = title)
}
```

```{r}
plot_volcano2 <- function(res2, padj_cutoff, nlabel = 10, label.by = "padj") {
  # assign significance to results based on padj
  res2 <- mutate(res2, significance = ifelse(res2$padj < padj_cutoff, paste0("padj < ", padj_cutoff), paste0("padj > ", padj_cutoff)))
  res2 <- res2[!is.na(res2$significance), ]
  significant_genes <- res2 %>% filter(significance == paste0("padj < ", padj_cutoff))
  
  # get labels for the highest or lowest genes according to either padj or log2FoldChange


  
  if (nrow(significant_genes) == 0) {
    stop("No significant genes found for the given padj_cutoff.")
  }
  
  if (label.by == "padj") { 
    top_genes <- significant_genes %>% arrange(padj) %>% head(nlabel)
    bottom_genes <- significant_genes %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% head(nlabel)
  } else if (label.by == "log2FoldChange") {
    top_genes <- head(arrange(significant_genes, desc(log2FoldChange)), nlabel)
    bottom_genes <- head(arrange(significant_genes, log2FoldChange), nlabel)
  } else {
    stop("Invalid label.by argument. Choose either padj or log2FoldChange.")
  }

ggplot(res2, aes(x=log2FoldChange, y=-log10(padj), color=Change)) + 
    geom_point(alpha=0.4, size=2) +  
    theme_bw(base_size = 12) + 
    xlab("Log2(Fold change)") + 
    ylab("-Log10(P.adj)") + 
    theme(plot.title = element_text(size=15,hjust = 0.5)) + 
    scale_colour_manual(values = c('UP' = 'red', 'DOWN' = 'steelblue')) + 
    geom_hline(yintercept = -log10(0.05), lty = 4) +
    geom_vline(xintercept = c(-0.5, 0.5), lty = 4) +
    ggrepel::geom_text_repel(data=top_genes, aes(label=head(GeneName,nlabel)), size = 3) +
    
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'))
}
```

read data
```{r}
data <- read.csv("C:/projects/TAIL/mouse/mouse_counts_all2.csv", header = T, row.names = "ensembl_id")
data <- data[,sort(colnames(data))]
head(data)
```

Creating the DESeq2 object
```{r}
# identify biological replicates
condition <- c(rep("WT", 2), rep("delE6_homo", 2), rep("insASAY_homo", 2),rep("delE6_hete", 2))
```
```{r}
# assign replicates to each sample name to construct colData
my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(data)
my_colData

```
```{r}
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = my_colData,
                              design = ~condition)
```
```{r}
dds <- DESeq(dds)
```
```{r}
normalized_counts <- counts(dds, normalized = T)

head(normalized_counts)
```

```{r}
normalized_counts <- rownames_to_column(as.data.frame(normalized_counts), var = "ensembl_id")
annotated_data <- right_join(annotation, normalized_counts, by = c("gene_id" = "ensembl_id"))
head(annotated_data)
```
```{r}
vsd <- vst(dds, blind = TRUE)
```

plot PCA(extended fig9a)
```{r}
plot_PCA(vsd)
```

generate DESEQ project

```{r}
colnames(data) <- c("WT1", "WT2", "homo_DelE61", "homo_DelE62", "insASAY1", "insASAY2", "hetero_DelE61","hetero_DelE62")
```


hete_delE6 vs WT differential gene analysis.
```{r}
sub_group <- c("hetero_DelE6","WT")
WT_delE6_hete <- generate_DESeq_object(data, sub_group)

```
```{r}
results(WT_delE6_hete, contrast = c("condition", "WT", "hete_DelE6"))
```
```{r}
WT_delE6_hete_vsd <- vst(WT_delE6_hete, blind = T)
variable_gene_heatmap(WT_delE6_hete_vsd, 1000, annotation = annotation, title = "WT_Y variable genes")

```

```{r}
my_annotation <- read.csv("C:/projects/genome_sequence/GRCm38.vM10.annotation.csv", header = T, stringsAsFactors = F,sep='\t')
```

```{r}
WT_delE6_hete_output <- generate_DE_results(WT_delE6_hete, c( "hetero_DelE6","WT"))
```

draw volcano plot
```{r}
### you should find the output differential gene file in your working directory with the suffix '*_allgenes.csv', and read it
res2 <- read.csv ('C:/projects/TAIL/hetero_DelE6_vs_WT_allgenes.csv', header = T)


```


```{r}
res2 <- res2 %>%
  mutate(
    Change = case_when(
      padj < 0.05 & abs(log2FoldChange) > 0.5 & log2FoldChange > 0.5 ~ "UP",
      padj < 0.05 & abs(log2FoldChange) > 0.5 & log2FoldChange <= 0.5 ~ "DOWN"
    )
  )

volcano <- plot_volcano2(res2, 0.05, nlabel = 30, label.by = "padj")
volcano

print(volcano)
dev.off()

```

and now we want to highlight the differential genes which are targeted by Tbxt from Lolas et al. (2014)
```{r}
pnas <- read.csv("./PNAS_genes.csv",header = T)
```

and attention: the red dots are the genes bound by Tbxt.so for those genes not differentially expressed in hete_delE6 vs WT are also highlighted in red, for those not differentially expressed are manually changed to grey color in Adobe Illustrator. The x-axis and label font are adjusted as well.

This replot the extended fig9b plot, for extended fig9c and extended fig9e plot you can easily change c( "hete_DelE6","WT") in line 352 to c( "insASAY","WT") or c( "homo_DelE6","WT"), and rerun the code from line 352 to line 401
```{r}
res2 <- res2 %>%
  mutate(
    Change = case_when(
      padj < 0.05 & abs(log2FoldChange) > 0.5 & log2FoldChange > 0.5 ~ "UP",
      padj < 0.05 & abs(log2FoldChange) > 0.5 & log2FoldChange <= 0.5 ~ "DOWN"
    )
  )

volcano <- plot_volcano3(res2, 0.05, nlabel = 30, label.by = "padj")
volcano

```
homo_delE6 vs WT differential gene analysis.
```{r}
sub_group <- c("homo_DelE6","WT")
WT_delE6_homo <- generate_DESeq_object(data, sub_group)
WT_delE6_homo_output <- generate_DE_results(WT_delE6_homo, sub_group)
### you should find the output differential gene file in your working directory with the suffix '*_allgenes.csv', and read it
res2 <- read.csv ('C:/projects/TAIL/homo_DelE6_vs_WT_allgenes.csv', header = T)


```
```{r}
res2 <- res2 %>%
  mutate(
    Change = case_when(
      padj < 0.05 & abs(log2FoldChange) > 0.5 & log2FoldChange > 0.5 ~ "UP",
      padj < 0.05 & abs(log2FoldChange) > 0.5 & log2FoldChange <= 0.5 ~ "DOWN"
    )
  )

volcano <- plot_volcano3(res2, 0.05, nlabel = 30, label.by = "padj")
volcano

```
insASAY VS WT differential gene analysis.
```{r}
sub_group <- c("insASAY","WT")
WT_insASAY <- generate_DESeq_object(data, sub_group)
WT_delE6_homo_output <- generate_DE_results(WT_insASAY, sub_group)

```
```{r}
### you should find the output differential gene file in your working directory with the suffix '*_allgenes.csv', and read it
res2 <- read.csv ('C:/projects/TAIL/insASAY_vs_WT_allgenes.csv', header = T)
res2 <- res2 %>%
  mutate(
    Change = case_when(
      padj < 0.05 & abs(log2FoldChange) > 0.5 & log2FoldChange > 0.5 ~ "UP",
      padj < 0.05 & abs(log2FoldChange) > 0.5 & log2FoldChange <= 0.5 ~ "DOWN"
    )
  )

volcano <- plot_volcano3(res2, 0.05, nlabel = 30, label.by = "padj")
volcano
```


Now replot extended fig9b, for this one, we need to get the tpm value of the count matrix
first we need to download mm10 gtf file from:  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz

```{r}
library(GenomicFeatures)
```
```{r}
txdbmm <- makeTxDbFromGFF(file = 'C://projects/genome_sequence/gencode.vM10.annotation.gtf')
```

get the transcript length
```{r}
mm10_transcript_len <- transcriptLengths(txdbmm, with.cds_len=TRUE,
                  with.utr5_len=FALSE, with.utr3_len=FALSE)
```

get the longest transcript length
```{r}
longest_transcripts <- list()

# Iterate over unique gene_id values
unique_genes <- unique(mm10_transcript_len$gene_id)
for (gene_id in unique_genes) {
  # Subset data for the current gene_id
  gene_data <- mm10_transcript_len[mm10_transcript_len$gene_id == gene_id, ]
  
  # Find the index of the longest transcript
  longest_index <- which.max(gene_data$tx_len)
  
  # Extract the information for the longest transcript
  longest_transcript <- gene_data[longest_index, c("tx_id", "tx_name", "gene_id", "nexon", "tx_len", "cds_len")]
  
  # Store the information in the list
  longest_transcripts[[gene_id]] <- longest_transcript
}

# Convert the list of data frames to a single data frame
result_df <- do.call(rbind, longest_transcripts)
```

now perform the tpm normalization

select the genes bound by Tbxt
```{r}
filtered_genes <- counts_annotated_len %>%
  filter(trimws(counts_annotated_len$GeneName) %in% pnas$gene)
```

now get all the differentially expressed genes in the three genotypes vs WT.
```{r}
###those files are output differnential genes cutted by padj < 0.05
  res1 <- read.csv ('C:/projects/TAIL/hetero_DelE6_vs_WT_padj_cutoff.csv', header = T)
res2 <- read.csv ('C:/projects/TAIL/homo_DelE6_vs_WT_padj_cutoff.csv', header = T)
res3 <- read.csv ('C:/projects/TAIL/insASAY_vs_WT_padj_cutoff.csv', header = T)
```
```{r}
res1_diff_gene <- res1[abs(res1$log2FoldChange) > 0.5, ]
res2_diff_gene <- res2[abs(res2$log2FoldChange) > 0.5, ]
res3_diff_gene <- res3[abs(res3$log2FoldChange) > 0.5, ]
```
### get raw counts with gene length and genename
```{r}
# Ensure row names are part of the data
data$RowName <- rownames(data)
result_df$RowName <- rownames(result_df)

# Merge the datasets based on this new RowName column
merged_data <- merge(data, result_df, by.x = "RowName", by.y = "RowName")

# If you no longer need the RowName column, you can remove it
merged_data$RowName <- NULL

```

```{r}
counts_annotated_len <- merge(merged_data, my_annotation, by = "gene_id", all = TRUE)
```

```{r}
# Assuming res1, res2, and res3 are your data frames
combined_diff_genes <- unique(c(res1_diff_gene$GeneName, res2_diff_gene$GeneName, res3_diff_gene$GeneName))
```
```{r}
filtered_genes <- filtered_genes %>%
  filter(trimws(filtered_genes$GeneName) %in% trimws(combined_diff_genes) )
```

```{r}
filtered_genes <- filtered_genes[, c("GeneName", "X23022FL.07.01.08_S8_L008", "X23022FL.07.01.09_S9_L008", 
                                     "X23022FL.07.01.12_S12_L008", "X23022FL.07.01.13_S13_L008", 
                                     "X23022FL.07.01.14_S14_L008", "X23022FL.07.01.15_S15_L008", 
                                     "X23022FL.07.01.10_S10_L008", "X23022FL.07.01.11_S11_L008", "tx_len")]

```

now calculate the tpm value
```{r}
rownames(filtered_genes) <- filtered_genes[,1]
counts <- filtered_genes[,2:9]
x <- counts/(filtered_genes$tx_len/1000)
```

```{r}
filtered_genes_tpm <- t( t(x) * 1e6 / colSums(x) )
```
```{r}
### now each column sum to 1M
colSums(filtered_genes_tpm)
```

```{r}
colnames(filtered_genes_tpm) <- c("WT1", "WT2", "insASAY1", "insASAY2","hete_DelE61","hete_DelE62","homo_DelE61", "homo_DelE62")
```
```{r}
write.csv(filtered_genes_tpm,"Diff_genes_bound_by_Tbxt_tpm.csv")
```


cluster genes across samples and plot heatmap


```{r}
target_tpm1 <- log10(filtered_genes_tpm)
```
```{r}
target_tpm <- scale(t(filtered_genes_tpm))
```
```{r}
tpm_z <-t(target_tpm)
```
```{r}
km <- kmeans(tpm_z,4)
```
```{r}
cluster <- km$cluster
```
```{r}
m2 <- cbind(tpm_z,cluster)
o <- order(m2[, 9])
m2 <- m2[o, ]
```

```{r}

pheatmap(m2[,1:8], cluster_rows=F,cluster_cols=F,border_color=TRUE,color=colorRampPalette(c( "dodgerblue4","lightblue3", "white", "tomato", "tomato4"))(200))
```

```{r}
pheatmap(m2[,1:8], cluster_rows=F,cluster_cols=F,border_color=TRUE)
```
