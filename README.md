# Tail-Loss-Primates
This repository contains scripts and processed data used for the project studying the genetic basis of tail-loss evolution in humans and apes (hominoids). 

The related manuscript information:  

Title: On the genetic basis of tail-loss evolution in humans and apes

Bo Xia1,2,3,4†, Weimin Zhang2*, Guisheng Zhao1,2*, Xinru Zhang3,5*, Jiangshan Bai3, Ran Brosh2, Aleksandra Wudzinska2, Emily Huang2, Hannah Ashe2, Gwen Ellis2, Maayan Pour1,2, Yu Zhao2, Camila Coelho2, Yinan Zhu2, Alexander Miller6, Jeremy S. Dasen6, Matthew T. Maurano2, Sang Y. Kim7, Jef D. Boeke2,8,9† and Itai Yanai1,2,8†

1 Institute for Computational Medicine, NYU Langone Health, New York, NY 10016, USA
2 Institute for Systems Genetics, NYU Langone Health, New York, NY 10016, USA
3 Gene Regulation Observatory, Broad Institute of MIT and Harvard, Cambridge, MA, 02142, USA
4 Society of Fellows, Harvard University, Cambridge, MA, 02138, USA
5 Department of Biology, Pennsylvania State University, University Park, PA, 16802, USA
6 Department of Neuroscience and Physiology, NYU Langone Health, New York, NY 10016, USA
7 Department of Pathology, NYU Langone Health, New York, NY 10016, USA
8 Department of Biochemistry and Molecular Pharmacology, NYU Langone Health, New York, NY 10016, USA
9 Department of Biomedical Engineering, NYU Tandon School of Engineering, Brooklyn,
NY, 11201, USA

*These authors contributed equally.

†Correspondence: xiabo@broadinstitute.org; Jef.Boeke@nyulangone.org; Itai.Yanai@nyulangone.org


## Overview
### 01_identify_mutations
In order to identify variants associated with the tail loss phenotype, we conducted a comparative phylogenetic study to scan for hominoid-specific variants in 140 genes (+/-10kb) related to tail development. The full gene list is attached in file gene140_location.csv (adapted from Ensembl BIOMART). Tbxt gene is used as an example for the pipeline.  
Taking Tbxt gene as an example:  
First extract the multiple sequence alignments for each gene. Whole genome alignment can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/. Eight sequences including six hominoid species ("hg38", "gorGor5", "panTro5", "panPan2", "ponAbe2", "nomLeu3") and  two non-hominoid species("macNem1","calJac3") are extracted. Example output is TBXT_new_8species.fasta
```
Rscript 140gene_fasta_new_8species.R "$chromosomeID"  
```
Then hominoid-specific variants can be identified by
```
python mutation.py TBXT_new_8species.fasta   
```
For 140 genes, after gathering fasta file for each gene, run:  
```
bash forloop_python.sh
```
### 02_classify_mutations
inputs/ contain 140 csv files generated from step 1 (identify mutations). For each of the 140 genes, variants are classified into insertions, deletions and SNPs. Each type of variant is stored as a bed file in the outputs/. 
```
Rscript mutation_classifier.R
```
### 03_predict_via_vep
First concatenate all the SNPs, deletions and insertions among 140 genes
```
cat 02_classify_mutations/*_snp.bed > total_snp.bed  
cat 02_classify_mutations/*_deletion.bed > total_del.bed  
cat 02_classify_mutations/*_insertion.bed > total_ins.bed
```
Outputs are saved in  03_predict_via_vep/.
### 04_filter_vep_results
```
Rscript filter_vep_visualization.R
```
Output files and plot are stored in 04_filter_vep_results/

### 05_bulk_RNA_seq
This folder contains files for bulk RNA-seq analysis presented in Extended Data Fig 9. 



