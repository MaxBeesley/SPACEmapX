# SPACEmapX Tutorial Page

This page is currently under development and will soon be published as a standard version.

# SPACE Lab
??????

# Introduction
Spatial transcriptomics technology provides in situ whole-transcriptome data localised to high resolution tissue regions. It is beginning to transform our understanding of somatic events underlying the development of cancer, compared with previous bulk or even single cell sequencing technology.

Copy number alterations (CNA) are a common mutational feature of cancer development, and can be used to interrogate clonal phylogenetics in pre-cancer and cancerous cells. inferCNV is a computational method to infer the copy number variations from transcriptomic data. It was initially developed for use in single cell RNAseq data (Patel et al., Science 2014 10.1126/science.1254257; github.com/broadinstitute/infercnv) and we have repurposed it for use with the 10x Genomics Visium spatial transcriptomics platform, previously developing a package called SpatialInferCNV (github.com/aerickso/SpatialInferCNV).
Here we further extend the SpatialInferCNV package to incorporate the latest Visium spatial transcriptomic chemistries (Visium v2 and Visium HD), providing a guide to do multi-section patient ST analyses through these sets:
    Step 1 - Preparation of input data
    Step 2 - Reference set selection
    Step 3 - Section-specific clone selection
    Step 4 - Global clone selection
    Step 5 - Generate phylogenetic tree
    Step 6 - Overlay with spatial visualisations 

For this demonstration, we have applied SPACEmapX to 10 patients including both prostate (P) and lymph node (LN) tissue. Each patient has approximately 10 sections (11mm2 or 6.5mm2), totalling 2 million transcriptomic data points (spots). 

# Quick Start
To analyse copy number variation (CNV) in spatial transcriptomics data, the two necessary inputs are gene expression matrix (filtered_feature_bc_matrix.h5) and an annotation file stating the cell type present at each spot as determined from histology (.csv file exported from the 10x Genomics Loupe browser as discussed below) or from scRNAseq-style transcriptomic analysis. 

## Software Requirements 
R environment (tutorial carried out on 4.4.3)
R package SpaceMapX (tutorial carried out on 1.1.0)
R package InferCNV (tutorial carried out on 1.22.0)
Loupe browser (tutorial carried out on 8.1.2).
!!!!! What RAM / nodes are needed for this analysis?? !!!!!

## Histological annotation of data spots
!!!!!!!!! describe out histology approach

## Downsampling
????









# Work-Flow
![workflow1](https://github.com/MaxBeesley/SPACEmapX/blob/main/Figures/SPACEmapX_github_workflow_MB.png)










## Step 1 - Preparation of input data
For analysis of multiple samples, a .csv file is required providing sample name (SampleID), location of gene expression data (SampleGeneExpression), location of Loupe browser annotation file (AnnotationFile) and tissue type (Type). Name this .csv file as "SpatialDataTable.csv"

![image](https://github.com/MaxBeesley/SPACEmapX/blob/main/Figures/dataset.png)

The relevant R code is shown below each step

``` r
## Example R script
# load the required packages and set paths
library(infercnv)
library(tidyverse)
library(dplyr)
library(Seurat)
library(tibble)
library(SPACEmapX)

working_directory="<working_directory>"
setwd(working_directory)
```

Step 1 - Input section data

Read the expression matrix into R using the automated SpatialDataTable.csv method and then merge the histology and count matrices.

```r

# load source data to analysis (example below from wiki)
SpatialDataTable <- SPACEmapX::loadSTinfo("SpatialDataTable.csv")

# load all expression matrix and barcode lists
SPACEmapX::loadSTdata(SpatialDataTable)


## Save histology and counts for each section to run initial analysis (ref=NULL) in order to pick reference barcodes
## We suggest that appropriate filtering is conducted to make cell types of interest equally represented in the dataset as no_ref inferCNV analysis relies on normalisation against dataset averages for specific gene expression (i.e. susceptible to being skewed by over-/under-representation)
counts_A2 <- SpatialInferCNV::MergingCountAndAnnotationData("A2", A2_histology, count_matrix_A2) # filter counts for barcodes shared with histology
annotations_A2 <- SpatialInferCNV::FinalAnnotations(A2_histology, count_matrix_A2) # filter annotations for barcodes shared with count matriix
write.table(counts_A2, "counts_A2.tsv", sep = "\t")
write.table(annotations_A2, "annotations_A2.tsv", 
            sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

counts_P4 <- SpatialInferCNV::MergingCountAndAnnotationData("P4", P4_histology, count_matrix_P4)
annotations_P4 <- SpatialInferCNV::FinalAnnotations(P4_histology, count_matrix_P4)
write.table(counts_P4, "counts_P4.tsv", sep = "\t")
write.table(annotations_P4, "annotations_P4.tsv", 
            sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

counts_LII1 <- SpatialInferCNV::MergingCountAndAnnotationData("LII1", LII1_histology, count_matrix_LII1)
annotations_LII1 <- SpatialInferCNV::FinalAnnotations(LII1_histology, count_matrix_LII1)
write.table(counts_LII1, "counts_LII1.tsv", sep = "\t")
write.table(annotations_LII1, "annotations_LII1.tsv", 
            sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

```

## Step 2 - Reference set selection
### a) Selecting section-specific reference clones
!!!!!!! needs work. & explan 'twig' terminology
a)	Analyse each individual section with SpatialInferCNV? without using a Reference.
b)	Use the "twig ID"? function to generate a section-specific heatmap
c)	Select the quietest? twig for each section.
### b) Selecting patient-specific reference clone
d)	Merge all of these ‘quietest twigs’ together and re-analyse with SpatialInferCNV? it without reference.
e)	Now assess this new heatmap and select the quietest twig overall which will now be used as the reference for all sections from this patient.

```r
## Step 2a - Section-specific reference set selection
## Run inferCNV with no reference to select the 'quiet twig' from this section
## siCNV_GeneOrderFile.tsv chromosome order reference file available at "raw.githubusercontent.com/aerickso/SpatialInferCNV/main/FigureScripts/siCNV_GeneOrderFile.tsv"

# Running inferCNV with no reference for each section, initially section A2
inferCNV_object <- infercnv::CreateInfercnvObject(
                                raw_counts_matrix = "counts_A2.tsv",
                                annotations = "annotations_A2.tsv",
                                gene_order_file = "siCNV_GeneOrderFile.tsv", ## Reference file - link available above
                                delim = "\t",
                                ref_group_names = NULL, ## No reference used on this run
                                chr_exclude = c("chrM")) ## Exclude M genes

inferCNV_object = infercnv::run(inferCNV_object,
                                cutoff=0.1, ## 0.1 for 10X Genomics data
                                out_dir="<out_dir>/A2_RefSelection", 
                                num_threads = 2, ## Number of computinng threads to be used
                                plot_chr_scale=TRUE, ## Scale size of chromosomes based on actual size
                                write_phylo = TRUE, ## Write phylogenetic tree as a file
                                window_length = 101, ## Smoothing window size od 101 genes (50 genes either side of query gene)
                                hclust_method='ward.D2',
                                denoise=TRUE,
                                cluster_by_groups=FALSE, 
                                analysis_mode = "samples",
                                tumor_subcluster_partition_method = "qnorm",
                                HMM=FALSE,
                                output_format = "png")


# Assessing phylogenetic dendogram and selecting barcodes from twigs of interest to select section-specific reference
clustering <- SPACEmapX::read.dendrogram(file="infercnv.preliminary.observations_dendrogram.txt")
clustering_phylo <- SPACEmapX::as.phylo(clustering)
my.subtrees = SPACEmapX::subtrees(clustering_phylo)  # subtrees() to subset

png("clustering_phylo.png",width=10000,height=10500, res = 300)
plot(clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:clustering_phylo$Nnode,node=1:clustering_phylo$Nnode+Ntip(clustering_phylo))
dev.off()

# Select quiet twig, e.g. 321 and generate a .csv file of the barcodes present in that twig
SPACEmapX::ShowTwigSpotsList(321)

## Repeat step 2 for each of the sections which results in list of barcodes representing the quiet barcodes across the patient
```
```r
## Step 2b - Patient-specific reference set selection
## Run inferCNV with no reference to select the 'quietest of the quiet twig' from this section
## Combine the data from the quiet twigs of each section and reanalyse with inferCNV without a reference to find the patient-wide quietest twig.
## This is shown for the 3 slices (A2, P4, LII1) after Step 2a has been conducted for each slice to pick the quiet twig and import the *_Ref.csv generated by ShowTwigSpotsList()

# Import A2_Ref and P4_Ref .csv files containing lists of quiet barcodes as generated above.
A2_Ref <- read.csv(A2_Ref.csv)
P4_Ref <- read.csv(P4_Ref.csv)
LII1_Ref <- read.csv(lII1_Ref.csv)

# Filter count matrices for selected barcodes
A2_ref_filtered <- SpatialInferCNV::MergingCountAndAnnotationData("A2", A2_Ref, count_matrix_A2)
P4_ref_filtered <- SpatialInferCNV::MergingCountAndAnnotationData("P4", P4_Ref, count_matrix_P4)
LII1_ref_filtered <- SpatialInferCNV::MergingCountAndAnnotationData("LII1", LII1_Ref, count_matrix_LII1)

# Merge histology and count data and save for use in inferCNV
ref_counts_joined <- merge(A2_ref_filtered, P4_ref_filtered, by="Genes") !!!!! check
ref_counts_joined <- merge(ref_counts_joined, LII1_ref_filtered, by="Genes") !!!!! check
ref_counts_joined <- ref_counts_joined %>% column_to_rownames(., var = "Genes")
write.table(ref_counts_joined, "ref_counts_joined.tsv", sep = "\t")
ref_histology_joined <- rbind(A2_Ref, P4_Ref, LII1_Ref)
ref_histology_joined <- SpatialInferCNV::FinalAnnotations(ref_histology_joined, ref_counts_joined)
write.table(ref_histology_joined, "ref_histology_joined.tsv", 
            sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

inferCNV_object <- infercnv::CreateInfercnvObject(
                                raw_counts_matrix="ref_counts_joined.tsv", 
                                annotations_file="ref_histology_joined.tsv",
                                gene_order_file="siCNV_GeneOrderFile.tsv",
                                delim="\t",
                                ref_group_names= NULL,
                                chr_exclude = c("ChrM"))

inferCNV_object = infercnv::run(inferCNV_object,
                                cutoff=0.1, ## 0.1 for 10X Genomics data
                                out_dir="<out_dir>/Patient_RefSelection", 
                                num_threads = 2, ## Number of computinng threads to be used
                                plot_chr_scale=TRUE, ## Scale size of chromosomes based on actual size
                                write_phylo = TRUE, ## Write phylogenetic tree as a file
                                window_length = 101, ## Smoothing window size od 101 genes (50 genes either side of query gene)
                                hclust_method='ward.D2',
                                denoise=TRUE,
                                cluster_by_groups=FALSE, 
                                analysis_mode = "samples",
                                tumor_subcluster_partition_method = "qnorm",
                                HMM=FALSE,
                                output_format = "png")


# Once again select the quietest twig to pick the patient-specific reference barcodes representing the quietest of the quiet 
# Assess phylogenetic dendogram and select barcodes from twigs of interest
clustering <- read.dendrogram(file="infercnv.preliminary.observations_dendrogram.txt")
clustering_phylo <- as.phylo(clustering)
my.subtrees = subtrees(clustering_phylo)  # subtrees() to subset

png("clustering_phylo.png",width=10000,height=10500, res = 300)
plot(clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:clustering_phylo$Nnode,node=1:clustering_phylo$Nnode+Ntip(clustering_phylo))
dev.off()

# Select quiet of the quietest twig, e.g. 123 and generate a ref.csv file of the barcodes present in that twig for use as patient-specific reference
ShowTwigSpotsList(123)


## Reference is now selected. Analysis can now continue with intitally each section and then patient-wide. Below analysis is shown on a global patient-wide scale as the pipeline is the same. 
```


## Step 3 - Section-specific clone selection
a)	Analyse each section using the reference identified during Step 2) [check below:pathological interested spots/cluster based spots]
b)	Identify clusters with common CNV features and annotate as a candidate clone.
c)	Run "twig ID"? function to generate global clone heatmap
d)	Select different clone using twig number?.
N.B. Key requirements to consider a region a 'clone'
Same copy number features (part manual visualisation; part dendogram)
Located in the same spatial location (can check with Loupe Browser)
Ten or more adjacent spots in one location


## Step 4 - Global clone selection
a)	SpaceMapX is used to analyze multiple sections of 10x ST slides. The overall number of spots can therefore be very large, making it challenging in terms of memory limitation and computation time. Some organisation-controlled workstations do not allow researchers to install DisplayX to remotely display images, and due to Linux package limitations (such as with the Cairo package which will only work with less than 40K spots on PNG output). It's important to find ways to reduce the time and spots size needed to combine all clones and generate a global clone heat map. 
b)	For large clones with more than 100 spots, we randomly select 100 spots + 10%. For clones of less than 100 spots, we take all the spots. This enables a reprresentaive sampling of the data. 
This is carried out in the R function: sample_n(CloneList, 100 + ceiling(total_rows * 0.1))
c)	Then combine all clones to plot a global clone heat map across each spatial slice from this patient.
d)	Identify each region of common inferCNV features and assign to a candidate clone 
> Naming convention: primary clone name = a-z; metastatic clone name = X1, X2, X3 etc.
e)	Re-generate inferCNV heatmap with clone order forced according to assigned clone names

```r
## Step 4 - Global clone selection

# Initially inferCNV is re-run, now using the chosen reference
ref <- read.csv(ref.csv)

# Generate patient-wide counts and histology matrices
histology_annotations <- rbind(A2_histology, P4_histology, LII1_histology))

## Select for cell types of interest
keep_types <- c("Benign","Cancer","Cancer 3+4") # Needs editing based on desired cell types for analysis
histology_annotations <- histology_annotations[histology_annotations$Histology %in% keep_types, ]

## Generate combnined count matrix
count_matrix <- merge(count_matrix_A2, count_matrix_P4, by = "genes")
count_matrix <- merge(count_matrix, count_matrix_LII1, by = "genes")
rownames(count_matrix) <- count_matrix$genes
count_matrix <- count_matrix[,-1]
count_matrix <- as.data.frame(count_matrix)
dim(count_matrix)

## Save combined counts for final analyses 
final_counts <- SpatialInferCNV::MergingCountAndAnnotationData(sample_name, histology_annotations, count_matrix)
write.table(final_counts, "final_counts.tsv", sep = "\t")

## Set reference barcodes in the combined patient histology and save for final analyses 
annotations_final <- SpatialInferCNV::FinalAnnotations(histology_annotations, count_matrix)
annotations_final[(annotations_final$Barcode %in% ref$Barcode),"Histology"] <- "reference"
write.table(annotations_final, "annotations_final.tsv", 
            sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)


# Now ref_group_names is set to "reference"
inferCNV_object <- infercnv::CreateInfercnvObject(
                                raw_counts_matrix="final_counts.tsv", 
                                gene_order_file="siCNV_GeneOrderFile.tsv",
                                annotations_file="annotations_final.tsv",
                                delim="\t",
                                ref_group_names= "reference", # set reference 
                                chr_exclude = c("ChrM"))


inferCNV_object = infercnv::run(inferCNV_object,
                                cutoff=0.1, ## 0.1 for 10X Genomics data
                                out_dir="<out_dir>/final_analyses", 
                                num_threads = 2, ## Number of computinng threads to be used
                                plot_chr_scale=TRUE, ## Scale size of chromosomes based on actual size
                                write_phylo = TRUE, ## Write phylogenetic tree as a file
                                window_length = 101, ## Smoothing window size od 101 genes (50 genes either side of query gene)
                                hclust_method='ward.D2',
                                denoise=TRUE,
                                cluster_by_groups=FALSE, 
                                analysis_mode = "samples",
                                tumor_subcluster_partition_method = "qnorm",
                                HMM=FALSE,
                                output_format = "png")


# Generate dendogram of clones and manually select barcodes in each clone for downstream analysis
clustering <- read.dendrogram(file="infercnv.preliminary.observations_dendrogram.txt")
clustering_phylo <- as.phylo(clustering)
my.subtrees = subtrees(clustering_phylo)  # subtrees() to subset

png("clustering_phylo.png",width=10000,height=10500, res = 300)
plot(clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:clustering_phylo$Nnode,node=1:clustering_phylo$Nnode+Ntip(clustering_phylo))
dev.off()

## Now conduct the manual selection of Clones from the CNV data in the inferCNV heatmap. Saving these using ShowTwigSpotsList() and naming these to what Clone they are. Clone X is used to designate clones that are metastatic.
# E.g. select the three twigs belonging to Clone C
ShowTwigSpotsList(112, 444, 677)

```
```r
## inferCNV is run a final time to produce the final heatmap whereby the barcodes are ordered based on Clones and each clone is annotated with a coloured bar on the Y-axis. For this, cluster_by_groups is changed to =TRUE

selected_clones <- c(cloneA, CloneB, CloneC, CloneD, CloneE, CloneF, CloneG, CloneH, CloneX) !!check
annotations_clonal <- histology_annotations[histology_annotations$Barcode %in% selected_clones$Barcode,]

annotations_clonal <- SpatialInferCNV::FinalAnnotations(annotations_clonal, count_matrix)
annotations_clonal[(annotations_clonal$Barcode %in% ref$Barcode),"Histology"] <- "reference"
write.table(annotations_clonal, "annotations_clonal.tsv", 
            sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

counts_clonal <- SpatialInferCNV::MergingCountAndAnnotationData("clonal", annotations_clonal, count_matrix)
write.table(final_counts, "counts_clonal.tsv", sep = "\t")

inferCNV_object <- infercnv::CreateInfercnvObject(
                                raw_counts_matrix="counts_clonal.tsv", 
                                gene_order_file="siCNV_GeneOrderFile.tsv",
                                annotations_file="annotations_clonal.tsv",
                                delim="\t",
                                ref_group_names= "reference", # set reference 
                                chr_exclude = c("ChrM"))

# Change cluster_by_groups to =TRUE, which will now reorder the heatmap to group clones on Y-axis and enable better visualisation
inferCNV_object = infercnv::run(inferCNV_object,
                                cutoff=0.1, ## 0.1 for 10X Genomics data
                                out_dir="<out_dir>/clonal_analyses", 
                                num_threads = 2, ## Number of computing threads to be used
                                plot_chr_scale=TRUE, ## Scale size of chromosomes based on actual size
                                write_phylo = TRUE, ## Write phylogenetic tree as a file
                                window_length = 101, ## Smoothing window size od 101 genes (50 genes either side of query gene)
                                hclust_method='ward.D2',
                                denoise=TRUE,
                                cluster_by_groups=TRUE, 
                                analysis_mode = "samples",
                                tumor_subcluster_partition_method = "qnorm",
                                HMM=FALSE,
                                output_format = "png")

# Next, assess the inferCNV outputs to produce phylogenetic clonal tree and spatial visualisations

```
## Step 5 - Generate phylogenetic clonal tree
This analysis enables a phylogenetic clonal tree to be generated by manual assessment of how each of the clones is related to one another by interpreting the CNV data. 


## Step 6 - Overlay with spatial visualisations 
Import inferCNV .csv file to Loupe Browser and visualise on spatial map.

Conducted using the Loupe browser by importing the inferCNV .csv file with two columns (Barcode and Clonal_annotation) and visualising this onto the spatial map. This overlay can also be achieved using Seurat or ggplot2 (using the annotation_raster function).

## Example of a final heatmap with different clones annotated
![image](https://github.com/yintz/SPACEmapX/blob/main/Figures/P15.preliminary.png)






# 10X Visium HD ST analysis

The improved 10x FFPE V2 HD contains more than 42 million barcodes in 2um squares. The barcode naming system has changed and also the size of the expression matrix, which often reaches a file size of 92GB - therefore a platform with around 200GB RAM is required. Broad analysis approach should still be applicable




# Citation
Will be added once article is published

# Funding 
Cancer Research UK (CRUK) #C57899/A25812 "Spatial Prostate Assessment and Circulating Environment – The SPACE Study"


