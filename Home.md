# SPACEmapX Tutorial Page

This page is currently under development and will soon be published as a standard version.


# Introduction
Spatial transcriptomics technology provides in situ whole-transcriptome data localised to high resolution tissue regions. It is beginning to transform our understanding of somatic events underlying the development of cancer, compared with previous bulk RNA or even single cell technology. 

Copy number alterations (CNA) are a common feature of cancer development, and can be used to interrogate clonal phylogenetics in pre-cancer and cancerous cells. inferCNV is a computational method to infer the copy number variations from transcriptomical data. It was initially developed for use in single cell RNAseq data (Patel et al, Science 2014; https://github.com/broadinstitute/infercnv) and we have repurposed it for use with visium spatial transcriptomics, developing a package called spatialinferCNV (https://github.com/aerickso/SpatialInferCNV).

Here we extend this package to incorporate the latest visium spatial transcriptomic chemistries (Visium v2 and Visium HD), providing a guide to undertake the following analyses
1. Selection of appropriate patient-specific reference set
2. Generation of a "global" inferCNV heatmap for multiple sections from a single patient
3. Construction of a unifying heatmap including suitable legends corresponding to section of origin, tissue type and clone membership
4. Deduction of a phylogenetic tree reprentative of the clones demonstrated in step 2 and 3
5. Visualisation of clones in a spatial summary of data

For illustration, we have applied SPACEmapX to 10 patients including both prostate and lymph node tissue. Each patients had approximately 10 sections (11mm or 6.5mmsq), totallying 2million transcriptomical data points (spots). Of note, inferCNV is limited to approximately 40,000 data points for plotting, so we demonstrate how to combine 100,000+ spots from a single patient by random selection of spots from each "clone" to facitate inclusion into a "global heatmap" prior to generation of phylogenetic trees and spatial data visulatisation.

# Quick Start
spatial transcriptomics copy number variation needs gene expression matrix and annotation file for each spots. 
they name as filtered_feature_bc_matrix.h5 and csv files(exported from loupe browser(all data points annotated by pathologist or other methods to tell the difference of the different dataset.)

# Software Requirements 
R,  
Loupe browser,  
workstation(either Windows/linux based)  


# Work-Flow
![workflow1](https://github.com/yintz/SPACEmapX/blob/d934cbc53e177f392432acbbc54da0ea0bd2abea/Figures/graph.png)





Preparation of data
* Construct matrix of all tissue sections to be analysed
* Include: filename; filesource; section name (e.g. A1, A2, A3, P1, etc); tissue type (e.g. prostate or lymph node)

# Data requirement
for multiple samples analysis, it needs a csv file build by users to have these information as input.
![](man/DataDemo.png)
![](https://github.com/yintz/SPACEmapX/tree/3aa863eb1cc5eb2c51a320bec94889f7e506eca9/man)
![](R/dataset.png)
![Alt text of the image](https://github.com/yintz/SPACEmapX/blob/aade49968c3acd8105eb3e1ad01468de1d30101b/Figures/dataset.png)

[
](https://github.com/yintz/SPACEmapX/blob/aade49968c3acd8105eb3e1ad01468de1d30101b/Figures/dataset.png)
Save this File as "SpatialDataTable.csv"


in Rstudio load SPACEmapX package
use LoadSTinfo('named.csv') to load all collective expression matrixes and related info.


1. Reference set selection
in order to run copy number variation, we need select a pure clean reference on diploid status. we use Run() to 
run all sections in original tissue.

* using dendrogram to select the quietest twigs(the twig contains the least red/blue pattern across all chromosome on figure.
save all the barcode list info. (this is the first round of quiet reference spots barcode lists)
* recreate a csv files to contain all matrix expression with the barcode selected.
use LoadSTinfo('named.csv') to load all data and merge all quite reference spots barcode together then use run() to generate a new heat map.
* on this heat map, to select the quietest twig and save the new selected barcode list as final reference for this whole dataset.


2. Section-specific clone selection
* run each section using reference set from 1) [check below:pathological interested spots/cluster based spots]
* identify each region with common iCNV features and annotate to candidate clone.
* run "twig ID" function for global clone heatmap
* select different clone using twig number.

N.B. Key requirements to consider a region a 'clone'
> * Same copy number features (part manual visualisation; part dendogram)
> * Located in the same place on Loupe Browser
> * Ten or more spots in that location

3. Global clone selection
* SpaceMapX is used to analyze multiple sections of 10x ST slides. The overall number of spots can be very large clones, making it challenging in terms of memory limitation and computation time. Some organisation-controlled workstations do not allow researchers to install DisplayX to remotely display images, and due to Linux package limitations (such as with the Cairo package(only work with less than 40K spots on PNG output). It's important to find ways to reduce the time and spots size needed to combine all clones and generate a global clone heat map. For large clones with more than 100 spots, we randomly select (**R** function **CloneList <- sample_n(CloneList, 100 + ceiling(total_rows * 0.1))** ) 100 spots and 10% of the total clone spots to represent the clone size differences in default setting.

* combine all clone to plot a global clone heat map.
* identify each region of common iCNV features and assign to candidate clone 
(naming convention: primary clone name = a-z; metastatic clone name = X1, X2, X3 etc)
* re-run iCNV heatmap with clone order forced according to assigned clone names

4. Generate phylogenetic tree

The SPACEmapX generates a phylogenetic file SPACEmapX.preliminary.observations_dendrogram.txt
 to further clone selection and refine.

5. Create spatial visualisations
* Use iCNV .csv file to import to Loupe Browser
* Create spatial visualisation of global clone map for figure generation. 


## Patient samples
![image](https://github.com/yintz/SPACEmapX/blob/main/Figures/results.png)


### Clusters spots based method
for 10x FFPE V2 slides, 11mm has around 14k spots,(It is a time-consuming job for pathologies to annotate all spots for benign, cancer stages, stroma, etc).and in this test, normally, it requires many sections, like several from prostate local cancer. several lymph nodes sections. it is even harder for most of lab. if you don't have such requirement, you can run analysis based on cluster and gene expression.



### Pathology annotation spots based method
if you have pathlogist to annotate all spots for each section. we recommended to have multiple pathologists to annotate and compare the different spots for double verification.

based on this pathologist annotation, it has more accurate spots and less number of spots to reduce the computation power.


# 10X ST HD slides copy number variation methods

10x FFPE V2 HD is the high resolution HD, it contains more than 42 million 2um square. the barcode naming system has changed and also the size of the expression matrix is around 92G high, so I recommend to use a workstation has at least 200G memory to handle the data.

because of the reading is relatively low in each 2um square now. I removed the criteria UMI count >500. 


I use 8um H5 expression matrix(contain loupe browser file) as input file. based on this, we can change to magnitude of 8um, like 16um 24um etc. I also created a function change the size of bin.

# Script

``` r
# load the required packages
install.packages("remotes")
remotes::install_github("yintz/SPACEmapX")
library(devtools)
library(ape)
library(phylogram)
library(tidyverse)
library(preprocessCore)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(remotes)
library(devtools)
library(ape)
library(phylogram)
library(tidyverse)
library(hdf5r)
library(Cairo)
devtools::install_github("broadinstitute/infercnv@RELEASE_3_16")

```

``` r
# set working directory
setwd("./Pt15")

``` 
Reads the expression matrix into R
if you have CSV file created as above.

```r

# load source data to analysis (example below from wiki)
infotable<-loadSTinfo("/data/SpatialDataTable.csv")

# load all expression matrix and barcode lists
loadSTdata(infotable)



histologylist<-rbind(A2_histology,P4_histology,LII1_histology)
histology_df <- as.data.frame(table(histologylist$Histology))
colnames(histology_df) <- c("Histology_Type", "Count")
print(histology_df)

keep_types <- c("Benign","Cancer","Cancer 3+4")
histologylist <- histologylist[histologylist$Histology %in% keep_types, ]


A2_Joined_Counts <- MergingCountAndAnnotationData("A2",A2_Histology, A2_ENSBMLID_Counts)
P4_Joined_Counts <- MergingCountAndAnnotationData("P4",P4_Histology, P4_ENSBMLID_Counts)
LII1_Joined_Counts <- MergingCountAndAnnotationData("LII1",LII1_Histology, LII1_ENSBMLID_Counts)

Counts_joined <- A1_Joined_Counts %>% replace(., is.na(.), 0)
Counts_joined <- Counts_joined %>% column_to_rownames(., var = "Genes")
write.table(Counts_joined, "countMatrix.tsv", sep = "\t")
MergeAllAnnotation<-A1_Joined_Counts
MergedAll_Final <- FinalAnnotations(MergeAllAnnotation, Counts_joined)
write.table(MergedAll_Final, "Annotations.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

Ref_SpaceMapX <- SpaceMapX::CreateSpaceMapXObject(raw_counts_matrix="countMatrix.tsv", 
                                             gene_order_file="./siCNV_GeneOrderFile.tsv",
                                             annotations_file="./Annotations.tsv",
                                             delim="\t",
                                             ref_group_names= NULL,
                                             chr_exclude = c("ChrM"))

Ref_SpaceMapXNoGroup = SpaceMapX::run(Ref_SpaceMapX,
                                  cutoff=0.1,
                                  out_dir="./CloneA1", 
                                  num_threads = 10,
                                  denoise=T,
                                  output_format = "png",
                                  hclust_method='ward.D2',
                                  cluster_by_groups=F,
                                  analysis_mode = "samples",
                                  tumor_subcluster_partition_method = "qnorm",
                                  HMM=F)


clustering <- read.dendrogram(file="./CloneA1/infercnv.preliminary.observations_dendrogram.txt")
clustering_phylo <- as.phylo(clustering)
my.subtrees = subtrees(clustering_phylo)  # subtrees() to subset
png("clustering_phylo.png",width=10000,height=10500, res = 300)
plot(clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:clustering_phylo$Nnode,node=1:clustering_phylo$Nnode+Ntip(clustering_phylo))
dev.off()


for example the quiet twig is 539
```
ShowTwigSpotsList(539)
```r
this function will generate a csv file with spots barcode.



Counts_joined <- P4_Joined_Counts %>% replace(., is.na(.), 0)
Counts_joined <- Counts_joined %>% column_to_rownames(., var = "Genes")
write.table(Counts_joined, "countMatrix.tsv", sep = "\t")
MergeAllAnnotation<-P4_Joined_Counts
MergedAll_Final <- FinalAnnotations(MergeAllAnnotation, Counts_joined)
write.table(MergedAll_Final, "Annotations.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

Ref_SpaceMapX <- SpaceMapX::CreateSpaceMapXObject(raw_counts_matrix="countMatrix.tsv", 
                                             gene_order_file="./siCNV_GeneOrderFile.tsv",
                                             annotations_file="./Annotations.tsv",
                                             delim="\t",
                                             ref_group_names= NULL,
                                             chr_exclude = c("ChrM"))

Ref_SpaceMapXNoGroup = SpaceMapX::run(Ref_SpaceMapX,
                                  cutoff=0.1,
                                  out_dir="./CloneP4", 
                                  num_threads = 10,
                                  denoise=T,
                                  output_format = "png",
                                  hclust_method='ward.D2',
                                  cluster_by_groups=F,
                                  analysis_mode = "samples",
                                  tumor_subcluster_partition_method = "qnorm",
                                  HMM=F)


clustering <- read.dendrogram(file="./CloneP4/infercnv.preliminary.observations_dendrogram.txt")
clustering_phylo <- as.phylo(clustering)
my.subtrees = subtrees(clustering_phylo)  # subtrees() to subset
png("clustering_phylo.png",width=10000,height=10500, res = 300)
plot(clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:clustering_phylo$Nnode,node=1:clustering_phylo$Nnode+Ntip(clustering_phylo))
dev.off()

for example the quiet twig is 324
```
ShowTwigSpotsList(324)
```r
this function will generate a csv file with spots barcode.

Now we have two quiet Ref. we need to merge them together to make the quietest ref as final ref.


load A2_Ref and P4_Ref csv files into Rstudo.

A2_Joined_Counts <- MergingCountAndAnnotationData("A2",A2_Ref, A2_ENSBMLID_Counts)
P4_Joined_Counts <- MergingCountAndAnnotationData("P4",P4_Ref, P4_ENSBMLID_Counts)

Counts_joined <- A2_Joined_Counts %>% P4_Joined_Counts
Counts_joined <- Counts_joined %>% replace(., is.na(.), 0)
Counts_joined <- Counts_joined %>% column_to_rownames(., var = "Genes")
write.table(Counts_joined, "countMatrix.tsv", sep = "\t")
MergeAllAnnotation<-rbind(A2_Ref,P4_Ref)
MergedAll_Final <- FinalAnnotations(MergeAllAnnotation, Counts_joined)
write.table(MergedAll_Final, "Annotations.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

Ref_SpaceMapX <- SpaceMapX::CreateSpaceMapXObject(raw_counts_matrix="countMatrix.tsv", 
                                             gene_order_file="./siCNV_GeneOrderFile.tsv",
                                             annotations_file="./Annotations.tsv",
                                             delim="\t",
                                             ref_group_names= NULL,
                                             chr_exclude = c("ChrM"))

Ref_SpaceMapXNoGroup = SpaceMapX::run(Ref_SpaceMapX,
                                  cutoff=0.1,
                                  out_dir="./CloneRef", 
                                  num_threads = 10,
                                  denoise=T,
                                  output_format = "png",
                                  hclust_method='ward.D2',
                                  cluster_by_groups=F,
                                  analysis_mode = "samples",
                                  tumor_subcluster_partition_method = "qnorm",
                                  HMM=F)


clustering <- read.dendrogram(file="./CloneRef/infercnv.preliminary.observations_dendrogram.txt")
clustering_phylo <- as.phylo(clustering)
my.subtrees = subtrees(clustering_phylo)  # subtrees() to subset
png("clustering_phylo.png",width=10000,height=10500, res = 300)
plot(clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:clustering_phylo$Nnode,node=1:clustering_phylo$Nnode+Ntip(clustering_phylo))
dev.off()

for example the quiet twig is 5543
```
ShowTwigSpotsList(5543)
```r
this function will generate a csv file with spots barcode, save this csv file / files  and load it back to Rstudio

# This is how the Reference created.


load Ref csv files into Rstudo. for example, the reference is only on A2 as A2Ref.csv

A2_Be_ENSBMLID_Counts <- ImportCountData("Be_A2", "./A2/filtered_feature_bc_matrix.h5")
Be_A2_Histology <- ImportHistologicalAnnotations("Be_A2","/A2Ref.csv")
Be_A2_Histology$Histology<-c("Pt15Ref")
benign_A2_Joined_Counts <- MergingCountAndAnnotationData("Be_A2",Be_A2_Histology, A2_Be_ENSBMLID_Counts )
selected_benign <- FinalAnnotations(Be_A2_Histology, benign_A2_Joined_Counts)

A2_Joined_Counts <- MergingCountAndAnnotationData("A2",A2_histology, A2_ENSBMLID_Counts)
P4_Joined_Counts <- MergingCountAndAnnotationData("P4",P4_histology, P4_ENSBMLID_Counts)
LII1_Joined_Counts <- MergingCountAndAnnotationData("LII1",LII1_histology, LII1_ENSBMLID_Counts)

AllEpitheliumEnriched <- rbind(A2_Histology,P4_Histology,LII1_Histology)
             
Counts_joined <- A2_Joined_Counts %>% P4_Joined_Counts
Counts_joined <- Counts_joined %>% LII1_Joined_Counts
Counts_joined <- Counts_joined %>% benign_A2_Joined_Counts
Counts_joined <- Counts_joined %>% replace(., is.na(.), 0)
Counts_joined <- Counts_joined %>% column_to_rownames(., var = "Genes")
write.table(Counts_joined, "countMatrix.tsv", sep = "\t")
MergeAllAnnotation<-rbind(AllEpitheliumEnriche,selected_benign)
MergedAll_Final <- FinalAnnotations(MergeAllAnnotation, Counts_joined)
write.table(MergedAll_Final, "Annotations.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

this time, change the ref_group_names to the name of the Ref "Pt15Ref"
Ref_SpaceMapX <- SpaceMapX::CreateSpaceMapXObject(raw_counts_matrix="countMatrix.tsv", 
                                             gene_order_file="./siCNV_GeneOrderFile.tsv",
                                             annotations_file="./Annotations.tsv",
                                             delim="\t",
                                             ref_group_names= c("Pt15Ref"),
                                             chr_exclude = c("ChrM"))


Ref_SpaceMapXNoGroup = SpaceMapX::run(Ref_SpaceMapX,
                                  cutoff=0.1,
                                  out_dir="./CloneFinalUnGroup", 
                                  num_threads = 10,
                                  denoise=T,
                                  output_format = "png",
                                  hclust_method='ward.D2',
                                  cluster_by_groups=F,
                                  analysis_mode = "samples",
                                  tumor_subcluster_partition_method = "qnorm",
                                  HMM=F)


clustering <- read.dendrogram(file="./CloneFinalUnGroup/infercnv.preliminary.observations_dendrogram.txt")
clustering_phylo <- as.phylo(clustering)
my.subtrees = subtrees(clustering_phylo)  # subtrees() to subset
png("clustering_phylo.png",width=10000,height=10500, res = 300)
plot(clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:clustering_phylo$Nnode,node=1:clustering_phylo$Nnode+Ntip(clustering_phylo))
dev.off()


separate twigs using CNV pattern. for example, it can be separated as 112, 444, 677. 3 twigs in total.
```
ShowTwigSpotsList(112, 444, 677)
```r
this function will generate a csv file with spots barcode, save this csv file as Ref.csv and load it back to Rstudio

set each csv file a clone name like clone0 for same as Ref, cloneA, cloneB, for more and more CNV patten changes, CloneX is the clone shown on both lymph nodes and also tissue. if there are more Clone on both tissues. name it as CloneX1, CloneX2, CloneX3.
load it back to Rstudio.


```r
AllEpitheliumEnriched <- rbind(cloneA,CloneB,CloneX)
MergeAllAnnotation<-rbind(AllEpitheliumEnriche,selected_benign)

MergedAll_Final <- FinalAnnotations(MergeAllAnnotation, Counts_joined)
write.table(MergedAll_Final, "Annotations.tsv", 
            sep = "\t",
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

#change the ref_group_names to the name of the Ref "Pt15Ref"
Ref_SpaceMapX <- SpaceMapX::CreateSpaceMapXObject(raw_counts_matrix="countMatrix.tsv", 
                                             gene_order_file="./siCNV_GeneOrderFile.tsv",
                                             annotations_file="./Annotations.tsv",
                                             delim="\t",
                                             ref_group_names= c("Pt15Ref"),
                                             chr_exclude = c("ChrM"))

# change the cluster_by_groups=T 
Ref_SpaceMapXGroup = SpaceMapX::run(Ref_SpaceMapX,
                                  cutoff=0.1,
                                  out_dir="./CloneFinalGroup", 
                                  num_threads = 10,
                                  denoise=T,
                                  output_format = "png",
                                  hclust_method='ward.D2',
                                  cluster_by_groups=T,
                                  analysis_mode = "samples",
                                  tumor_subcluster_partition_method = "qnorm",
                                  HMM=F)

```
# this is the final SPACEmapX results.

## Results Demo
![image](https://github.com/yintz/SPACEmapX/blob/main/Figures/P15.preliminary.png)

# Citation

# Funding 



