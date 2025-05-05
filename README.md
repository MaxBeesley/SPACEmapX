<img src="https://www.nds.ox.ac.uk/images/logos/secondary-logo" height="75" /> <img src="https://www.nds.ox.ac.uk/images/logos/primary-logo" height="75"/> 
<img src="https://github.com/user-attachments/assets/3615186f-6b02-4df3-a012-afe77714044a" height="75" />


<a href="https://twitter.com/intent/follow?screen_name=Wien_Yin">
<img src="https://img.shields.io/twitter/follow/Wien_Yin?style=social&logo=X",alt="follow on Twitter"></a>
<a href="https://twitter.com/intent/follow?screen_name=lambalastair">
<img src="https://img.shields.io/twitter/follow/lambalastair?style=social&logo=X",alt="follow on Twitter"></a>
<a href="https://twitter.com/intent/follow?screen_name=OxPCaBiol">
<img src="https://img.shields.io/twitter/follow/OxPCaBiol?style=social&logo=X",alt="follow on Twitter"></a>


[![](https://img.shields.io/badge/SPACEmapX-version0.99-blue.svg)](https://github.com/yintz/SPACEmapX/releases)[![](https://img.shields.io/github/last-commit/yintz/SPACEmapX.svg)](https://github.com/yintz/SPACEmapX/commits/main)


# SPACEMapX - Spatial Phylogenetic Analysis and Clonal Evolution: MAPping the lethal X clone

This package provides guidance for analysis of spatially resolved tissue (Visium ST), focusing on clonal dynamics derived from inferred copy number status.

Our particular goal has been identification of the "lethal clone", defined as the clone which metasises from primary cancer tissue to lymph nodes. However, this package can be used for any analysis of clonal evolution and/or somatic mosaicism in heterogeneous tissue.  

If you need any assistance with running this package, please feel free to contact us via GitHub "issue" messaging



## Purpose
SPACEmapX is based on inferCNV and done in collaboration with the authors of SpatialInferCNV. This package is designed to simplify complex workflows when using inferCNV on multiple sections of Spatial Transcriptomic (ST) data from an individual patient.




## Overview

This packages is a wrapper to provide I) support for importing and processing spatial transriptomics data, II)  improved Heatmap Plots and III) a tool to select dendrogram barcodes. These are additions to the core inferCNV package, which is used to infer Copy Number Variation (CNV) from transcriptomic data and provides the design of the objects used throughout. 
inferCNV: https://github.com/broadinstitute/infercnv






# Software Requirements
JAGS
R (tutorial completed on R4.4.1)
R libraries: inferCNV, Seurat, SPACEmapX

## JAGS Pre-requirement
JAGS
In order to run infercnv, JAGS (Just Another Gibbs Sampler) must be installed.

Install JAGS
Install the JAGS package. This can be installed for Mac, Windows, or Linux. Be sure to download and install the corresponding package matching your OS.


# Installation
The package can be installed through GitHub using;
``` r
install.packages("remotes")
remotes::install_github("yintz/SPACEmapX")
```


# Tutorial 
https://github.com/MaxBeesley/SPACEmapX/wiki/SPACEmapX-tutorial





## Novel functions introduced in addition to SpatialInferCNV
``` r
LoadSTinfo(summary_csv_file)
LoadSTdata(loadedSTinfoData)
ShowTwigSpotsList(TwigID1,TwigID2,TwigID3,...)
ShowTwig(TwigID)
ShowTwigSectionName(TwigID)
```





# Data

This package has been designed to work with Visium Fresh Frozen, FFPE V1,V2 data.

If you would like a representative ST test dataset, please try these data [from Erickson et al, Nature 2022; https://www.nature.com/articles/s41586-022-05023-2]. 
[MendeleyRepository](https://data.mendeley.com/v1/datasets/svw96g68dv/draft?a=3f263217-2bd3-4a3c-8125-8c517c3a9e29).

For more ST data, you can also go to 10X website to download assorted ST datasets [https://www.10xgenomics.com/datasets?query=&page=1&configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&refinementList%5Bproduct.name%5D%5B0%5D=Spatial%20Gene%20Expression&refinementList%5Bspecies%5D%5B0%5D=Human&refinementList%5Bchemistry.version%5D%5B0%5D=2].


# Funding 
Cancer Research UK (CRUK) #C57899/A25812 "Spatial Prostate Assessment and Circulating Environment – The SPACE Study"
