# A data repository for rheumatic disease gene expression data

This repository contains data and processing code for use in a project examining gene expression patterns in autoimmune/rheumatic diseases. 

## Datasets

### A note on processing 

The processing code for each dataset (or compendium in the case of [`sle-wb`](https://github.com/greenelab/rheum-plier-data/tree/master/sle-wb)) is contained within each subdirectory (if applicable). 
For more information on our data processing strategy, see [`sle-wb/README.md`](https://github.com/greenelab/rheum-plier-data/blob/master/sle-wb/README.md).

### [recount2](https://jhubiostatistics.shinyapps.io/recount/)

Within this repository, we obtain recount2 data through the [`recount` bioconductor package](http://bioconductor.org/packages/release/bioc/html/recount.html), further process it, and apply [`PLIER`](https://github.com/wgmao/PLIER).

The recount2 data and results are too large to be stored with Git LFS, so we have placed them on [figshare](https://figshare.com/). 
1. [Output](https://doi.org/10.6084/m9.figshare.5716033.v1) from [`recount2/get.all.recount.dataset.R`](https://github.com/greenelab/rheum-plier-data/blob/ee962291521bbe88ae8b92fddddfb9e7ab3ee667/recount2/get.all.recount.dataset.R)
2. [Output](https://doi.org/10.6084/m9.figshare.5712124.v3) from [`recount2/recount.plier.R`](https://github.com/greenelab/rheum-plier-data/blob/ee962291521bbe88ae8b92fddddfb9e7ab3ee667/recount2/recount.plier.R)

Citations:

> Collado-Torres L, Nellore A, Kammers K, et al. [Reproducible RNA-seq analysis using _recount2_](https://http://doi.org/10.1038/nbt.3838). _Nature Biotechnology_, 2017. doi: 10.1038/nbt.3838. 

> Mao W, Chikina M. [Pathway-Level Information ExtractoR (PLIER): a generative model for gene expression data](http://dx.doi.org/10.1101/116061). _bioRxiv_, 2017. doi: 10.1101/116061

### Granulomatosis with polyangiitis

Two GPA (Wegener's) datasets are included in this repository:

* [NARES](https://github.com/greenelab/rheum-plier-data/tree/master/NARES) -- a dataset that consists of nasal brushings from patients with GPA with or without a history of nasal disease.
* [GSE18885](https://github.com/greenelab/rheum-plier-data/tree/master/gpa-blood) -- a blood (fractions) dataset; we use submitter-processed data from [GEO](https://www.ncbi.nlm.nih.gov/geo/).

Citations:

> Grayson PC, Steiling K, Platt M, et al. [Defining the Nasal Transcriptome in Granulomatosis with Polyangiitis](https://dx.doi.org/10.100art.39185). _Arthritis & Rheumatology_, 2015. doi: 10.1002/art.39185.

> Cheadle C, Berger AE, Andrade F, et al. [Transcription of PR3 and Related Myelopoiesis Genes in Peripheral Blood Mononuclear Cells in Active Wegener’s Granulomatosis](https://dx.doi.org/10.1002/art.27398). _Arthritis & Rheumatism_, 2010. doi: 10.1002/art.27398.

### Systemic lupus erythematosus whole blood

See [`sle-wb`](https://github.com/greenelab/rheum-plier-data/tree/master/sle-wb) for more information (including citations).

## License 

This repository is dual licensed as [BSD 3-Clause](https://github.com/greenelab/rheum-plier-data/blob/master/LICENSE_BSD-3.md) (source code) and [CC0 1.0](https://github.com/greenelab/rheum-plier-data/blob/master/LICENSE_CC0.md) (figures, documentation, and our arrangement of the facts contained in the underlying data).
