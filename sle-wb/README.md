

## A systemic lupus erythematosus whole blood compendium

Systemic lupus erythematosus (SLE) is an autoimmune disease that affects multiple organs. Unlike some other conditions that fall under the rheumatic/autoimmune disease umbrella, SLE is not a rare disease and there are multiple gene expression studies (and/or cohorts) of SLE that are publicly available.

Herein, we choose to study SLE whole blood (WB) samples because methods around analyzing these data, including [modular transcriptional analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4118927/), are mature. This tissue—a mixture of adaptive and innate immune cells—is also highly relevant to this condition and allows us to explore whether or not we can capture or retain cell type-specific expression patterns in our processed data.

In this repository, we compare several normalization and processing methods for the purpose of integrating 7 microarray datasets into a "harmonized" compendium for further analyses. 

#### Pipeline

To follow our SLE WB pipeline, run the numbered R scripts in this directory (`sle-wb`) from the top directory of this repository sequentially. For example, the first step of our pipeline can be run from the top directory with the following command:

```shell
Rscript sle-wb/1-process_affy_data.R
```

### Datasets

We include the following datasets (retrieved from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/)):

| Accession                                |     # Samples      |          Platform           | Raw Available? | Citation (if available)                  |
| ---------------------------------------- | :----------------: | :-------------------------: | :------------: | :--------------------------------------- |
| [E-GEOD-65391](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-65391/) |        996         |  Illumina HumanHT-12 V4.0   |       N        | [Banchereau, et al. _Cell._ 2016.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5426482/) |
| [E-GEOD-11907](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-11907/) | 41 (SLE & Control) |   Affymetrix HG-U133A & B   |       Y        | [Chaussabel, et al. _Immunity._ 2008.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2727981/) |
| [E-GEOD-39088](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-39088/) |        142         | Affymetrix HG-U133 plus 2.0 |       Y        | [Lauwerys, et al. _Arthritis Rheum._ 2013.](https://doi.org/10.1002/art.37785) |
| [E-GEOD-72747](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-72747/) |         30         | Affymetrix HG-U133 plus 2.0 |       Y        | [Ducreux, et al. _Rheumatology (Oxford)._ 2016.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5034220/) |
| [E-GEOD-61635](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-61635/) |        129         | Affymetrix HG-U133 plus 2.0 |       Y        |                                          |
| [E-GEOD-49454](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-49454/) |        177         |  Illumina HumanHT-12 V4.0   |       N        | [Chiche, et al. _Arthritis Rheumatol._ 2014.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4157826/) |
| [E-GEOD-78193](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-78193/) |        125         |        Agilent 4x44K        |       N        | [Welcher, et al. _Arthritis Rheumatol._ 2015.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5054935/) |

These selected data sets are somewhat representative of the larger body of publicly available in a few ways: 1) the majority of assays are run on Affymetrix (for which raw data are typically available) and 2) some accessions only have submitter-processed data available on ArrayExpress (like E-GEOD-78193).

### Correlation between normalization methods 

First, we compare methods for normalizing Affymetrix arrays: [Robust Multiarray Average (RMA)](https://doi.org/10.1093/biostatistics/4.2.249), [Single Channel Array Normalization (SCAN)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3508193/) and its "fast" counterpart (SCANfast), and find that they are generally well correlated. (See [`sle-wb/plots/norm_method_correlation`](https://github.com/greenelab/rheum-plier-data/tree/master/sle-wb/plots/norm_method_correlation) for the full output of this analysis.)

![](https://github.com/greenelab/rheum-plier-data/raw/master/sle-wb/plots/norm_method_correlation/SLE-WB_affy_norm_correlation_hgu133plus2_RMA_v_SCAN.png)

We're interested in normalizing single samples (SCAN/SCANfast), rather than batches of multiple samples (RMA), because that removes any influence of the other samples processed jointly on the single sample's output. This is particularly attractive when thinking about assembling [a large, consistently updated compendium](http://www.ccdatalab.org/blog/data-refinery-one/) of gene expression samples.

### Platform and/or batch effect 

