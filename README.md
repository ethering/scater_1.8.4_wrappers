# Wrappers for Scater v1.8.4

This code wraps a number of [scater](https://bioconductor.org/packages/release/bioc/html/scater.html) functions, ultimately for their use in Galaxy wrappers. Each script can be used in isolation. Briefly, the `scater-create-qcmetric-ready-sce.R` script takes a sample x gene expression matrix (usually read-counts) and cell annotation file, creates a [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object and runs scater's calculateQCMetrics function (using other supplied files such as ERCC's and mitochondrial gene features).
Various filter scripts are provided, along with some plotting functions for QC.





## Install
The scripts can be used direclty from the command line, althought the following R packages are required.

```
>source("https://bioconductor.org/biocLite.R")
Bioconductor version 3.7 (BiocInstaller 1.30.0), ?biocLite for help
>biocLite("scater")           #Should install >= scater-1.8.4
>devtools::install_github("ebi-gene-expression-group/workflow-scripts-common")
>install.packages(c("ggpubr", "optparse", "mvoutlier"))

```
R-3.5.1 was used to develop the code. Earlier versions have not been tested.


## Commands

For help with any of the following scripts, run:
 `<script-name> --help`

---

`scater-create-qcmetric-ready-sce.R`
Takes an expression matrix (usually read-counts) of samples (columns) and gene/transcript features (rows), along with other annotation information, such as cell metadata, control genes (mitochondrail genes, ERCC's), creates a [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) object and runs scater's calculateQCMetrics. Returns a SingleCellExperiment RDS object.


```
./scater-create-qcmetric-ready-sce.R -a test-data/counts.txt -c test-data/annotation.txt -f test-data/mt_controls.txt  -o test-data/scater_qcready_sce.rds
```
---

`scater-plot-dist-scatter.R`
Takes SingleCellExperiment object (from RDS file) and lots a panel of read and feature graphs, including the distribution of library sizes, distribution of feature counts, a scatterplot of reads vs features, and % of mitochondrial genes in library.

```
./scater-plot-dist-scatter.R -i test-data/scater_qcready_sce.rds -o test-data/scater_reads_genes_dist.pdf
```

![scater-plot-dist-scatter](images/scater_reads_genes_dist.png?raw=true)
---
`scater-plot-exprs-freq.R`
Plots mean expression vs % of expressing cells and provides information as to the number of genes expressed in 50% and 25% of cells.


![scater-plot-exprs-freq](images/scater_exprs_freq.png?raw=true)
---

`scater-pca-filter.R`
Takes SingleCellExperiment object (from RDS file) and automatically removes outliers from data using PCA. Returns a filtered SingleCellExperiment RDS object.

```
./scater-pca-filter.R -i test-data/scater_qcready_sce.rds -o test-data/scater_pca_filtered.rds
```

___

`scater-manual-filter.R`
Takes SingleCellExperiment object (from RDS file) and filters data using user-provided parameters. Returns a filtered SingleCellExperiment RDS object.

```
./scater-manual-filter.R -i test-data/scater_qcready_sce.rds -l 10000 -d 4 -m 33 -o test-data/scater_manual_filtered.rds
```
---
`scater-normalize.R`
Compute log-normalized expression values from count data in a SingleCellExperiment object, using the size factors stored in the object. Returns a normalised SingleCellExperiment RDS object.

```
./scater-normalize.R -i test-data/scater_manual_filtered.rds -o test-data/scater_man_filtered_normalised.rds
```
---
`scater-plot-pca.R`
PCA plot of a normalised SingleCellExperiment object (use `scater-normalize.R` before use). The options `-c`, `-p`, and `-s` all refer to cell annotation features. These are the column headers of the `-c` option used in `scater-create-qcmetric-ready-sce.R`.

```
./scater-plot-pca.R -i test-data/scater_man_filtered_normalised.rds -c Treatment -p Mutation_Status -o test-data/scater_pca_plot.pdf
```

## Typical workflow

#### 1. Read in data

```
./scater-create-qcmetric-ready-sce.R -a test-data/counts.txt -c test-data/annotation.txt -f test-data/mt_controls.txt  -o test-data/scater_qcready_sce.rds
```

#### 2. Visualise it
Take a look at the distribution of library sizes, expressed features and mitochondrial genes.
```
./scater-plot-dist-scatter.R -i test-data/scater_qcready_sce.rds -o test-data/scater_reads_genes_dist.pdf
```
Then look at the distibution of genes across cells

```
./scater-plot-exprs-freq.R -i test-data/scater_qcready_sce.rds -o test-data/scater_exprs_freq.pdf
```


#### 3. Use the plots to decide on filtering parameters
Option 1.
Automatically filter on outliers

```
./scater-pca-filter.R -i test-data/scater_qcready_sce.rds -o test-data/scater_filtered.rds
```

Option 2.
Manually filter

```
./scater-manual-filter.R -i test-data/scater_qcready_sce.rds -l 10000 -d 4 -m 33 -o test-data/scater_filtered.rds
```


#### 4. Visualise data again to see how the filtering performed

```
./scater-plot-dist-scatter.R -i test-data/scater_filtered.rds -o test-data/scater_filtered_genes_dist.pdf
```
Decide if you're happy with the data. If not, try increasing or decreasing the filtering parameters.


#### 5. Normalise data

```
./scater-normalize.R -i test-data/scater_filtered.rds -o test-data/scater_filtered_normalised.rds
```

#### 6. Investigate other confounding factors
Plot the data (using PCA) and display various annotated properties of the cells

```
./scater-plot-pca.R -i test-data/scater_filtered_normalised.rds -c Treatment -p Mutation_Status -o test-data/scater_pca_plot.pdf
```
