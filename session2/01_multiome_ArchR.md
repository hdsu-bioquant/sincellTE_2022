

# Multiome scRNA-seq/scATAC-seq analysis with ArchR

## Data download

In this tutorial we will use the scRNA-seq/scATAC-seq multiome example data provided by 10x Genomicsfor human PBMCs.

The data was downloaded using the following commands:

```
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
```

## Create ArchR Arrow file

The main input to create an ArchR project are Arrow files created from the raw alignments, this files can be created from the fragments files resulting from the cellranger-atac pipeline, or from a bam file.

Here we are going to create an arrow file from the fragments file: *pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz*

```r
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                    Load package and global settings                        ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## Setting default genome to Hg38.
library(ArchR)
addArchRGenome("hg38")

## Setting default number of Parallel threads to 5.
addArchRThreads(5)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                    Load package and global settings                        ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Get fragment file
inputFiles <- list.files("data", pattern="fragments.tsv.gz$", full.names=TRUE)
names(inputFiles) <- "PBMC_10k"

# Create Arrow file
createArrowFiles(inputFiles  = inputFiles, 
                 sampleNames = "PBMC_10k", 
                 QCDir       = "data/QualityControl",
                 logFile     = createLogFile(name = "createArrows", 
                                             logDir = "data/ArchRLogs"),
                 force       = TRUE)

```

<details>
<summary><b>Click for Answer</b></summary>

```
Using GeneAnnotation set by addArchRGenome(Hg38)!
Using GeneAnnotation set by addArchRGenome(Hg38)!
ArchR logging to : data/ArchRLogs/ArchR-createArrows-130626bf3d8f6-Date-2021-12-30_Time-21-13-35.log
If there is an issue, please report to github with logFile!
2021-12-30 21:13:35 : Batch Execution w/ safelapply!, 0 mins elapsed.
2021-12-30 21:13:35 : (PBMC_10k : 1 of 1) Reading In Fragments from inputFiles (readMethod = tabix), 0.001 mins elapsed.
2021-12-30 21:13:35 : (PBMC_10k : 1 of 1) Tabix Bed To Temporary File, 0.001 mins elapsed.
2021-12-30 21:15:50 : (PBMC_10k : 1 of 1) Successful creation of Temporary File, 2.253 mins elapsed.
2021-12-30 21:15:50 : (PBMC_10k : 1 of 1) Creating ArrowFile From Temporary File, 2.253 mins elapsed.
2021-12-30 21:17:02 : (PBMC_10k : 1 of 1) Successful creation of Arrow File, 3.449 mins elapsed.
2021-12-30 21:18:08 : (PBMC_10k : 1 of 1) CellStats : Number of Cells Pass Filter = 11582 , 4.544 mins elapsed.
2021-12-30 21:18:08 : (PBMC_10k : 1 of 1) CellStats : Median Frags = 13610 , 4.544 mins elapsed.
2021-12-30 21:18:08 : (PBMC_10k : 1 of 1) CellStats : Median TSS Enrichment = 13.6245 , 4.544 mins elapsed.
2021-12-30 21:18:12 : (PBMC_10k : 1 of 1) Adding Additional Feature Counts!, 4.617 mins elapsed.
2021-12-30 21:18:34 : (PBMC_10k : 1 of 1) Removing Fragments from Filtered Cells, 4.981 mins elapsed.
2021-12-30 21:18:34 : (PBMC_10k : 1 of 1) Creating Filtered Arrow File, 4.982 mins elapsed.
2021-12-30 21:19:34 : (PBMC_10k : 1 of 1) Finished Constructing Filtered Arrow File!, 5.99 mins elapsed.
2021-12-30 21:19:35 : (PBMC_10k : 1 of 1) Adding TileMatrix!, 5.991 mins elapsed.
2021-12-30 21:22:37 : (PBMC_10k : 1 of 1) Adding GeneScoreMatrix!, 9.03 mins elapsed.
2021-12-30 21:24:23 : (PBMC_10k : 1 of 1) Finished Creating Arrow File, 10.8 mins elapsed.
ArchR logging successful to : data/ArchRLogs/ArchR-createArrows-130626bf3d8f6-Date-2021-12-30_Time-21-13-35.log
[1] "PBMC_10k.arrow"
```

</details>


## Create ArchR project

An ArchR project is created from a list of previously computed Arrow files

```r
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          Create ArchR project                              ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
archrproj <- ArchRProject(ArrowFiles = "PBMC_10k.arrow", 
                     outputDirectory = "results/ArchROutput")
```

<details>
<summary><b>Click for Answer</b></summary>

```
Using GeneAnnotation set by addArchRGenome(Hg38)!
Using GeneAnnotation set by addArchRGenome(Hg38)!
Validating Arrows...
Getting SampleNames...
1 
Copying ArrowFiles to Ouptut Directory! If you want to save disk space set copyArrows = FALSE
1 
Getting Cell Metadata...
1 
Merging Cell Metadata...
Initializing ArchRProject...

                                                   / |
                                                 /    \
            .                                  /      |.
            \\\                              /        |.
              \\\                          /           `|.
                \\\                      /              |.
                  \                    /                |\
                  \\#####\           /                  ||
                ==###########>      /                   ||
                 \\##==......\    /                     ||
            ______ =       =|__ /__                     ||      \\\
        ,--' ,----`-,__ ___/'  --,-`-===================##========>
       \               '        ##_______ _____ ,--,__,=##,__   ///
        ,    __==    ___,-,__,--'#'  ==='      `-'    | ##,-/
        -,____,---'       \\####\\________________,--\\_##,/
           ___      .______        ______  __    __  .______      
          /   \     |   _  \      /      ||  |  |  | |   _  \     
         /  ^  \    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\  \   |      /     |  |     |   __   | |      /     
       /  _____  \  |  |\  \\___ |  `----.|  |  |  | |  |\  \\___.
      /__/     \__\ | _| `._____| \______||__|  |__| | _| `._____|
    
```

</details>


## Add gene expression to ArchR project

Once the ArchR project is created the next step is to add the gene expression data:

First we have to import the feature matrix from the 10x feature hdf5 file:

```r
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          Adding scRNA-seq data                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
scRNA <- import10xFeatureMatrix(
    input = "data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
    names = "PBMC_10k")
scRNA
```

<details>
<summary><b>Click for Answer</b></summary>

```
class: RangedSummarizedExperiment 
dim: 36578 11909 
metadata(0):
assays(1): counts
rownames(36578): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
rowData names(5): feature_type genome id interval name
colnames(11909): PBMC_10k#AAACAGCCAAGGAATC-1 PBMC_10k#AAACAGCCAATCCCTT-1 ... PBMC_10k#TTTGTTGGTTGGTTAG-1
  PBMC_10k#TTTGTTGGTTTGCAGA-1
colData names(0):
    
```

</details>


Overlap of cells between the scATAC-seq data and the scRNA-seq data:

```r
ggVennDiagram::ggVennDiagram(list(scATAC_cells = ArchR::getCellNames(archrproj),
                                  scRNA_cells  = colnames(scRNA)))
```

<details>
<summary><b>Click for Answer</b></summary>

<img src="figs/scATAC_scRNA_overlap.png" width="90%" />

</details>

Adding scRNA-seq data to the ArchR project:

```r
archrproj <- addGeneExpressionMatrix(input = archrproj, seRNA = scRNA, force = TRUE)
```

<details>
<summary><b>Click for Answer</b></summary>

```
ArchR logging to : ArchRLogs/ArchR-addGeneExpressionMatrix-462713e88a37d-Date-2021-12-31_Time-05-58-35.log
If there is an issue, please report to github with logFile!
2021-12-31 05:58:35 : Overlap w/ scATAC = 0.981
2021-12-31 05:58:35 : Overlap Per Sample w/ scATAC : PBMC_10k=11361
2021-12-31 05:58:38 : Batch Execution w/ safelapply!, 0 mins elapsed.
2021-12-31 05:58:40 : Adding PBMC_10k to GeneExpressionMatrix for Chr (1 of 23)!, 0.022 mins elapsed.
2021-12-31 05:58:42 : Adding PBMC_10k to GeneExpressionMatrix for Chr (2 of 23)!, 0.066 mins elapsed.
2021-12-31 05:58:45 : Adding PBMC_10k to GeneExpressionMatrix for Chr (3 of 23)!, 0.101 mins elapsed.
2021-12-31 05:58:47 : Adding PBMC_10k to GeneExpressionMatrix for Chr (4 of 23)!, 0.148 mins elapsed.
2021-12-31 05:58:49 : Adding PBMC_10k to GeneExpressionMatrix for Chr (5 of 23)!, 0.179 mins elapsed.
2021-12-31 05:58:51 : Adding PBMC_10k to GeneExpressionMatrix for Chr (6 of 23)!, 0.213 mins elapsed.
2021-12-31 05:58:54 : Adding PBMC_10k to GeneExpressionMatrix for Chr (7 of 23)!, 0.259 mins elapsed.
2021-12-31 05:58:56 : Adding PBMC_10k to GeneExpressionMatrix for Chr (8 of 23)!, 0.294 mins elapsed.
2021-12-31 05:58:58 : Adding PBMC_10k to GeneExpressionMatrix for Chr (9 of 23)!, 0.325 mins elapsed.
2021-12-31 05:59:01 : Adding PBMC_10k to GeneExpressionMatrix for Chr (10 of 23)!, 0.37 mins elapsed.
2021-12-31 05:59:03 : Adding PBMC_10k to GeneExpressionMatrix for Chr (11 of 23)!, 0.403 mins elapsed.
2021-12-31 05:59:05 : Adding PBMC_10k to GeneExpressionMatrix for Chr (12 of 23)!, 0.437 mins elapsed.
2021-12-31 05:59:07 : Adding PBMC_10k to GeneExpressionMatrix for Chr (13 of 23)!, 0.484 mins elapsed.
2021-12-31 05:59:09 : Adding PBMC_10k to GeneExpressionMatrix for Chr (14 of 23)!, 0.514 mins elapsed.
2021-12-31 05:59:11 : Adding PBMC_10k to GeneExpressionMatrix for Chr (15 of 23)!, 0.546 mins elapsed.
2021-12-31 05:59:14 : Adding PBMC_10k to GeneExpressionMatrix for Chr (16 of 23)!, 0.591 mins elapsed.
2021-12-31 05:59:16 : Adding PBMC_10k to GeneExpressionMatrix for Chr (17 of 23)!, 0.624 mins elapsed.
2021-12-31 05:59:18 : Adding PBMC_10k to GeneExpressionMatrix for Chr (18 of 23)!, 0.657 mins elapsed.
2021-12-31 05:59:20 : Adding PBMC_10k to GeneExpressionMatrix for Chr (19 of 23)!, 0.701 mins elapsed.
2021-12-31 05:59:22 : Adding PBMC_10k to GeneExpressionMatrix for Chr (20 of 23)!, 0.734 mins elapsed.
2021-12-31 05:59:24 : Adding PBMC_10k to GeneExpressionMatrix for Chr (21 of 23)!, 0.765 mins elapsed.
2021-12-31 05:59:27 : Adding PBMC_10k to GeneExpressionMatrix for Chr (22 of 23)!, 0.807 mins elapsed.
2021-12-31 05:59:29 : Adding PBMC_10k to GeneExpressionMatrix for Chr (23 of 23)!, 0.838 mins elapsed.
ArchR logging successful to : ArchRLogs/ArchR-addGeneExpressionMatrix-462713e88a37d-Date-2021-12-31_Time-05-58-35.log
```

</details>

Filtering out low quality cells and doublets:

```r
# Low quality cells
archrproj <- archrproj[archrproj$TSSEnrichment > 6 & archrproj$nFrags > 2500 & !is.na(archrproj$Gex_nUMI)]

# Filtering doublets
archrproj <- addDoubletScores(archrproj)
archrproj <- filterDoublets(archrproj)

archrproj
```

<details>
<summary><b>Click for Answer</b></summary>

```
ArchR logging to : ArchRLogs/ArchR-addDoubletScores-49f5a3d2e03a0-Date-2021-12-31_Time-06-21-52.log
If there is an issue, please report to github with logFile!
2021-12-31 06:21:52 : Batch Execution w/ safelapply!, 0 mins elapsed.
2021-12-31 06:21:52 : PBMC_10k (1 of 1) :  Computing Doublet Statistics, 0 mins elapsed.
PBMC_10k (1 of 1) : UMAP Projection R^2 = 0.9982
ArchR logging successful to : ArchRLogs/ArchR-addDoubletScores-49f5a3d2e03a0-Date-2021-12-31_Time-06-21-52.log

Filtering 1185 cells from ArchRProject!
	PBMC_10k : 1185 of 10887 (10.9%)



           ___      .______        ______  __    __  .______      
          /   \     |   _  \      /      ||  |  |  | |   _  \     
         /  ^  \    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\  \   |      /     |  |     |   __   | |      /     
       /  _____  \  |  |\  \\___ |  `----.|  |  |  | |  |\  \\___.
      /__/     \__\ | _| `._____| \______||__|  |__| | _| `._____|
    
class: ArchRProject 
outputDirectory: /home/bq_aquintero/projects/sincell_2022/results/ArchROutput 
samples(1): PBMC_10k
sampleColData names(1): ArrowFiles
cellColData names(19): Sample TSSEnrichment ... DoubletScore DoubletEnrichment
numberOfCells(1): 9702
medianTSS(1): 13.598
medianFrags(1): 13571
```

</details>

