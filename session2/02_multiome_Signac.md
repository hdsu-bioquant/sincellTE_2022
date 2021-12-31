

# Multiome scRNA-seq/scATAC-seq analysis with Signac

## Data download

In this tutorial we will use the scRNA-seq/scATAC-seq multiome example data provided by 10x Genomics for human PBMCs.

The data was downloaded using the following commands:

```
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
```

## Load the gene expression and chromatin accessibility data


```r
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                              Load packages                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                              Load data                                     ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# load the RNA and ATAC data
counts <- Read10X_h5("data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
# If the previous line fails, please use the following:
# counts <- readRDS("data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.RDS")

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            Load annotation                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
# If the previous line fails, please use the following:
# ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
# seqlevels(annotation) <- ucsc.levels

# Check annotation seqnames
unique(seqnames(annotation))

# Check counts
lapply(counts, dim)
lapply(counts, class)

```

<details>
<summary><b>Click for Answer</b></summary>

```
> unique(seqnames(annotation))
 [1] chrX  chr20 chr1  chr6  chr3  chr7  chr12 chr11 chr4  chr17 chr2  chr16 chr8  chr19 chr9  chr13 chr14 chr5  chr22 chr10
[21] chrY  chr18 chr15 chr21 chrM 
25 Levels: chrX chr20 chr1 chr6 chr3 chr7 chr12 chr11 chr4 chr17 chr2 chr16 chr8 chr19 chr9 chr13 chr14 chr5 ... chrM

> lapply(counts, dim)
$`Gene Expression`
[1] 36601 11909

$Peaks
[1] 108377  11909

> lapply(counts, class)
$`Gene Expression`
[1] "dgCMatrix"
attr(,"package")
[1] "Matrix"

$Peaks
[1] "dgCMatrix"
attr(,"package")
[1] "Matrix"


```

</details>




## Create a Signac object from the counts

To create a Signac object, the first step is to create a standard Seurat object including only the gene expression data. Once this object is created, a Chromatin Assay object is added as one of the Seurat object assays.

Once both datasets are loaded in the object, we can perform a QC analysis to remove low quality cells:

```r
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                              Create Signac                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# First create a Seurat object containing the gene expression data alone
signacobj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# Then, create a ChromatinAssay and add it to the previous object
signacobj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = "data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
  annotation = annotation
)

DefaultAssay(signacobj) <- "ATAC"

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                                    QC                                      ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
signacobj <- NucleosomeSignal(signacobj)
signacobj <- TSSEnrichment(signacobj)

VlnPlot(
  object = signacobj,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

```

<details>
<summary><b>Click for Answer</b></summary>

<img src="figs/signac_QC.png" width="90%" />

</details>


Low quality cells filtering:

```r
signacobj <- subset(
  x = signacobj,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
signacobj

```

<details>
<summary><b>Click for Answer</b></summary>

```
An object of class Seurat 
144978 features across 11331 samples within 2 assays 
Active assay: ATAC (108377 features, 0 variable features)
 1 other assay present: RNA

```

</details>

## Call peaks and add peak counts matrix


```r
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                    Call peaks and make feature matrix                      ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Call peaks using MACS2
peaks <- CallPeaks(signacobj, macs2.path = "/shared/software/miniconda/envs/macs2-2.2.7.1/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(signacobj),
  features = peaks,
  cells = colnames(signacobj)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
signacobj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = "data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
  annotation = annotation
)

peaks
```

<details>
<summary><b>Click for Answer</b></summary>

```
GRanges object with 131364 ranges and 6 metadata columns:
           seqnames              ranges strand |                   name     score fold_change neg_log10pvalue_summit
              <Rle>           <IRanges>  <Rle> |            <character> <integer>   <numeric>              <numeric>
       [1]     chr1         10032-10322      * |  SeuratProject_peak_44       142     4.99234                16.3379
       [2]     chr1       180709-181030      * |  SeuratProject_peak_45       149     5.12372                17.1108
       [3]     chr1       181296-181600      * |  SeuratProject_peak_46       291     7.35714                31.7325
       [4]     chr1       191304-191914      * |  SeuratProject_peak_47       142     4.99234                16.3379
       [5]     chr1       267874-268087      * |  SeuratProject_peak_48       134     4.86097                15.5761
       ...      ...                 ...    ... .                    ...       ...         ...                    ...
  [131360]     chrX 155880631-155881911      * | SeuratProject_peak_1..       824     8.67288                87.4204
  [131361]     chrX 155891339-155891781      * | SeuratProject_peak_1..       105     4.30809                12.5625
  [131362]     chrX 155966929-155967163      * | SeuratProject_peak_1..       134     4.86097                15.5761
  [131363]     chrX 155997247-155997787      * | SeuratProject_peak_1..       263     6.17155                28.7770
  [131364]     chrX 156029849-156030260      * | SeuratProject_peak_1..       106     4.33546                12.6467
           neg_log10qvalue_summit relative_summit_position
                        <numeric>                <integer>
       [1]                14.2177                      126
       [2]                14.9684                      124
       [3]                29.1857                      137
       [4]                14.2177                      145
       [5]                13.4782                      134
       ...                    ...                      ...
  [131360]                82.4206                      637
  [131361]                10.5583                      258
  [131362]                13.4782                      112
  [131363]                26.3134                      342
  [131364]                10.6377                      245
  -------
  seqinfo: 24 sequences from an unspecified genome; no seqlengths
```

</details>
