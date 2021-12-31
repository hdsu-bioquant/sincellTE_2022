

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
