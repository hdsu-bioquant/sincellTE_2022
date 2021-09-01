# [BC]<sup>2</sup> Tutorial - Defining genomic signatures with Non-Negative Matrix Factorization

![](./The_end.png)

Welcome to the **[BC]<sup>2</sup> Tutorial - Defining genomic signatures with Non-Negative Matrix Factorization** This tutorial will guide you through all the necessary steps to extract genomic signatures from high-dimensional data using Non-Negative Matrix Factorization (NMF). 

The tutorial will run over one complete day **(Monday, 13 September 2021)** from 9:00 am to 4:00 pm.

******
## Tutors

* Carl Herrmann, [Health Data Science Unit](https://www.hdsu.org/) Medical Faculty Heidelberg and BioQuant (carl.herrmann@bioquant.uni-heidelberg.de)
* Andres Quintero, [Health Data Science Unit](https://www.hdsu.org/) Medical Faculty Heidelberg and BioQuant (andres.quintero@bioquant.uni-heidelberg.de)


********

## Is this tutorial for me?

The aim of this tutorial is to learn how to use the R package ButchR to **perform signature identification in different types of genomic data using NMF**. To explore the results of an NMF analysis, we will provide a ready to use Docker image with RStudio, ButchR, and pre-loaded publicly available datasets, including bulk and single-cell RNA-seq data, as well as an interactive application. The tutorial will show how to run an NMF-based analysis from start to end.

If you are a computational biologist dealing with large scale omics datasets (e.g. RNA-seq, ATAC-seq, …) looking for solutions to reduce the dimensionality of the data to a small set of informative signatures, this tutorial will be perfect for you.


### IMPORTANT NOTE! 

 While we will start at a very basic level, we would **strongly encourage absolute beginners**, who have never ever worked with R, to complete a very simple online R intro course on DataCamp (["Introduction to R"](https://learn.datacamp.com/courses/free-introduction-to-r)), which will give you the very basic first concepts on what R is, and how to do some very simple operations with it.


********


## Schedule

<!-- <style>
.heatMap {
    width: 70%;
    text-align: center;
}
.heatMap th {
background: grey;
word-wrap: break-word;
text-align: center;
}
.heatMap tr:nth-child(1) { background: firebrick; }
.heatMap tr:nth-child(4) { background: cadetblue; }
.heatMap tr:nth-child(5) { background: firebrick; }
.heatMap tr:nth-child(9) { background: cadetblue; }
.heatMap tr:nth-child(10) { background: firebrick; }
.heatMap tr:nth-child(15) { background: firebrick; }
.heatMap tr:nth-child(17) { background: darkcyan; }

.heatMap th:first-of-type { width: 50%; }
.heatMap th:nth-of-type(2) { width: 20%; }

</style> -->

<div class="heatMap">

| Activity | Time |
| -- | ----------- |
| Session 1 - Introduction |  |
| Ice breaker: Course expectations | 9:00 - 9:30 | 
| Introduction to Non-Negative Matrix Factorization (NMF) and its usage in genomics | 9:30 - 10:15 | 
| Cofee break and discussion | 10:15 - 10:45| 
| Session 2 - Matrix decomposition |  |
| How to use ButchR with Docker | 10:45 - 11:15 | 
| Pre-processing data to use with NMF | 11:15 - 11:45 | 
| Matrix decomposition with ButchR | 11:45 - 12:15 | 
| Lunch break | 12:15 - 13:30 | 
| Session 3 - Result interpretation |  |
| Selection of optimal factorizatin rank | 13:30 - 14:00 | 
| Signature identification | 14:00 - 14:30 | 
| Feature extraction and enrichment analysis | 14:30 - 15:00 | 
| Interactive analysis with ShinyButchR | 15:00 - 15:30| 
| Session 4 - Discussion | | 
| Discussion and concluding remarks | 15:30 - 16:00| 
| [BC]<sup>2</sup> Welcome lecture | 17:00 | 

</div>




********
<!-- ## Slides

Here are the links to the slides

* Day 1 : [introduction](./irtg2021_intro.pdf)
* Day 1 : [R markdown](./irtg2021_rmarkdown.pdf)
* Day 1 : [Data types](./irtg2021_datatypes.pdf)
* Day 1 : [Statistical tests](./irtg2021_tests.pdf)

* Day 2 : [Introduction to single-cell analysis](https://docs.google.com/presentation/d/1DSC6gUIbO6PzrqLCt1jp-sIx1U31TvMdDGgKdhohCIY/edit?ts=60c8bafb#slide=id.gdf238a40cf_0_5) -->


## Practical parts


**You need to install a couple of R packages; you can download the 
[following script](./install_packages.R). Load it into RStudio, and then hit the *Run* button at the top to execute it. It should run smoothly!**

**Please document your progress in this [Google Sheet](https://docs.google.com/spreadsheets/d/1rFcWJJD-qOqeRWZvhqPEqMCt_ddtinvdTlLPl2Syomw/edit?usp=sharing)**

### Day 1: General introduction - (almost) first steps in R!                                        


On the first day, we will guide you through the first steps of working with R, from reading data to exploratory analysis and basic statistics.

0. [Goals of Day 1](./day1/00_Objectives.md)
1. [Getting started with RStudio](./day1/01_rstudio.md)
2. [Reading in a data table](./day1/02_dataframe.md)
3. [Cleaning the dataset](./day1/03_cleanup.md)
4. [Making plots](./day1/04_plotting.md)
5. [Statistical tests](./day1/05_test.md)
                                                   

### Day 2: a simple single-cell RNA-seq analysis workflow

On the second day, we will go through a step by step simple analysis of a small scRNA-seq dataset using the Seurat toolkit. **Don't expect to be able to carry a full scRNA-seq analysis after this!** This is meant to give you an idea of a typical workflow rather.

0. [Intro on scRNA-seq analysis](./day2/00.md)
1. [Creating a SEURAT object](./day2/01.md)
2. [Basic QC and normalization ](./day2/02.md)
3. [Feature selection](./day2/03.md)
4. [Normalization and Dimensionality reduction](./day2/04.md)
5. [Cluster visualization](./day2/05.md)
6. [Differential gene expression analysis](./day2/06.md)
7. [Profiling cells](./day2/07.md)




*********
## Organisation

The course will be online only! 
* Lectures will be over **Zoom** (we will send the link via email prior to the course)
* We will use **Discord channels** for the practical sessions (register [using this link](https://discord.gg/gPXJDukGfQ))

**********
## Technical pre-requisites

 The attendees are expected to bring their own laptop with Docker pre-installed. To avoid any delay in setting up the container during the practice sessions, the Docker image for the workshop should be downloaded beforehand. This can be done by opening a command-line terminal (e.g., Powershell and Terminal) and running the command “docker pull hdsu/butchr”. A complete overview of how to install Docker can be found here: https://docs.docker.com/desktop/. 


Please check that you can run the hdsu/butchr image without error message!