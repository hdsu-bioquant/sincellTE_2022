# IRTG Course - Introduction to R for genomics

![circos](./circos.png)

Welcome to the **IRTG Course - Introduction to R** This workshop is meant for individuals with little previous knowledge R. 

The course will run over 2 days **(Wednesday, 16.06 and Thursday, 17.06)** from 10 am - 12.30 pm and 1.30 pm - 5.30 pm.


******
## Tutors

* Carl Herrmann, [Health Data Science Unit](https://www.hdsu.org/) Medical Faculty Heidelberg and BioQuant (carl.herrmann@bioquant.uni-heidelberg.de)
* Carlos Ramirez, [Health Data Science Unit](https://www.hdsu.org/) Medical Faculty Heidelberg and BioQuant (carlos.ramirez@bioquant.uni-heidelberg.de)


********

## Is this course for me?

In this two-day course, we want to give you an **introduction to working with R in simple data analysis tasks**; you will learn the basic principles of reading in a data table, doing some descriptive statistics, making nice plots.
On the second day, we will focus on a simple single-cell analysis workflow, which will guide you through the first steps of this kind of analysi!

### IMPORTANT NOTE! 

 While we will start at a very basic level, we would **strongly encourage absolute beginners**, who have never ever worked with R, to complete a very simple online R intro course on DataCamp (["Introduction to R"](https://learn.datacamp.com/courses/free-introduction-to-r)), which will give you the very basic first concepts on what R is, and how to do some very simple operations with it.
 We will send you a link so that you can freely register to DataCamp and follow this course. 

********
## Practical parts

### Day 1: General introduction - (almost) first steps in R!                                        


**You need to install a couple of R packages; you can download the 
[following script](./install_packages.R). Load it into RStudio, and then hit the *Run* button at the top to execute it. It should run smoothly!**


On the first day, we will guide you through the first steps of working with R, from reading data to exploratory analysis and basic statistics.

1. [Getting started with RStudio](./01_rstudio.md)
2. [Reading in a data table](./02_dataframe.md)
3. [Cleaning the dataset](./03_cleanup.md)
4. [Making plots](./04_plotting.md)
5. [Statistical tests](./05_test.md)
                                                   

### Day 2: a simple single-cell RNA-seq analysis workflow

On the second day, we will go through a step by step simple analysis of a small scRNA-seq dataset using the Seurat toolkit.

1. Reading the data
2. Creating a SEURAT object
3. Basic QC and normalization 
4. Visualization (UMAP) and cell clustering
5. Differential expression analysis




*********
## Organisation

The course will be online only! 
* Lectures will be over **Zoom** (we will send the link via email prior to the course)
* We will use **Discord channels** for the practical sessions (register [using this link](https://discord.gg/gPXJDukGfQ))

**********
## Technical pre-requisites

Every participant will work on her/his own laptop. The easiest way to work with R is using the **RStudio** interface.
Please install RStudio Desktop prior to the start of the course:

1. first install R for your operating system; you will find the correct version [on this website](https://cran.rstudio.com/) 
2. once R is installed, you can isntall the RStudio Desktop version, which you find [here](https://www.rstudio.com/products/rstudio/download/#download)

Please check that you can open RStudio without error message!