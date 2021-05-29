# General instructions for the practical session

## 1. bash console / R console: where am I?


When you connect to the VM using the ssh command, you are in the **bash console**. You can recognize that you are in the bash console by looking at the **prompt** (i.e. the left part of the input line); it should look line this:

```bash
user5@powerfulrutherford-30fb1:~$
```

Each participant has a home directory on the virtual machine, for example:

```bash
/home/user5/
```

Each user can jump directly into his home directory using the simple command

```bash
cd
````



Sometimes, we will use **R** to process the data; to activate R, simply type the following into the bash console:

```bash
/usr/bin/R
```

Now, you will be in the **R console**, which you can recognize by looking at the prompt, which now looks like this:

```
> 
```

To get back to the **bash console**, simply type the following in the R console (and answer `n` to any question):

```r
q()
```

You should be back to the **bash console** (check the prompt)!

## 2. Do I really need to type all the commands?

In the tutorial, we have written the commands in the grey boxes in such a way that you can copy/paste the command into your bash console, and execute it. If you do so, make sure that you understand how the command is structured and what it does!

You can also type the command yourself in the console; it that case, make sure to respect the blank spaces inside the command!

## 3. Where is the data?

Each user has in his home directory 2 folders:
* `data`: this folder contains all the data (fastq / bam / ...) that you will need for the analysis
* `analysis`: this folder will be used to store all the outputs of your analysis.

The `data` folder has the following structure:

```
data
├── ext_data
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── hg38.genome
│   ├── MA0139.1.jaspar
│   └── motifs.jaspar
├── fastqdata
│   ├── ATACseq
│   └── ChIPseq
└── processed
    ├── ATACseq
    ├── CTCF
    └── H3K4me3
```

* `ext_data`: contains files  needed for the analysis
* `fastqdata`: contains some of the raw and trimmed fastq files
* `processed`: contains some of the pre-processed files (like the aligned bam files)

During the analysis, we will create further subdirectories into the `analysis` directory.

**+++ GOOD PRACTICE ADVICE++**
 
1. structure you analysis directory such that you can easily trace back what type of file/results is stored in which folder!
2. often, some input files are processed and generate new output files. Some tools give automatically names to the output file, otherwise it is a good practice to give names that make it clear what the file contains...
Example:
```
## original fastq file
CTCF_rep1_IP.fq.gz

## aligned file

### bad name
IP.bam 

### good name
CTCF_rep1_IP.bam

## after filtering

### bad name
filtered_file.bam

### good name
CTCF_rep1_IP.mapq_filtered.dup_removed.bam
```

## 4. How do I access the tools needed for the analysis?

We will use a number of software tools for the analysis; we have prepared a **virtual environment** using **conda** containing all required tools. You first need to activate this environment, by running the simple command (in the bash console):

```bash
conda activate chipatac
```

You should see the prompt changing from 

```
user5@powerfulrutherford-30fb1:~$
```
to
```
(chipatac) user5@powerfulrutherford-30fb1:~$
```

Now you can use all tools described in the tutorial!

## 5. Why is it so slow??

You are working on a server with 26 CPU cores and 512 Gb memory. Some steps of the analysis are quite computationally intensive, and can lead to delays, especially if run simultaneously by several users! In that case, take a break while it is running, and get a fresh cup of coffee, you might need it ;-)