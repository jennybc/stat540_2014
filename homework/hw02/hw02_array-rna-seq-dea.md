Homework 02
======================================================================

In this assignment you will work with data originally published in "A comprehensive comparison of RNA-Seq-based transcriptome analysis from reads to differential gene expression and cross-comparison with microarrays: a case study in _Saccharomyces cerevisiae_." by Nookaew et al (Nucleic Acids Res. 2012 Nov 1;40(20):10084-97. PMID 22965124). The article is available here: [doi: 10.1093/nar/gks804](http://dx.doi.org/10.1093%2Fnar%2Fgks804).

The authors used two different platforms -- microarrays and RNA-Seq -- to obtain gene expression data for yeast grown in two conditions: batch medium and chemostat. They then compared the results of differential expression analysis (DEA) across platforms and under various combinations of aligners and methods for identifying differential expression. We will do the same. Some of the DEA methods covered in this class were included in their analysis, `edgeR` and `DESeq`, while others were not, such as `limma` and, in particular, the `voom` function.

We will work with their data obtained from the NCBI SRA and GEO repositories. Because of some differences in pre-processing your results will differ from theirs.

## Evaluation

The assignment has 25 points and represents 25% of your overall course grade. Overall presentation and mechanics is worth 5 points, with full points awarded if it's exceptional. The remaining 20 points are spread among the specific questions, as detailed below. We will try to give partial credit when/if we can figure out what you're trying to do.

Each student must submit their own work: no group efforts. It's fine to talk to fellow students and have discussion on the Google Group. If someone is really helpful, give them some credit. You have definitely crossed the line if there's any copying and pasting of code or wording.

## Late policy

Late submission will cost 20% (5 points) per day or fraction thereof.

If you expect that you will be late on the deadline, contact us before missing it, and submit whatever partial work you have. 

## How to submit the work

+ We will soon provide specific guidelines on how to submit homework 2.
+ In the meantime, write your homework in R Markdown (.rmd)
+ Once you finish your homework, you will add the R Markdown, Markdown and HTML to your private GitHub as you did for homework 1.




## Your mission

## Q1) Microarray Analysis

The six samples in this study, 3 replicates for each of the two conditions, were analysed on the Affymetrix Yeast Genome Array 2.0 platform. We have already downloaded the raw CEL files from GEO and normalised them. The normalized data is saved in the file [`GSE37599-data.tsv`](../../examples/yeastPlatforms/data/GSE37599-data.tsv).

### a) (1pt) Load Microarray Data

Load the normalized data.   
  
What are dimensions of the dataset? In addition to reporting number of rows and columns, make it clear what rows and columns represent and how you're interpreting column names.




### b) (2pt) Identify Sample Swap

The labels on two of the samples have been swapped, that is one of the batch samples has been labelled as chemostat and vice-versa. Produce the plots described below and explain how they allow you to identify the swapped samples.
  
i. (High volume) scatter plot matrix. 




ii. A heatmap of the first 100 genes (you can try more but it gets slow).




iii. Compute the Pearson correlation of the samples and plot the results using a heatmap.




iv. Scatterplot the six data samples with respect to the first two principal components and label the samples.




> Hint: If using the base graphics function `plot()`, a subsequent call to the `text()` function can be used to label each point with the corresponding sample name.

### c) (2pt) Microarray Differential Expression

Fix the label swap identified in question 1b. We want to swap b1 <--> c2. Revisit one or more elements of question 1b to sanity check before proceeding. 




Now use this data to do a differential expression analysis with `limma`.




Package these results in a data frame with six columns:

* probe.id - The array probe id.

* gene.id - The id of the gene which the probe overlaps (see below).

* p.value - The raw p-value for the probe.

* q.value - The BH corrected p-value, aka the q-value.

* log.fc - The log fold change which is the column called "logFC" in the limma results table.

* test.stat - The test statistics which for limma is the moderated t statstic. This is the column called "t" in the limma results table.

>  The gene id can be retrieved using the `yeast2.db` package from Bioconductor. In particular, the `yeast2ORF` object available after loading `yeast2.db` contains the mapping between probe IDs and yeast gene ids. Assuming you have a character vector of probes called `probe.ids`, the gene IDs can be retrieved using `gene.ids <- unlist(mget(probe.ids, yeast2ORF))`.

Remove any rows with probes which don't map to genes. You'll be able to find these because they will have `NA` as their gene id. Work with this data.frame to answer the questions below.




i. How many probes did we start with and how many remain after removing probes without gene ids?




ii. Illustrate the differential expression between the batch and the chemostat samples for the top hit (i.e., probe with the lowest p- or q-value).




iii. How many probes are identified as differentially expressed at a false discovery rate (FDR) of 1e-5 (note: this is a FDR cutoff used in the original paper)?




iv. Save your results for later with `write.table()`.

> When using write.table to save the data, you need to pass the arguments `row.names = TRUE, col.names = NA` to make sure the headers are written correctly.


__MORE QUESTIONS COMING SHORTLY. WE ARE RE-RUNNING THROUGH EVERYTHING OURSELVES.__
