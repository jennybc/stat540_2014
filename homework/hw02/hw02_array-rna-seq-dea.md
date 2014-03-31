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

See the [homework submission instructions](../hw_submission-instructions.html).

## Coaching

Some helpful tips can be found in [this coaching document](hw02_coaching.html).




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

* test.stat - The test statistics which for limma is the moderated t statistic. This is the column called "t" in the limma results table.

>  The gene id can be retrieved using the `yeast2.db` package from Bioconductor. In particular, the `yeast2ORF` object available after loading `yeast2.db` contains the mapping between probe IDs and yeast gene ids. Assuming you have a character vector of probes called `probe.ids`, the gene IDs can be retrieved using `gene.ids <- unlist(mget(probe.ids, yeast2ORF))`.

Remove any rows with probes which don't map to genes. You'll be able to find these because they will have `NA` as their gene id. Work with this data.frame to answer the questions below.




i. How many probes did we start with and how many remain after removing probes without gene ids?




ii. Illustrate the differential expression between the batch and the chemostat samples for the top hit (i.e., probe with the lowest p- or q-value).




iii. How many probes are identified as differentially expressed at a false discovery rate (FDR) of 1e-5 (note: this is a FDR cutoff used in the original paper)?




iv. Save your results for later with `write.table()`.

> When using write.table to save the data, you need to pass the arguments `row.names = TRUE, col.names = NA` to make sure the headers are written correctly.




## Q2) RNA-Seq Analysis

We have aligned the RNA-Seq library using the [Stampy](http://www.well.ox.ac.uk/project-stampy) aligner and generated count data. The data file is available as [stampy.counts.tsv](../../examples/yeastPlatforms/data/stampy.counts.tsv). In this question you will use this data to do a differential expression analysis using different packages from Bioconductor.

### a) (1pt) Load RNA Count Data and Sanity Check

Load the count data using `read.table`; you will need to pass the arguments `header=TRUE` and `row.names=1`. 

i) What are dimensions of the dataset? In addition to reporting number of rows and columns, make it clear what rows and columns represent. What is the difference between the rows of this dataset versus rows of the array data in question 1a?




ii) Do a sanity check to make sure there is no sample swap by plotting a heatmap of the sample correlations.




### b) (2pt) `edgeR` Differential Expression Analysis

Now you will use `edgeR` to identify differentially expressed genes between the batch medium vs. chemostat conditions.

i)  Recall that `edgeR` needs to estimate the dispersion parameter in the negative binomial model using an empirical Bayes method. Estimate the dispersion parameters using `estimateGLMCommonDisp`, `estimateGLMTrendedDisp` and `estimateGLMTagwiseDisp`. Plot the tagwise dispersion against log2-CPM (counts per million).  

> Seminar 7 was now corrected so that the design matrix is used as an argument of `estimateGLMTrendedDisp`.




ii)  Use the glm functionality of `edgeR`, i.e. use the `glmFit` function, to identify differentially expressed genes between conditions. 




Package these results in a data.frame called 'edger.results' with five columns:

* gene.id - The id of the gene which reads were aligned to.

* p.value - The raw p-value for the gene.

* q.value - The BH corrected p-value, aka the q-value.

* log.fc - The log fold change which is the column called "logFC" in the `edgeR` results table.

* test.stat - The test statistic, which for `edgeR` is a likelihood ratio. This is the column called "LR" in the `edgeR` results table.

Save your results for later with `write.table()` in file called `stampy.edger.results.tsv`.




iii) How many genes are differentially expressed between conditions at a false discovery rate (FDR) of 1e-5?




iv) How many genes are differentially over-expressed in chemostat compared to batch medium samples at a false discovery rate (FDR) of 1e-5?




### c) (2pt) `DESeq` Differential Expression Analysis

Now you will use `DESeq` to identify differentially expressed genes between the batch medium vs. chemostat conditions.

i)  `DESeq` also needs to estimate the dispersion. Use `estimateSizeFactors` and `estimateDispersions` to normalize the data. Plot the estimated dispersions against the mean normalized counts.




ii)  Use the negative binomial test of `DESeq`, i.e. use the `nbinomTest` function, to identify differentially expressed genes between conditions. Note that the output of this function does not return results ordered by p-values or logged fold-changes. You can manually reorder the results if you want (not required for this homework).

Package these results in a data.frame called 'deseq.results' with four columns:

* gene.id - The id of the gene which reads were aligned to.

* p.value - The raw p-value for the gene.

* q.value - The BH corrected p-value, aka the q-value.

* log.fc - The log fold change which is the column called "log2FoldChange" in the `DESeq` results table.

Save your results for later with `write.table()` in file called `stampy.deseq.results.tsv`.




iii) How many genes are differentially expressed between conditions at a false discovery rate (FDR) of 1e-5?




iv) How many differentially expressed genes are identified by both 'edgeR' and 'DESeq'?

> The function `intersect` -- which finds the elements common to two sets (vectors) -- will be helpful.




### d) (2pt) `voom` Differential Expression Analysis

Now you will use `voom+limma` to identify differentially expressed genes between the batch medium vs. chemostat conditions.

i)  `voom` normalizes the counts before it converts counts to log2-cpm. Use `calcNormFactors` to normalize counts.




ii)  Use `voom' to convert count data into logged CPM data and then use 'limma' to identify differentially expressed genes between conditions. 

Package these results in a data.frame called 'voom.limma.results' with five columns:

* gene.id - The id of the gene which reads were aligned to.

* p.value - The raw p-value for the gene.

* q.value - The BH corrected p-value, aka the q-value.

* log.fc - The log fold change which is the column called "logFC" in the `edgeR` results table.

* test.stat - The test statistic, which is the column called "t".

Save your results for later with `write.table()` in file called `stampy.limma.results.tsv`.




iii) How many genes are differentially expressed between conditions at a false discovery rate (FDR) of 1e-5?




iv)  What fraction of the genes identified using `voom+limma` are also found by `edger` and `DESeq` methods? For example if the DE analysis using `voom+limma` found 1000 genes and both `edgeR` and `DESeq`  found 500 of these, the fraction of genes found would be $\frac{500}{1000}=0.5$.




### e) (3pt) Comparison of Differential Expression Analyses
 
Now that we have the results of the differential expression analysis performed by three popular methods, we are going to compare and illustrate the results.

i) In previous questions, we noticed that different methods identified different differentially expressed genes. Create a Venn diagram showing all genes identified as differentially expressed by `edgeR`, `DESeq`, and `voom+limma`. Check your if your answers to questions 2c-iv, and 2d-iv are correct.

> The Venn diagram can be drawn using the `VennDiagram` package. It's a little obtuse if you want to plot to screen (or embed image using knitr), but the following code should get you started. Also be aware there is an argument called `force.unique`, which defaults to TRUE, that determines how elements that appear more than once in a set are handled when forming the Venn counts (in particular useful in question 3). 


```r
library(VennDiagram)
```

```
## Loading required package: grid
```

```r

# Fake some gene names for 4 different methods.  Note that in this example,
# I'm comparing 4 different sets so that you can see how to handle more
# complex cases.

method1.de.genes <- c("A", "B", "C")

method2.de.genes <- c("A", "B", "D", "E", "F")

method3.de.genes <- c("A", "B", "D", "E")

method4.de.genes <- c("A", "V", "E", "F")

# Put the things you want to plot in a list. The names in the list will be
# put on the plot.
de.genes <- list(Method1 = method1.de.genes, Method2 = method2.de.genes, Method3 = method3.de.genes, 
    Method4 = method4.de.genes)

# Start a new plot
plot.new()

# Draw the Venn diagram. Note the argument `filename=NULL` tells it to
# create a plot object instead of outputting to file.
venn.plot <- venn.diagram(de.genes, filename = NULL, fill = c("red", "blue", 
    "green", "yellow"))

# Draw the plot on the screen.
grid.draw(venn.plot)
```

![plot of chunk vennDiagram](figure/vennDiagram.png) 





ii) Using the function `plotSmear` function from `edgeR`, you can look at a scatterplot of observed differential expression (y-axis) against overall abundance (x-axis), both axes logarithmically transformed -- to check that putative DE genes seem plausible. Create a smear plot. Within this plot, identify the set of genes which are differentially expressed at an FDR of 1e-5 using all three methods (i.e., the q-values estimated by `edgeR`, `DESeq`, and `voom` are below 1e-5). Explain how you interpret this plot. Do you see reasonable results?

> Use the output of `DGEList` as the object of `plotSmear`. Use de.tags to highlight genes selected by all methods.




iii) There are two genes identified by `edgeR` and `voom+limma` but not by `DESeq`. Illustrate the logged counts of them. Compare the (log) counts of these two genes with those of two genes identified by the three methods (see example below)

> As an example, I illustrate two gene that were identified by all three methods. The function `setdiff` is helpful to find differences between two sets.









```r
library(lattice)
featureMe <- c("YDR384C", "YDR345C")
(featureCounts <- counts[featureMe, ])
```

```
##           b1   b2   b3   c1   c2   c3
## YDR384C  176  243  182 3332 3778 4531
## YDR345C 6155 8629 6357  322  345  462
```

```r
featureDat <- data.frame(gene.id = factor(rep(rownames(featureCounts), ncol(featureCounts))), 
    cond = factor(rep(groups, each = nrow(featureCounts))), log.count = log2(unlist(featureCounts)))
stripplot(gene.id ~ log.count, featureDat, groups = cond, auto.key = TRUE, jitter = TRUE)
```

![plot of chunk two-de-gene-demo](figure/two-de-gene-demo.png) 

```r

# Using the example created before to illustrate the `setdiff` function
setdiff(method1.de.genes, method2.de.genes)  #'C' is present in Method1 but not in Method2
```

```
## [1] "C"
```





## Q3) Compare DEA results between RNA-Seq and array

In question 1, you performed a DEA of array data using `limma`. In question 2, you used different methods to perform DEA of RNA-Seq data. In particular, in this question, you will focus on the DEA using `edgeR` (question 2b) to compare the results of RNA-Seq DEA with those of the array DEA . 

> Remember that you've packaged your results and saved them using `write.table()`. If the data.frames containing the results of these analyses are no longer in your workspace, load them using `read.table()`.

i) Use a Venn diagram to display the overlap and non-overlap of the __genes__ identified as differentially expressed at an FDR of 1e-5 by these analyses (i.e., array vs `edgeR` differentially expressed genes).

> Note that the number of __probes__ that you identified as differentially expressed at an FDR of 1e-5 in question 1c(iii) may be different from the number of __genes__ identified as differentially expressed at the same FDR! Why? The use of the argument `force.unique` of `venn.diagram()` is crucial here.




ii) As expected, more genes were identified as differentially expressed using RNA-Seq data. In this question, you will examine the difference between the q-values from both analyses (i.e., array and `edgeR`) by overlaying density plots of the q-values from each analysis.

  * To respond to this question, make two plots. One plot that includes the densities of q-values of the genes analyzed by both platforms (i.e., genes shared by both data frames), and another plot that includes the densities of q-values of ALL genes analyzed by at least one of the platforms.
  
Make some observations about the strengths of these two platforms.
  






iii) We provide a data set with array expression and count data for 5 interesting genes; below is also code to load it and a figure depicting it. Consult the array and the `edgeR` DEA results from your previous analyses for these genes. For each gene, state its status with respect to these analyses, i.e. where it falls in those Venn diagrams. Comment on the results and plots.

The data set is on GitHub and in the course webpage:

  * <https://raw.github.com/jennybc/stat540_2014/master/examples/yeastPlatforms/data/featGenesData-q3-DPUT.txt>
  * <http://www.ugrad.stat.ubc.ca/~stat540/examples/yeastPlatforms/data/featGenesData-q3-DPUT.txt>
  
Load it with a line like this:  
```
jDat <- dget("featGenesData-q3-DPUT.txt")
```




![plot of chunk stripplot-five-interesting-genes](figure/stripplot-five-interesting-genes.png) 


