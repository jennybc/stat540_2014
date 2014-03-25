Homework 2 coaching
======================================================================

Contributor: Gabriela Cohen Freue

## Question 1: Microarray Analysis

Nothing new here, just a simple limma analysis of microarray data! This is very similar to what you have done in HW1.

### Q1a) Load the data
Throughout the homework, you will need to load and use different data sets and data frames with results. 
- Choose good names for your objects. For example, I call these data `array.dat`.
- Pass the arguments `header=TRUE` and `row.names=1` to `read.table`.

### Q1b) Identify Sample Swap
Data exploration and quality control! This should always be the first step of your analysis (even when nobody tells you to do it). You have done this many times in previous seminars and homework. 
- Although you can `plot` to draw a scatter plot matrix, this is slow and hides some detail. Check `lattice` graphics in Seminar 3 for better graphic alternatives.

### Q1c) Microarray Differential Expression
You want to swap b1 <--> c2. Yes, this is what you have to prove in question 1b. But if something went wrong there, you can start fresh here. Question 3 depends on getting this part correct. So: 
- make sure you fix the column names before doing a differential expression analysis (DEA)
- you may want to reorder the columns in the data so that you have all batch samples together and all the chemostat together. This is not really needed. Just be careful about the order of the samples if you write the design matrix by hand
- check that things are ok now (redo some of the plots from Q1b) 
- we did not ask for this, but it is good practice to write the fixed version of the data to file for future use. If you don't do this then you would need to run the code to fix the swap everytime you work with the data

The DEA is a standard two group comparison with limma:
- create a factor for the two groups (either manually or using the column names). In any case, check that the factor you created reflects the order of the samples in the data
- define the design matrix. You don't need to download it, just use the factor you've created
- follow the hint given in the HW to identify the gene matching to each probe (i.e., unlist(mget(probe.ids, yeast2ORF))). It is important to use unlist here instead of as.character, in terms of preserving NAs as proper NAs not the string 'NA'.
- save these results to file using `write.table`. Specify `col.names=NA` in order to get the column headers to line up correctly in the file.

## Question 2: RNA-Seq Analysis

This is just [Seminar 7](http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar07_RNA-seq.html)! You will analize the same count data with three different methods. Use good names to store your results. Avoid calling your data frame `results`. Instead, you can use `edger.results` for example. Use consistent naming scheme.

### Q2b) Edge R
In an early version of Seminar 7, the design matrix was not passed as an argument of `estimateGLMTrendedDisp`. Note that the design is needed in all functions estimating the dispersion. This was corrected in Seminar 7 and the current version contains the correct analysis.

### Q2e) Comparison of DEA
ii) The hint of the Venn Diagrams says: 

> Since multiple probes in the array data can map to the same gene, you will obtain different answers depending on the value of this parameter.

Note that this part of the hint is not really useful in question 2, which is just about RNA-Seq data analysis. However, you need this part of the hint in question 3i.

iii) In this question, you need to use `setdiff`. Play with this function and note that, in general, for two sets A and B, setdiff(A,B) is not the same as setdiff(B, A).

## Question 3: Compare RNA-Seq with array DEA
Here you will compare one of the results from RNA-Seq DEA with the array DEA. You should have two data frames containing the results from the RNA-Seq and the array DEA, respectively. If these are not in the memory of your R session, you can load the results that you saved in previous questions. 

### Q3ii) Density plots
There are many ways to answer this question. One option is to use the function `merge` passing the argument `by="gene.id"` (gene.id are in both data frames). If you use this function, note that most of the columns in these data frames have the same name. Thus, it will also be helpful to add a *suffix* to the column names from each data frame. For example if you had the results stored in data frames called "array.results" and "edger.results", the merge command would be `merge(array.results, edger.results, by="gene.id", suffixes=c(".array", ".rnaseq"))`. Also, note that the RNASeq data contain more genes than those analyzed in the array data. The argument `all=TRUE` controls which genes are kept in the resulting data frame.

Why are we asking for two plots in this questions? From the previous analysis, you’ve noted that there are more DE genes in the RNA-Seq analysis. However, it is also true that more genes were measured and analyzed by RNA-Seq data. Thus, it is interesting to look at two sets of genes: those genes shared by both platforms (plot 1), and all studied genes (i.e., by at least one platform) (plot 2). Both plots illustrate 2 density curves, one for array's q-values and one for RNA-Seq's q-values.

### Q3iii) Interesting genes
In this question we provide the data and the plot of 5 interesting genes. You don't need to make any additional plot here. You need to interpret the plot given and you need to retreive information about these genes from your previous analysis.

