Homework 01
======================================================================

In this assignment you will be analyzing a publicly-available expression study of mouse brain tissue, run on the Affymetrix MGU74Av2 platform. You can read about the study at <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7191>.

The motivation is to study a gene [S1P2, sphingosine-1-phosphate receptor 2](http://www.ncbi.nlm.nih.gov/pubmed/11553273), which, when mutated, results in seizures. In contrast, mutation of the related gene S1P3 does not do this. The study includes two brain regions, hippocampus and neocortex, and three mouse strains (wild type, S1P2 mutant, and S1P3 mutant). In addition the authors provide information on the gender of the mice. Paul was also able to extract "processing date" information from the raw data files (the authors did not explicitly provide this).
<!---Obviously the researchers were looking for genes that change expression in S1P2 mice but not in S1P3 mice, relative to the wild types.--->

You can find two data files for this study [here](../examples/mouseBrain/), specifically the [expression data](../examples/mouseBrain/data/GSE7191-data.txt) and the [design](../examples/mouseBrain/data/GSE7191-design.txt) of the experiment. The data have been re-preprocessed with RMA from the Bioconductor `affy` package, so it is on a log2 scale. Some probes that are considered "controls" were also removed, so overall it's a little different than the processed data provided via GEO. Our focus here is on quality control, exploration, and differential expression analysis.

## Evaluation

The assignment has 25 points and represents 25% of your overall course mark. Overall presentation and mechanics is worth 5 points, with full points awarded if it's exceptional. The remaining 20 points are spread among the specific questions, as detailed below. We will try to give partial credit when/if we can figure out what you're trying to do.

Each student must submit their own work: no group efforts. It's fine to talk to fellow students and have discussion on the Google Group. If someone is really helpful, give them some credit. You have definitely crossed the line if there's any copying and pasting of code or wording. 

## How to submit the work

*to be written ... write it in R markdown, compile to markdown and HTML, put all that in private repo on Github, share links with TA*




## Your mission

### Q0 **(0 pts)** Intake

Load the data. Smell test it. Fiddle with factor levels and do any other light cleaning you deem necessary. You can set `include = FALSE` for some/all of this code, but DO IT.
















### Q1 **(2 points)** What are the basic characteristics of the data and meta-data? 

#### Q1a: How many probes? How many samples?




#### Q1b: What is the breakdown of samples for Genotype, BrainRegion, Sex, and DateRun?

For starters, cross-tabulate Genotype and BrainRegion. Then do some cross-tabulation involving DateRun and, optionally, Sex to investigate if these "nuisance" factors are confounded with the experimental factors. Make sure you've included counts for all 4 individual factors somewhere, even if it's just in the margins of a cross-tabulation. How do you feel about the experimental design? Hint: `table()` and `addmargins()` will help.




#### Q1c: Create a plot showing the gene expression data for one probe.

Use position, panels or facets, and color to convey the Genotype, BrainRegion, and Sex of the samples.




#### Q1d: Report the average expression of the selected probe for all possible combinations of Genotype, BrainRegion and Sex.




### Q2 **(4 points)** Examine the sample correlation matrix.

#### Q2a: Depict the sample-to-sample correlations in a heatmap.

At the very least, order the samples by run date; within run date, you could also sort on other factors. What do you see about batch effects, outliers, etc.? Hint: `cor()` and `heatmap()` are helpful.




#### Q2b: Identify the outlier sample.

Try to go beyond merely "eyeballing" it, i.e. make a quantitative statement about how it "sticks out" from the other samples.




#### Q2c: Examine the outlier in the context of its experimental group.

Which group is the outlier in, in terms of Genotype, BrainRegion, and Sex? Scatter plot these samples against each other and comment on the plot. Try to address overplotting! Hint: `splom()` from `lattice` is helpful.




### Q3 **(4 points)** Normalization, what it can and can't do.

#### Q3a: Make boxplots.

Make a pair of plots: each plot will show gene expression boxplots, one per sample. Identify the outlier sample with a different color, if you can. Comment on these plots, in general and with respect to the outlier.

Here are the 2 plots:

  * The data as provided
  * The data after you apply quantile normalization.

Hint: `normalize.quantiles()` in the Bioconductor package `preprocessCore` is helpful. So is `boxplot()`.




#### Q3b: Did normalization fix the outlier?

With the quantile-normalized data, re-make the heatmap of the sample correlation matrix and re-compare the outlier to its experimental group. Is everything OK now?




#### Q3c: Form a dataset that omits the outlier and quantile normalize it.




#### Q3d Re-make the heatmap of the sample correlation matrix, now that the worst outlier is gone. Interpret what you see.




#### Q3e: Remake the expression boxplots for all samples before moving on to differential expression analysis.

How many samples remain? 




### Q4 **(5 points)** For each probe, test if it has differential expression across the three genotypes *within the neocortex brain region only*.

You should be using data with the worst outlier removed, after quantile normalization and, for this question, restricted to neocortex.

#### Q4a: Write out, in an equation or English or, ideally, both, the model you are fitting.

In the context of that model, what statistical test(s) are you conducting? (You can gloss over the neocortex restriction here.)







#### Q4b: Explore your hits.

Display the expression data for the top 50 probes in a heat map. Order the probes by p-value; order the samples by genotype.




#### Q4c: Count your hits.

How many probes have unadjusted p-value < 10e-4? If we took the top 50 probes as our "hits", what is the (estimated) false discovery rate? How many of these hits do we expect to be false discoveries?




#### Q4d: Plot the gene expression data for a few top hits and a few boring probes.

What are the p-values associated with these probes? Comment on the plots.







#### Q4e: Find probes where the expression in the S1P3 knockout is the same as that of the wild type.

How many hits do you find if you control FDR at 0.10?




### Q5 **(5 points)** Differential expression analysis for Genotype * BrainRegion.

You should be using data with the worst outlier removed, after quantile normalization and for both brain regions.

#### Q5a: Fit a 3x2 full factorial model.

Test for any effect of Genotype and/or BrainRegion, i.e. test your model against a model with just an intercept. How many probes have a BH-adjusted p-value, a.k.a. q-value, less than 10e-4?




#### Q5b: Test the null hypothesis that BrainRegion doesn't matter, i.e. that all terms involving BrainRegion are zero.

How many probes have a BH-adjusted p-value less than 0.1?




#### Q5c: Highlight some probes where BrainRegion does and does not matter.

Using the results from Q5b, plot and describe some results for a couple of hits and non-hits.




#### Q5d: Test the null hypothesis that Genotype doesn't matter, i.e. that all terms involving Genotype are zero.

How many probes have a BH-adjusted p-value less than 0.1? Compare these results to those obtained for BrainRegion in Q5b, either based on this count or based on all the p-values. What do you conclude about the relative magnitude of the influence of brain region vs. genotype on gene expression?



