The aim of the study was to "generate gene expression profiles of
purified photoreceptors at distinct developmental stages and from
different genetic backgrounds". The experimental units were mice and
the microarray platform was Affymetrix mouse genomic expression array
430 2.0.

For more information on this study, please refer to the 2006 publication:
http://www.ncbi.nlm.nih.gov/pubmed/16505381

The data is also directly accessible from GEO:
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4051

TO DO: can we track down provenance of the files we start with? And/or
can we start over with data from a link above and preserve the entire
process?

There are two main sources of information. Look for data files in
"data/photoRec/".

1. The gene expression data itself. The GSE4051_data_RAW.txt file
contains expression values of 29949 probes from photoreceptor cells in
39 mice samples. A "cleaned" version, where the columns = variables
have been rearranged rationally, is given in GSE4051_data.txt. See the
script 'caseStudies/photoRecPreprocess/02-cleanData.R' for details.

2. The metadata file GSE4051_design_RAW.txt describes the experimental
condition for each sample. Gene expression was studied at 5 different
developmental stages: day 16 of embryonic development (E16), postnatal
days 2,6 and 10 (P2, P6 and p10) as well as 4_weeks. Each of these 5
experimental conditions are tested in wild type Nrl mice and knockout
Nrl mice. A "cleaned" version, with variables renamed and factors and
row order rationalized, is given in GSE4051_design.txt. The same info
is preserved in two R-specific formats in files
GSE4051_design_DPUT.txt and GSE4051_design.robj. These files have the
advantage that, upon import, e.g. the levels of the developmental
stage factor will be in chronological order (not alphabetical). See
the script 'caseStudies/photoRecPreprocess/01-cleanDesign.R' for
details.

There are several derived datasets, created from processing the above.

GSE4051_MINI.txt holds data for 3 randomly selected probesets, renamed
for fun, transposed into convenient column variables and stored
together in a data.frame with the experimental condition information
for each sample. See the script `03-createMiniDataset.R` for how it
was created.

TO DO: simple differential expression analysis

TO DO: use diff exp analysis to pick 3 more interesting probesets for
a new mini dataset

HOW TO READ DATA/DESIGN
------------------------

WARNING: It is your responsibility to make sure the working directory
is set to where these files live or to edit paths accordingly below!

Raw data and design:

prDat <- read.table("GSE4051_data_RAW.txt", sep = "\t", header = T, row.names = 1)
prDes <- read.table("GSE4051_design_RAW.txt", sep = "\t", header = T, row.names = 1)

After the completion of data import, both "data" and "design" objects should be of class "data.frame".

> class(prDat)
[1] "data.frame"
> str(prDat)
'data.frame':	29949 obs. of  39 variables:
 $ Sample_35: num  7.15 9.22 10.06 8.35 8.45 ...
 $ Sample_32: num  7.54 9.53 9.92 8.78 8.57 ...
< ... snip, snip ... >
 $ Sample_3 : num  7.16 9.55 9.84 8.33 8.5 ...
 $ Sample_14: num  7.09 9.56 9.88 8.57 8.59 ...

> class(prDes)
[1] "data.frame"
> str(prDes)
'data.frame':	39 obs. of  2 variables:
 $ developmentStage   : Factor w/ 5 levels "4_weeks","E16",..: 1 1 1 1 1 1 1 1 2 2 ...
 $ genotypeOrVariation: Factor w/ 2 levels "Nrl_deficient",..: 1 1 1 1 2 2 2 2 1 1 ...


Cleaned data and design (saved in various formats):

prDat <- read.table("photoRec/GSE4051_data.txt")
str(prDat)
## 'data.frame':	29949 obs. of  39 variables:
##  $ Sample_20: num  7.24 9.48 10.01 8.36 8.59 ...
##  $ Sample_21: num  7.41 10.02 10.04 8.37 8.62 ...
## ...
##  $ Sample_2 : num  7.35 9.66 9.91 8.4 8.37 ...
##  $ Sample_9 : num  7.32 9.8 9.85 8.4 8.46 ...

prDes <- read.table("GSE4051_design.txt")
str(prDes)
## 'data.frame':	39 obs. of  3 variables:
##  $ sample  : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType   : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...

In the above case, note that the factor levels for devStage and gType
may not be as you want. Wild type ('wt') is not the reference level
and the developmental stages are not in chronological order. Set
explicitly or import from a format that preserves factor levels (see
below).

prDes <- dget("GSE4051_design_DPUT.txt")
load("GSE4051_design.robj")

Both the dput/dget and the save/load approaches will leave the prDes
data.frame like so:
str(prDes)
## 'data.frame':	39 obs. of  3 variables:
##  $ sample  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...


Mini dataset:

The usual problem with factor level order occurs here:
> read.table("GSE4051_MINI.txt")
> str(kDat)
'data.frame':	39 obs. of  6 variables:
 $ sample    : int  20 21 22 23 16 17 6 24 25 26 ...
 $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
 $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
 $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
 $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
 $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...

Importing from these R-specific formats preserves factor levels:
dget("GSE4051_MINI_DPUT.txt")
load("GSE4051_MINI.robj")
str(kDat)
'data.frame':	39 obs. of  6 variables:
 $ sample    : num  20 21 22 23 16 17 6 24 25 26 ...
 $ devStage  : Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
 $ gType     : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
 $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
 $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
 $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...

