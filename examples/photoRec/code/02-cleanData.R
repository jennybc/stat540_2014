whereAmI <- "/Users/jenny/teaching/2012-2013/STAT540"

prDat <- read.table(file.path(whereAmI,
                              "rmd/data/photoRec/GSE4051_data_RAW.txt"),
                    sep = "\t", header = T, row.names = 1)

str(prDat)
## 'data.frame':	29949 obs. of  39 variables:
##  $ Sample_35: num  7.15 9.22 10.06 8.35 8.45 ...
##  $ Sample_32: num  7.54 9.53 9.92 8.78 8.57 ...
## ...
##  $ Sample_3 : num  7.16 9.55 9.84 8.33 8.5 ...
##  $ Sample_14: num  7.09 9.56 9.88 8.57 8.59 ...

## bring the design data.frame in

## showing various options, given all the different ways we preserved
## it after cleaning

## write.table / read.table
prDes <- read.table(file.path(whereAmI,
                              "rmd/data/photoRec/GSE4051_design.txt"))

str(prDes)
## 'data.frame':	39 obs. of  3 variables:
##  $ sample  : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType   : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...

## PROBLEM: OUR FACTOR LEVELS ARE NOT IN THE ORDER WE WANTED!
prDes$devStage <- factor(prDes$devStage,
                         levels = c("E16", "P2", "P6", "P10", "4_weeks"))
prDes$gType <- factor(prDes$gType, levels = c("wt", "NrlKO"))
str(prDes)
## 'data.frame':	39 obs. of  3 variables:
##  $ sample  : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
## PROBLEM SOLVED.

## dput / dget
prDes <-
    dget(file = file.path(whereAmI,
         "rmd/data/photoRec/GSE4051_design_DPUT.txt"))
str(prDes)
## 'data.frame':	39 obs. of  3 variables:
##  $ sample  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
## NOTICE the factor levels set during cleaning are preserved.

## save / load
load(file = file.path(whereAmI,
     "rmd/data/photoRec/GSE4051_design.robj"))
str(prDes)
## 'data.frame':	39 obs. of  3 variables:
##  $ sample  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...

## I want to rearrange the columns / variables of the gene expression
## data to match the order of the rows in the design data.frame

## see they are not the same!
cbind(design = rownames(prDes), data = names(prDat))

## find out where the rownames of jDat appear in those of prDes
(alpha <- match(names(prDat), rownames(prDes)))

## if I order variables of prDat with alpha, will that work?
## see if these two columns are "harmonized"

## visual check   ..... NO!
cbind(names(prDat)[alpha], rownames(prDes))

## automated check
identical(names(prDat)[alpha], rownames(prDes))
## [1] FALSE

## let's try it the other way around
beta <- match(rownames(prDes), names(prDat))

## if I order jDat with beta, will that work?

## visual check   ..... YES!
cbind(names(prDat)[beta], rownames(prDes))

## automated check
identical(names(prDat)[beta], rownames(prDes))
## [1] TRUE

prDat <- prDat[, beta]

## I consider this data.frame clean and ready for downstream analyses
## .... time to save our work

## there is no compelling reason to use any strategy other than
## write.table here (e.g. no factor levels to preserve)
write.table(prDat,
            file = file.path(whereAmI,
            "rmd/data/photoRec/GSE4051_data.txt"),
            quote = FALSE)

prDat <- read.table(file.path(whereAmI,
                              "rmd/data/photoRec/GSE4051_data.txt"))


peek(prDat)

kDat <- data.frame(prDes, jDat)
str(kDat)
## 'data.frame':	39 obs. of  5 variables:
##  $ developmentStage   : Factor w/ 5 levels "4_weeks","E16",..: 1 1 1 1 1 1 1 1..
##  $ genotypeOrVariation: Factor w/ 2 levels "Nrl_deficient",..: 1 1 1 1 2 2 2 2..
##  $ crabHammer         : num  9.68 9.13 9.74 9.82 9.96 ...
##  $ eggBomb            : num  7.2 7.17 7.11 6.56 7.87 ...
##  $ poisonFang         : num  6.98 7.35 7.08 7.04 6.99 ...

peek(kDat)
##           developmentStage genotypeOrVariation crabHammer eggBomb poisonFang
## Sample_12          4_weeks       Nrl_deficient      9.129   7.165      7.350
## Sample_36          4_weeks           wild_type      9.960   7.866      6.993
## Sample_39          4_weeks           wild_type     10.200   7.003      7.320
## Sample_14               P2       Nrl_deficient      9.572   6.138      7.250
## Sample_5                P2       Nrl_deficient      8.925   6.269      7.405
## Sample_24               P2           wild_type      8.869   6.587      7.508
## Sample_7                P6       Nrl_deficient      8.803   6.188      7.754



set.seed(5782)

y <- prDat[sample(1:nrow(prDat), 1), ]

set.seed(5782)
(getMe <- sample(1:nrow(prDat), 3))
## [1]  2668 20534 17554

## I really wanted to demonstrate the use of subset() but the subset
## argument expects a logical, which is awkward for me here ...
jDat <- prDat[getMe, ]
str(jDat)
## 'data.frame':	3 obs. of  39 variables:
##  $ Sample_35: num  8.45 6.16 7.2
##  $ Sample_32: num  10.25 8.17 7

## renaming the probesets from the Pokemon move list to make them more
## catchy!
## http://pokemondb.net/move/all

jDat <- data.frame(t(jDat))
str(jDat)
## 'data.frame':	39 obs. of  3 variables:
##  $ X1419037_at: num  8.45 10.25 8.52 9 8.21 ...
##  $ X1446244_at: num  6.16 8.17 6.76 7.08 6.53 ...
##  $ X1440765_at: num  7.2 7 8.58 8.09 7.43 ...
names(jDat) <- c("crabHammer", "eggBomb", "poisonFang")
str(jDat)
## 'data.frame':	39 obs. of  3 variables:
##  $ crabHammer: num  8.45 10.25 8.52 9 8.21 ...
##  $ eggBomb   : num  6.16 8.17 6.76 7.08 6.53 ...
##  $ poisonFang: num  7.2 7 8.58 8.09 7.43 ...


## let's make the names a bit shorter
(foo <- names(kDat))
## [1] "developmentStage"    "genotypeOrVariation" "crabHammer"
## [4] "eggBomb"             "poisonFang"          "sample"

names(kDat) <- c("devStage", "gType", foo[-(1:2)])

peek(kDat)
##           devStage         gType crabHammer eggBomb poisonFang sample
## Sample_9   4_weeks Nrl_deficient      9.822   6.558      7.043      9
## Sample_38  4_weeks     wild_type      9.767   6.608      7.329     38
## Sample_39  4_weeks     wild_type     10.200   7.003      7.320     39
## Sample_20      E16     wild_type     10.220   7.462      7.370     20
## Sample_22      E16     wild_type      9.642   6.720      7.350     22
## Sample_13      P10 Nrl_deficient      9.838   7.228      7.459     13
## Sample_18      P10 Nrl_deficient     10.140   7.438      7.363     18

## finally, rearrange the variables rationally
(yo <- match("sample", names(kDat)))    # sample is 6th variable
names(kDat)[c(yo, seq_along(kDat)[-yo])]

kDat <- kDat[c(yo, seq_along(kDat)[-yo])]

str(kDat)
## 'data.frame':	39 obs. of  6 variables:
##  $ sample    : num  11 12 2 9 36 37 38 39 16 17 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 1 1 1 1 1 1 1 1 2 2 ...
##  $ gType     : Factor w/ 2 levels "Nrl_deficient",..: 1 1 1 1 2 2 2 2 1 1 ...
##  $ crabHammer: num  9.68 9.13 9.74 9.82 9.96 ...
##  $ eggBomb   : num  7.2 7.17 7.11 6.56 7.87 ...
##  $ poisonFang: num  6.98 7.35 7.08 7.04 6.99 ...

peek(kDat)
##           sample devStage         gType crabHammer eggBomb poisonFang
## Sample_2       2  4_weeks Nrl_deficient      9.744   7.107      7.075
## Sample_9       9  4_weeks Nrl_deficient      9.822   6.558      7.043
## Sample_17     17      E16 Nrl_deficient     10.140   7.065      7.005
## Sample_13     13      P10 Nrl_deficient      9.838   7.228      7.459
## Sample_35     35      P10     wild_type      8.449   6.155      7.201
## Sample_3       3       P2 Nrl_deficient      9.414   6.166      7.200
## Sample_7       7       P6 Nrl_deficient      8.803   6.188      7.754


peek(kDat)
##    sample devStage         gType crabHammer eggBomb poisonFang
## 4       9  4_weeks Nrl_deficient      9.822   6.558      7.043
## 11      6      E16 Nrl_deficient     10.340   7.017      6.735
## 13     21      E16     wild_type     10.020   6.890      7.177
## 18     18      P10 Nrl_deficient     10.140   7.438      7.363
## 31     27       P2     wild_type      8.613   6.800      7.843
## 32      1       P6 Nrl_deficient      8.920   6.286      7.378
## 35      7       P6 Nrl_deficient      8.803   6.188      7.754

stripplot(~ crabHammer, kDat)

stripplot(~ crabHammer, kDat,
          groups = devStage, auto.key = TRUE)

stripplot(devStage ~ crabHammer, kDat,
          groups = gType, auto.key = TRUE)
