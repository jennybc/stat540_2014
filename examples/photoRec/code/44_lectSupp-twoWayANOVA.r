library(limma)
library(car)                            # Anova()

whereAmI <- "/Users/jenny/teaching/2012-2013/STAT540"

whereLectureLives <-
  "/Users/jenny/teaching/2012-2013/STAT540-jennyLocal/classMeetings/"

## code Rick developed for working with mlm objects, ie fitted multivariate
## regression
source(file.path(whereAmI, "rmd/caseStudies/photoRecDiffExp/anova.mlm.R"))

## design data.frame
load(file = file.path(whereAmI,
     "rmd/data/photoRec/GSE4051_design.robj"))
str(prDes)
## 'data.frame':	39 obs. of  3 variables:
##  $ sample  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...

## gene expression data.frame
prDat <- read.table(file.path(whereAmI,
                              "rmd/data/photoRec/GSE4051_data.txt"))
str(prDat, max.level = 0)
## 'data.frame':	29949 obs. of  39 variables:
##  $ Sample_20: num  7.24 9.48 10.01 8.36 8.59 ...
##  $ Sample_21: num  7.41 10.02 10.04 8.37 8.62 ...
## ...
##  $ Sample_2 : num  7.35 9.66 9.91 8.4 8.37 ...
##  $ Sample_9 : num  7.32 9.8 9.85 8.4 8.46 ...

##########################################################
## homegrown mlm code
##########################################################

## responses must be a matrix, one *row* per response
prMat <- t(as.matrix(prDat))
gType <- prDes$gType              # lesser of two evils
devStage <- prDes$devStage        # lesser of two evils
rFit <- lm(prMat ~ gType * devStage)

## summary(rFit) hangs unless you've loaded our homegrown code
## source()ed at top!
rfSumm <- summary(rFit)
rfSumm                                  # by default, info on first 2
                                        # fitted models
print(rfSumm, show = c(2, 4555, 29403)) # show gives more flexibility

str(rfSumm, max.level = 1)
str(rfSumm$Coef)

##########################################################
## testing for interaction
##########################################################

rfAnova <- anova(rFit)
plot(rfAnova)                           # things whizz by ....
rfAnova                                 # first 2 fits, by default

peek(rfAnova)

##########################################################
## make it easy to find various types of 'hits'
##########################################################

## get p-values for the the main effects and interaction
pVals <- rfAnova[ , c("gType:Pval", "devStage:Pval",
                      "gType:devStage:Pval")]
colnames(pVals) <- c("gType", "devStage", "gType:devStage")
splom(~ pVals, panel = panel.smoothScatter)
pValsTall <- data.frame(pVals = as.vector(pVals),
                        what = rep(colnames(pVals),
                        each = nrow(pVals)))
densityplot(~ pVals, pValsTall,
            group = what, auto.key = TRUE,
            plot.points = FALSE)

## note from the future: having this global p-value cutoff doesn't
## turn out to be as useful as you hope. I end up threshholding on the
## fly often below, so maybe eliminate this step here?
pValCutoff <- 0.15
pHits <- apply(pVals, 2, function(yo) yo <= pValCutoff)
table(gType = pHits[ , "gType"],
      devStage = pHits[ , "devStage"],
      interaction = pHits[ , "gType:devStage"])

##########################################################
## helper functions for finding and featuring genes
##########################################################

## handy for isolating and reshaping the data for a small set of
## probesets
extractIt <- function(myGenes) {
    miniDat <- t(prDat[myGenes, ])
    miniDat <- data.frame(gExp = as.vector(miniDat),
                          gene = rep(colnames(miniDat), each = nrow(miniDat)))
    miniDat <- data.frame(prDes, miniDat)
    miniDat
}

## handy for stripplotting the data for a small set of probesets with
## our 2x2 factorial mentality in mind
graphIt <- function(jDat, yFree = FALSE) {
    thePlot <- stripplot(gExp ~ devStage | gene, jDat,
                         type = c('p', 'a'), grid = TRUE,
                         group = gType, auto.key = TRUE,
                         jitter.data = TRUE)
    if(yFree) {
        thePlot <- update(thePlot,
                          scales = list(y = list(relation = "free")))
    }
    print(thePlot)
}

##########################################################
## finding some genes to feature in lecture
##########################################################

## boring
getMe <- which(rowSums(pHits) == 0)
length(getMe)
set.seed(111)
(boring <- sample(getMe, size = 3))
graphIt(extractIt(boring))
dev.print(pdf,
          paste0(whereLectureLives, "figs/gTypeByDevStage/stripplot-boring.pdf"),
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMat) == names(boring)[2]))
print(rfAnova, show = which(colnames(prMat) == names(boring)[2]))


## devStage YES NrlKO NO interaction NO
getMe <- which(pVals[ , "devStage"] < 0.0005 & # demand more here!
               !pHits[ , "gType"] &
               !pHits[ , "gType:devStage"])
length(getMe)                           # 1243
## set.seed(555)
## (devStageOnly <- sample(getMe, size = 12))
## through various random looks, I settled on these 3:
foo <- c("1456219_at", "1455007_s_at", "1445613_at")
(devStageOnly <- which(rownames(prDat) %in% foo))
names(devStageOnly) <- foo
graphIt(extractIt(devStageOnly))
dev.print(pdf,
          paste0(whereLectureLives, "figs/gTypeByDevStage/stripplot-devStage.pdf"),
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMat) == names(devStageOnly)[2]))
print(rfAnova, show = which(colnames(prMat) == names(devStageOnly)[2]))

## devStage NO NrlKO YES interaction NO
getMe <- which(pVals[ , "devStage"] > 0.8 &
               pVals[ , "gType"]  < 0.05 & # demand more here!
               !pHits[ , "gType:devStage"])
length(getMe)                           # 775
set.seed(111)
(gTypeOnly <- sample(getMe, size = 12))
## through various random looks, I settled on these 3:
foo <- c("1452142_at", "1429993_s_at", "1456468_x_at")
(gTypeOnly <- which(rownames(prDat) %in% foo))
names(gTypeOnly) <- foo
graphIt(extractIt(gTypeOnly))
dev.print(pdf,
          paste0(whereLectureLives, "figs/gTypeByDevStage/stripplot-gType.pdf"),
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMat) == names(gTypeOnly)[2]))
print(rfAnova, show = which(colnames(prMat) == names(gTypeOnly)[2]))
print(rfAnova, show = which(colnames(prMat) == names(gTypeOnly)[1]))
print(rfAnova, show = which(colnames(prMat) == names(gTypeOnly)[3]))

## devStage YES NrlKO YES interaction NO
getMe <- which(pVals[ , "devStage"] < 0.05 &
               pVals[ , "gType"]  < 0.05 &
               pVals[ , "gType:devStage"] > 0.7)
length(getMe)                           # 515
set.seed(111)
(mainEffOnly <- sample(getMe, size = 12))
## through various random looks, I settled on these:
foo <- c("1438815_at", "1438786_a_at", "1448690_at")
## more candidates:
## "1433716_x_at", "1434564_at", , "1427347_s_at", "1435176_a_at"
(mainEffOnly <- which(rownames(prDat) %in% foo))
names(mainEffOnly) <- foo
graphIt(extractIt(mainEffOnly))
dev.print(pdf,
          paste0(whereLectureLives, "figs/gTypeByDevStage/stripplot-mainYESinterNO.pdf"),
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMat) == names(mainEffOnly)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(mainEffOnly)))

## devStage YES NrlKO YES interaction YES
getMe <- which(pVals[ , "devStage"] < 0.05 &
               pVals[ , "gType"]  < 0.05 &
               pVals[ , "gType:devStage"] < 0.05)
length(getMe)                           # 1570
set.seed(111)
(exciting <- sample(getMe, size = 12))
## through various random looks, I settled on these:
foo <- c("1451763_at", "1419655_at", "1429225_at")
## more candidates:
## "1452601_a_at", "1433577_at"
(exciting <- which(rownames(prDat) %in% foo))
names(exciting) <- foo
graphIt(extractIt(exciting))
dev.print(pdf,
          paste0(whereLectureLives, "figs/gTypeByDevStage/stripplot-mainYESinterYES.pdf"),
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMat) == names(exciting)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(exciting)))

##########################################################
## detailed look at a boring gene
##########################################################

## http://goanna.cs.rmit.edu.au/~fscholer/anova.php

## https://stat.ethz.ch/pipermail/r-help/2011-April/273863.html

## http://psych.wisc.edu/moore/Rpdf/610-R12_UnbalFactorialBetw.pdf

## referred to multiple times ... server not responding last I checked
## http://prometheus.scp.rochester.edu/zlab/sites/default/files/InteractionsAndTypesOfSS.pdf

(luckyGene <- boring[2])
jDat <- data.frame(gType = prDes$gType,
                   devStage = prDes$devStage,
                   gExp = unlist(prDat[luckyGene, ]))
lmRes <- lm(gExp ~ gType * devStage, jDat)
summary(lmRes)

deviance(lmRes)

anova(lmRes)
##                Df  Sum Sq  Mean Sq F value Pr(>F)
## gType           1 0.02985 0.029848  1.5657 0.2208
## devStage        4 0.10365 0.025914  1.3594 0.2722
## gType:devStage  4 0.07191 0.017977  0.9430 0.4532
## Residuals      29 0.55283 0.019063

lmResAlt <- lm(gExp ~ devStage * gType, jDat)

anova(lmResAlt)
##                Df  Sum Sq  Mean Sq F value Pr(>F)
## devStage        4 0.10328 0.025819  1.3544 0.2739
## gType           1 0.03022 0.030225  1.5855 0.2180
## devStage:gType  4 0.07191 0.017977  0.9430 0.4532
## Residuals      29 0.55283 0.019063

## evidence of the slightly weird and undesirable nature of the Type I
## sums of squares approach ... see how the entries for gType and
## devStage depend on the order in which they appear in the model
## formula?

Anova(lmRes)
##                 Sum Sq Df F value Pr(>F)
## gType          0.03022  1  1.5855 0.2180
## devStage       0.10365  4  1.3594 0.2722
## gType:devStage 0.07191  4  0.9430 0.4532
## Residuals      0.55283 29

Anova(lmResAlt)
##                 Sum Sq Df F value Pr(>F)
## devStage       0.10365  4  1.3594 0.2722
## gType          0.03022  1  1.5855 0.2180
## devStage:gType 0.07191  4  0.9430 0.4532
## Residuals      0.55283 29

## see these Type II sums of squares results are invariant to the
## order in which the factors are given in the model

Anova(lm(gExp ~  gType * devStage, jDat,
         contrasts=list(gType = contr.sum, devStage = contr.sum)),
         type=3)
##                 Sum Sq Df    F value Pr(>F)
## (Intercept)    2753.78  1 1.4446e+05 <2e-16 ***
## gType             0.03  1 1.6862e+00 0.2043
## devStage          0.10  4 1.3679e+00 0.2693
## gType:devStage    0.07  4 9.4300e-01 0.4532
## Residuals         0.55 29

aovRes <- aov(gExp ~ gType * devStage, jDat)
summary(aovRes)
##                Df Sum Sq Mean Sq F value Pr(>F)
## gType           1 0.0298 0.02985   1.566  0.221
## devStage        4 0.1037 0.02591   1.359  0.272
## gType:devStage  4 0.0719 0.01798   0.943  0.453
## Residuals      29 0.5528 0.01906
## same table as anova(lmRes)

aovResAlt <- aov(gExp ~ devStage * gType, jDat)
summary(aovResAlt)
##                Df Sum Sq Mean Sq F value Pr(>F)
## devStage        4 0.1033 0.02582   1.354  0.274
## gType           1 0.0302 0.03022   1.585  0.218
## devStage:gType  4 0.0719 0.01798   0.943  0.453
## Residuals      29 0.5528 0.01906
## same table as anova(lmResAlt)

anova(lm(gExp ~ gType * devStage, jDat))
anova(lm(gExp ~ devStage * gType, jDat))

Anova(lm(gExp ~ gType * devStage, jDat))
Anova(lm(gExp ~ devStage * gType, jDat))


lmBig <- lm(gExp ~ gType * devStage, jDat)
deviance(lmBig)

lmSmall <- lm(gExp ~ gType + devStage, jDat)
deviance(lmSmall)

deviance(lmSmall) - deviance(lmBig)
lmBig$rank - lmSmall$rank

(numerator <- (deviance(lmSmall) - deviance(lmBig)) /
 (lmBig$rank - lmSmall$rank))

(denominator <- deviance(lmBig) / lmBig$df.residual)

(Fstat <- numerator / denominator)

pf(Fstat, lmBig$rank - lmSmall$rank, lmBig$df.residual, lower.tail =
   FALSE)

## Type I, also called "sequential" sum of squares:
## SS(A) for factor A.
## SS(B | A) for factor B.
## SS(AB | B, A) for interaction AB.

## This tests the main effect of factor A, followed by the main effect
## of factor B after the main effect of A, followed by the interaction
## effect AB after the main effects.

lmFull <- lm(gExp ~ gType * devStage, jDat)
lmNoInter <- lm(gExp ~ gType + devStage, jDat)
lmGtype <- lm(gExp ~ gType, jDat)
lmInter <- lm(gExp ~ 1, jDat)

myFtest <- function(lmBig, lmSmall) {
    numerator <- (deviance(lmSmall) - deviance(lmBig)) /
        (lmBig$rank - lmSmall$rank)
    denominator <- deviance(lmBig) / lmBig$df.residual
    Fstat <- numerator / denominator
    pVal <- pf(Fstat, lmBig$rank - lmSmall$rank, lmBig$df.residual,
    lower.tail = FALSE)
    return(c(num = numerator, den = denominator, Fstat, pVal))
}

myOtherTest <- function(lmHuge, lmBig, lmSmall) {
    numerator <- (deviance(lmSmall) - deviance(lmBig)) /
        (lmBig$rank - lmSmall$rank)
    denominator <- deviance(lmHuge) / lmHuge$df.residual
    Fstat <- numerator / denominator
    pVal <- pf(Fstat, lmBig$rank - lmSmall$rank, lmHuge$df.residual,
    lower.tail = FALSE)
    return(c(num = numerator, den = denominator, Fstat, pVal))
}

deviance(lmFull)
deviance(lmFull) / lmFull$df.resid

deviance(lmNoInter) - deviance(lmFull)
(deviance(lmNoInter) - deviance(lmFull)) / (lmFull$rank -
                                            lmNoInter$rank)
myFtest(lmFull, lmNoInter)

deviance(lmGtype) - deviance(lmNoInter)
(deviance(lmGtype) - deviance(lmNoInter)) / (lmNoInter$rank -
                                             lmGtype$rank)
myFtest(lmNoInter, lmGtype)
myOtherTest(lmFull, lmNoInter, lmGtype)
myOtherTest(lmFull, lmGtype, lmInter)


anova(lmFull, lmNoInter, lmGtype, lmInter)
anova(lmInter, lmGtype, lmNoInter, lmFull)

anova(lmInter, lmGtype)

