library(limma)
library(hexbin)

##########################################################
## load data and design
##########################################################

## load design data.frame
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes) # 39 obs. of  4 variables

## load gene expression data.frame
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)
## 'data.frame':  29949 obs. of  39 variables:

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

##########################################################
## helper functions for finding and featuring genes
##########################################################

source("98_lectSupp-prepareDataFunction.r")
source("99_lectSupp-plotItFunction.r")

##########################################################
## source in-house functions for post-processing mlm objects
##########################################################

## code written by Rick White
source("80_anova-mlm.r")

## was helpful when printing stuff to screen to insert into lecture slides
options(width = 110)

##########################################################
## limma
##########################################################
jDesMat <- model.matrix(~ gType * devStage, prDes)

## make some stuff to show in lecture slides
## ridiculous machination to print a version to screen with small
## variable names
str(jDesMat)
cbind(colnames(jDesMat))
foo <- jDesMat
dimnames(jDesMat)
colnames(foo) <- paste0("X", sprintf("%02d", seq_len(ncol(jDesMat))))
cbind(prDes, foo)

## fitting two-way ANOVA to all probesets at once
system.time(jFit <- lmFit(prDat, jDesMat))
str(jFit, max.level = 3)
summary(jFit)

## can we get hits yet?
topTable(jFit)
topTable(jFit, coef = 2)
topTable(jFit, coef = "gTypeNrlKO")
## NO

ebFit <- eBayes(jFit)
str(ebFit, max.level = 3)
summary(ebFit)

## can we get hits now?
topTable(ebFit)
## YES

##########################################################
## lm for multivariate regression
##########################################################
## prepare to use homegrown mlm code
## responses must be a matrix, one *row* per response

prMat <- t(as.matrix(prDat))
twAnova <- lm(prMat ~ gType * devStage, prDes)
twSumm <- summary(twAnova)
twAnovaNoInter <- lm(prMat ~ gType + devStage, prDes)
twNiSumm <- summary(twAnovaNoInter)
dsOnly <- lm(prMat ~ devStage, prDes)
dsSumm <- summary(dsOnly)

##########################################################
## testing gTypeNrlKO
##########################################################

## this is equivalent: topTable(ebFit, coef = 2)
## but much more self-documenting; USE NAMES!
## you will be glad you did when you revisit/read your code
topTable(ebFit, coef = "gTypeNrlKO")

hits <- topTable(ebFit, coef = "gTypeNrlKO")

stripplotIt(hitDat <- prepareData(head(rownames(hits), 6)),
            scales = list(rot = c(45, 0)))

dev.print(pdf,
          "../figs/limma/stripplot-NrlKO.E16-hits.pdf",
          width = 10, height = 6)

##########################################################
## testing for interaction
##########################################################
colnames(coef(ebFit))                 # parameters in the linear model
grep(":", colnames(coef(ebFit)))      # isolate the interaction terms
hits <- topTable(ebFit, coef = grep(":", colnames(coef(ebFit))))

stripplotIt(hitDat <- prepareData(head(rownames(hits), 6)),
            scales = list(rot = c(45, 0)))

dev.print(pdf,
          "../figs/limma/stripplot-interaction-hits.pdf",
          width = 10, height = 6)

(huh <- rownames(hits)[1])
stripplotIt(huhDat <- prepareData(huh))

huhFitBig <- lm(gExp ~ gType * devStage, huhDat)
huhFitSmall <- lm(gExp ~ gType + devStage, huhDat)
anova(huhFitSmall, huhFitBig)
summary(huhFitBig)

hits <- topTable(ebFit, coef = grep(":", colnames(coef(ebFit))),
                 n = Inf, sort = "none")

foo <- anova(twAnova, twAnovaNoInter)[ , "M1:Pval"]
str(foo)

xyplot(hits$P.Value ~ foo, xlab = "interaction p-value, lm",
       ylab = "interaction p-value, limma", aspect = 1, xbins = 50,
       panel = function(x, y, ...) {
           panel.hexbinplot(x, y, ...)
           panel.abline(a = 0, b = 1, col = "orange")
       })

dev.print(pdf,
          "../figs/limma/smoothScatter-interactionPvalsLimmaVsLm.pdf",
          width = 7, height = 7)

densityplot(~ foo + hits$P.Value, plot.points = FALSE, n = 300, xlab = "p-value")

dev.print(pdf,
          "../figs/limma/densityplot-interactionPvalsLimmaVsLm.pdf",
          width = 7, height = 7)

table(hits$P.Value > foo)
## FALSE  TRUE
## 17702 12247

table(hits$P.Value > foo)/nrow(prDat)
##     FALSE      TRUE
## 0.5910715 0.4089285

sum(foo < 0.01) # 1874
sum(hits$P.Value < 0.01) # 1888

##########################################################
## does genotype matter?
##########################################################
colnames(coef(ebFit))                 # parameters in the linear model
grep("gType", colnames(coef(ebFit)))  # isolate the genotype terms
hits3 <- topTable(ebFit,
                  coef = grep("gType", colnames(coef(ebFit))),
                  n = Inf)

stripplotIt(hitDat <- prepareData(c(head(rownames(hits3), 3),
                                    tail(rownames(hits3), 3))),
            scales = list(rot = c(45, 0)))

dev.print(pdf,
          "../figs/limma/stripplot-gType-hits.pdf",
          width = 10, height = 6)

huh <- twSumm$Full[ , "Pval"]
foo <- anova(twAnova, dsOnly)[ , "M1:Pval"]
str(foo)

hits3 <- topTable(ebFit,
                  coef = grep("gType", colnames(coef(ebFit))),
                  n = Inf, sort = "none")

xyplot(foo ~ hits3$P.Value)
densityplot(~ foo + hits3$P.Value, plot.points = FALSE, n = 300, xlab = "p-value")

## need this for multiple testing lecture to talk about multiple
## testing across contrasts
topTable(ebFit, coef = grep("gType", colnames(coef(ebFit))))

hits3 <- topTable(ebFit,
                  coef = grep("gType", colnames(coef(ebFit))))

##########################################################
## looking at variance
##########################################################

## 'raw' residual variance
str(jFit$sigma)

## df associated with the prior
ebFit$df.prior
## [1] 3.102649

## location associated with the prior
ebFit$s2.prior
## [1] 0.06166446

## post-eb residual variance
str(ebFit$s2.post)

summary(jFit$sigma)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## 0.0704  0.1993  0.2723  0.3305  0.4219  1.4260

summary(jFit$sigma ^ 2)
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
## 0.004956 0.039710 0.074170 0.140700 0.178000 2.033000

table(ebFit$s2.post > jFit$sigma^2)
## FALSE  TRUE
## 17129 12820

table(ebFit$s2.post > jFit$sigma^2)/nrow(prDat)
##    FALSE     TRUE
## 0.571939 0.428061

xyplot(ebFit$s2.post ~ jFit$sigma ^ 2,
       aspect = 1,
       prepanel = function(x, y, ...) {
           rng <- range(x, y, finite = TRUE)
           ## comment next line out to get all vars!
           #rng = c(0, 0.1)
           list(xlim = rng, ylim = rng)
       },
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#          panel.smoothScatter(x, y, ...)
           panel.abline(a = 0, b = 1, col = "grey40")
           with(ebFit,
                panel.abline(a = (df.prior/median(df.total)) * s2.prior,
                             b = (df.residual/df.total),
                             col = "orange"))
      })

dev.print(pdf,
          #"../figs/limma/xyplot-rawAndModeratedResidVarZoom.pdf",
          "../figs/limma/xyplot-rawAndModeratedResidVar.pdf",
          width = 7, height = 7)


str(jFit$sigma)

## df associated with the prior
ebFit$df.prior
## [1] 3.102649

## location associated with the prior
ebFit$s2.prior
## [1] 0.06166446

## post-eb residual variance
str(ebFit$s2.post)


## 2014 refresh has gotten to here ... refresh below as we approach multiple testing lecture

##########################################################
## trying to understand the utility, if any, of a 'sum to zero' parametrization
##########################################################
miniDes <- with(prDes, expand.grid(levels(gType), levels(devStage)))
names(miniDes) <- c("gType", "devStage")

## "the usual"
model.matrix(~ gType * devStage, miniDes)

## "sum to zero"
jDesMat2 <- model.matrix(~ gType * devStage, prDes,
                         contrasts.arg = list(gType = "contr.sum",
                         devStage = "contr.sum"))
system.time(jFit2 <- lmFit(prDat, jDesMat2))
ebFit2 <- eBayes(jFit2)
colnames(coef(ebFit2))                # parameters in the linear model
grep("gType1", colnames(coef(ebFit2)))  # isolate the genotype terms
hits2 <- topTable(ebFit2, coef = grep("gType1",
                          colnames(coef(ebFit2))),
                  n = Inf)
tail(hits2)

colnames(coef(ebFit))                 # parameters in the linear model
grep("gType", colnames(coef(ebFit)))  # isolate the genotype terms
hits3 <- topTable(ebFit,
                  coef = grep("gType", colnames(coef(ebFit))),
                  n = Inf)

head(hits2)
head(hits3)

stripplotIt(hitDat <- prepareData(tail(rownames(hits2), 6)),
            scales = list(rot = c(45, 0)))

##########################################################
## experimenting with contrasts.fit and decideTests for the multiple
## testing lecture
##########################################################

## get my head straight with just 1 probeset
stripplotIt(tinyDat <- prepareData("1448904_at"))
#dev.print(pdf,
 #         file.path(whereLectureLives, "figs/limma/stripplot-1448904_at.pdf"),
  #        width = 7, height = 5)

tinyDat$grp <- with(tinyDat, interaction(gType, devStage))
summary(lm(gExp ~ grp + 0, tinyDat))

## set up design matrix to get "cell means" fit for devStage, within
## gType
prDes$grp <- with(prDes, interaction(gType, devStage))
kDesMat <- model.matrix(~ grp + 0, prDes)

## fitting two-way ANOVA to all probesets at once
kFit <- lmFit(prDat, kDesMat)
colnames(coef(kFit))                  # parameters in the linear model
                                      # = columns of the design matrix
(cont.matrix <- makeContrasts(
    wt.P10vsE16 = grpwt.P10 - grpwt.E16,
    NrlKO.P10vsE16 = grpNrlKO.P10 - grpNrlKO.E16,
    levels = kDesMat))
kFitCont <- contrasts.fit(kFit, cont.matrix)
kebFit <- eBayes(kFitCont)

cutoff <- 1e-4

## use topTable to get separate hits for each of the contrasts
wtHits <-
    topTable(kebFit, coef = "wt.P10vsE16", p.value = cutoff, n = Inf)
NrlKOHits <-
    topTable(kebFit, coef = "NrlKO.P10vsE16", p.value = cutoff, n = Inf)

## same topTable -- same order -- but don't enforce a cutoff
wtHitsAll <- topTable(kebFit, coef = "wt.P10vsE16", n = Inf)
NrlKOHitsAll <- topTable(kebFit, coef = "NrlKO.P10vsE16", n = Inf)

## same topTable -- same order -- but don't enforce a cutoff
wtAllUnsrtd <- topTable(kebFit, coef = "wt.P10vsE16", n = Inf, sort = "none")
NrlKOAllUnsrtd <- topTable(kebFit, coef = "NrlKO.P10vsE16", n = Inf,
                           sort = "none")

nrow(wtHits)                            # 349
nrow(NrlKOHits)                         # 196
nrow(wtHits) + nrow(NrlKOHits)          # 545 total if no overlap
sepHits <- union(wtHits$ID, NrlKOHits$ID)
length(sepHits)                         # 482 actual total

nrow(wtHitsAll)                         # yeah 29949 as expected
## make sure order is same
head(wtHits, 4)
head(wtHitsAll, 4)

nrow(NrlKOHitsAll)                      # yeah 29949 as expected
## make sure order is same
head(NrlKOHits, 4)
head(NrlKOHitsAll, 4)

## show the top hit for wt and NrlKO
stripplotIt(tinyDat <- prepareData(c("1425171_at", "1425530_a_at")))

#dev.print(pdf,
 #         file.path(whereLectureLives, "figs/limma/stripplot-topHitsE16vsP10ByBtype_at.pdf"),
  #        width = 7, height = 5)

## use decideTests, optionally with method = "global"
kRes <- decideTests(kebFit, p.value = cutoff)
kResGlobal <- decideTests(kebFit, p.value = cutoff, method = "global")

summary(kRes)
##    wt.P10vsE16 NrlKO.P10vsE16
## -1         256             52
## 0        29600          29753
## 1           93            144

summary(kResGlobal)
##    wt.P10vsE16 NrlKO.P10vsE16
## -1         246             61
## 0        29615          29735
## 1           88            153

## verifying that the separate hits obtained via two topTable calls
## and via one decideTests call are the same
length(kHitsSep <- which(rowSums(abs(kRes)) > 0)) # 482
all(sort(rownames(kRes[kHitsSep, ])) == sort(sepHits)) # TRUE

peek(kRes[kHitsSep, ])

## how many "global" hits are there?
length(kHitsGlobal <- which(rowSums(abs(kResGlobal)) > 0)) # 485

peek(kResGlobal[kHitsGlobal, ])

## I would have expected the "global" hit list to be *shorter* than
## the separate but apparently my intuition is faulty
## https://stat.ethz.ch/pipermail/bioconductor/2008-July/023593.html
## http://thread.gmane.org/gmane.science.biology.informatics.conductor/17847

## how do these relate to the separate hits?

## first we've got to parse the global hits more finely and track if
## they are hits w.r.t. wt or NrlKO
length(globHitsWt <-
       names(which(abs(kResGlobal[ , "wt.P10vsE16"]) == 1))) # 334
length(globHitsNrlKO <-
       names(which(abs(kResGlobal[ , "NrlKO.P10vsE16"]) == 1))) # 214

## overlap between separate and global results for wt
sum(abs(kResGlobal[ ,"wt.P10vsE16"]))
vennCounts(cbind(separate = wtAllUnsrtd$adj.P.Val < cutoff,
                 global = kResGlobal[ ,"wt.P10vsE16"]),
           include = "both")

## overlap between separate and global results for NrlKO
sum(abs(kResGlobal[ ,"NrlKO.P10vsE16"]))
vennCounts(cbind(separate = NrlKOAllUnsrtd$adj.P.Val < cutoff,
                 global = kResGlobal[ ,"NrlKO.P10vsE16"]),
           include = "both")

## I want to see the global hit status for wtHits
blah <- subset(wtHitsAll, select = c(ID, adj.P.Val))
blah$sepHit <- blah$adj.P.Val < cutoff
blah$globHit <- blah$ID %in% globHitsWt
## checking that sepHit is TRUEs followed by FALSEs
max(which(blah$sepHit))                 # 349
min(which(!blah$sepHit))                # 350
## YES that holds
## checking that globHit is TRUEs followed by FALSEs
max(which(blah$globHit))                # 334
min(which(!blah$globHit))               # 335
## YES that holds
## where does the 'crossover' happen
min(which(blah$sepHit & !blah$globHit)) # 335
max(which(blah$sepHit & !blah$globHit)) # 349
blah[333:353, ]

## I want to see the global hit status for NrlKOHits
blah <- subset(NrlKOHitsAll, select = c(ID, adj.P.Val))
blah$sepHit <- blah$adj.P.Val < cutoff
blah$globHit <- blah$ID %in% globHitsNrlKO
## checking that sepHit is TRUEs followed by FALSEs
max(which(blah$sepHit))                 # 196
min(which(!blah$sepHit))                # 197
## YES that holds
## checking that globHit is TRUEs followed by FALSEs
max(which(blah$globHit))                # 214
min(which(!blah$globHit))               # 215
## YES that holds
## where does the 'crossover' happen
min(which(blah$sepHit & !blah$globHit)) # oops! does NOT happen
min(which(!blah$sepHit & blah$globHit)) # 197
max(which(!blah$sepHit & blah$globHit)) # 214
blah[195:217, ]

## combine all hit info together so I can choose some genes to
## spotlight
blah <-
    data.frame(subset(wtAllUnsrtd, select = c(ID, adj.P.Val)))
blah$sepHit <- blah$adj.P.Val < cutoff
blah$globHit <- blah$ID %in% globHitsWt
names(blah) <- paste0("wt.", names(blah))
blahblah <-
    data.frame(subset(NrlKOAllUnsrtd, select = c(ID, adj.P.Val)))
blahblah$sepHit <- blahblah$adj.P.Val < cutoff
blahblah$globHit <- blahblah$ID %in% globHitsNrlKO
names(blahblah) <- paste0("NrlKO.", names(blahblah))
jRes <- merge(blah, blahblah, by.x = "wt.ID", by.y = "NrlKO.ID")
str(jRes)
with(jRes,
     table(wt.sepHit, wt.globHit, NrlKO.sepHit, NrlKO.globHit))

## picking probes to feature
set.seed(123)
(boring <- with(jRes,
                sample(which(!wt.sepHit & !wt.globHit &
                             !NrlKO.sepHit & ! NrlKO.globHit),
                       size = 1)))      # 8620

(yesYesYes <- with(jRes,
                   sample(which(wt.sepHit & wt.globHit &
                                NrlKO.sepHit & NrlKO.globHit),
                          size = 1)))      # 22171

(NrlKOdisc <- with(jRes,
                   sample(which(!wt.sepHit &
                                !NrlKO.sepHit & NrlKO.globHit),
                          size = 1)))      # 15969

(wtdisc <- with(jRes,
                sample(which(wt.sepHit & !wt.globHit &
                !NrlKO.globHit & !NrlKO.globHit),
                       size = 1)))      # 28080

getMe <- rownames(prDat)[c(boring, yesYesYes, NrlKOdisc, wtdisc)]
#getMe <- rownames(prDat)[yesYesYes]
hitDat <- prepareData(getMe)

stripplot(gExp ~ devStage | gene, hitDat,
          group = gType, subset = devStage %in% c("E16", "P10"),
          jitter.data = TRUE,
          auto.key = TRUE, type = c('p', 'a'), grid = TRUE,
          scales = list(x = list(rot = c(45, 0)),
          y = list(relation = "free")))

dev.print(pdf,
          file.path(whereLectureLives, "figs/limma/stripplot-fourExamplesE16vsP10ByBtype_at.pdf"),
          width = 10, height = 6)

##########################################################
## inspecting results for Nrl = 1450946_at
##########################################################

##########################################################
## trying to understand limma manual around page 43 (different ways to
## do 2x2 factorial
##########################################################

limdat <- data.frame(strain = factor(c("wt", "wt", "mut", "mut", "mut"),
                     levels = c("wt", "mut")),
                     tx = factor(c("untx", "stim", "untx", "stim",
                     "stim"),
                     levels = c("untx", "stim")))

model.matrix(~ strain + strain:tx, limDat)
##   (Intercept) strainmut strainwt:txstim strainmut:txstim
## 1           1         0               0                0
## 2           1         0               1                0
## 3           1         1               0                0
## 4           1         1               0                1
## 5           1         1               0                1

model.matrix(~ strain/tx, limDat)
##   (Intercept) strainmut strainwt:txstim strainmut:txstim
## 1           1         0               0                0
## 2           1         0               1                0
## 3           1         1               0                0
## 4           1         1               0                1
## 5           1         1               0                1

## Yep, the same result.

(huh <- rownames(hits)[1])
stripplotIt(huhDat <- prepareData(huh))

huhFitBig <- lm(gExp ~ gType * devStage, huhDat)
huhFitSmall <- lm(gExp ~ gType + devStage, huhDat)
anova(huhFitSmall, huhFitBig)

## Yes these are the same ... I trust
head(cbind(coef(ebFitGtypeStarDevStage),
           coef(fitGtypeStarDevStage)))

## this is the mean of the inverse Chisquare prior for the
## gene-specific variances
ebFitGtypeStarDevStage$s2.prior
## [1] 0.06166446

## this is the degrees of freedom associated with the inverse
## ChiSquare prior for gene-specific variances
ebFitGtypeStarDevStage$df.prior
## [1] 3.102649



library(preprocessCore)
s1pNorm <- as.data.frame(normalize.quantiles(as.matrix(s1pDat)))
dimnames(s1pNorm) <- dimnames(s1pDat)
str(s1pNorm, list.len = 6)
s1pNormTall <- # this takes a moment; horribly inefficient
  data.frame(sample = factor(rep(rownames(s1pNorm), each = nrow(s1pNorm)),
                             levels = rownames(s1pNorm)),
             gExp = unlist(s1pNorm))
