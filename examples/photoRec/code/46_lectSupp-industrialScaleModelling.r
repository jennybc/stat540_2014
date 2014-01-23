library(car)                            # recode()

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
## adding quantitative version of devStage
##########################################################
## recode() is from add-on package 'car'
prDes$age <-
    recode(prDes$devStage,
           "'E16'=-2; 'P2'=2; 'P6'=6; 'P10'=10; '4_weeks'=28",
           as.factor.result = FALSE)
peek(prDes)
str(prDes)

##########################################################
## prepare to use homegrown mlm code
##########################################################

## responses must be a matrix, one *row* per response
prMat <- t(as.matrix(prDat))
gType <- prDes$gType              # lesser of two evils
devStage <- prDes$devStage        # lesser of two evils
age <- prDes$age                  # lesser of two evils

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

## handy for scatterplotting the data for a small set of probesets
graphIt <- function(jDat, yFree = FALSE, ...) {
    thePlot <- stripplot(gExp ~ devStage | gene, jDat,
                      type = c('p', 'a'), grid = TRUE,
                      group = gType, auto.key = TRUE,
                      jitter.data = TRUE, ...)
    if(yFree) {
        thePlot <- update(thePlot,
                          scales = list(y = list(relation = "free")))
    }
    print(thePlot)
}

## handy for scatterplotting the data for a small set of probesets
graphIt2 <- function(jDat, yFree = FALSE, ...) {
    thePlot <- xyplot(gExp ~ age | gene, jDat,
                      type = c('p', 'a'), grid = TRUE,
                      group = gType, auto.key = TRUE,
                      jitter.data = TRUE, ...)
    if(yFree) {
        thePlot <- update(thePlot,
                          scales = list(y = list(relation = "free")))
    }
    print(thePlot)
}

## handy for scatterplotting the data for a small set of probesets
graphIt3 <- function(jDat, yFree = FALSE, ...) {
    thePlot <- xyplot(gExp ~ age | gene, jDat,
                      type = c('p', 'r'), grid = TRUE,
                      subset = gType == "wt",
                      jitter.data = TRUE, ...)
    if(yFree) {
        thePlot <- update(thePlot,
                          scales = list(y = list(relation = "free")))
    }
    print(thePlot)
}

##########################################################
## re-fitting some of the models we've worked with
##########################################################

twAnova <- lm(prMat ~ gType * devStage)
dsOnly <- lm(prMat ~ devStage)
twAnovaNoInter <- lm(prMat ~ gType + devStage)
quadBoth <- lm(prMat ~ gType * (age + I(age^2)))
quadAge <- lm(prMat ~ age + I(age^2))

twSumm <- summary(twAnova)
dsSumm <- summary(dsOnly)
twNiSumm <- summary(twAnovaNoInter)
qBothSumm <- summary(quadBoth)
qAgeSumm <- summary(quadAge)

##########################################################
## error variance
##########################################################
str(twSumm, max.level = 1)
str(twSumm$Full)
str(twSumm$Full[ , "Res. SE"])

stdDevn <-
    data.frame(twAnova = twSumm$Full[ , "Res. SE"],
               dsOnly = dsSumm$Full[ , "Res. SE"],
               quadBoth = qBothSumm$Full[ , "Res. SE"],
               quadAge = qAgeSumm$Full[ , "Res. SE"])
summary(stdDevn)
sdTall <-
    data.frame(fittedModel =
               factor(rep(names(vars), each = nrow(vars)),
                      levels = names(vars)),
               estResSd = unlist(stdDevn))

splom(stdDevn,
      prepanel.limits = function(x) c(0, 2),
      panel = function(x, y, ... ) {
          panel.smoothScatter(x, y, ...)
          panel.abline(a = 0, b = 1, col = "orange")
      })

dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/splom-errorStdDevn.pdf"),
          width = 7, height = 7)

densityplot(~ estResSd, sdTall,
            group = fittedModel, auto.key = TRUE,
            plot.points = FALSE, n = 300)

dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/densityplot-errorStdDevn.pdf"),
          width = 7, height = 7)

sum(highVar <- with(stdDevn,
                    twAnova > quantile(twAnova, 0.95))) # 1498

sum(lowVar <- with(stdDevn,
                    twAnova < quantile(twAnova, 0.05))) # 1498

set.seed(111)
graphIt2(extractIt(yo <- c(sample(which(highVar), size = 4),
                           sample(which(lowVar), size = 4))),
         layout = c(4, 2))
dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/scatterplot-lowHighVar.pdf"),
          width = 7, height = 7)

##########################################################
## interaction between gType and devStage, ANOVA
##########################################################
huh <- twSumm$Full[ , "Pval"]
foo <- anova(twAnova, twAnovaNoInter)[ , "M1:Pval"]
str(foo)

densityplot(~ foo, plot.points = FALSE, n = 300, xlab = "p-value")

dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/densityplot-interactionYesAnova.pdf"),
          width = 7, height = 7)

sum(interactionYES <- foo < 0.005 & huh < 0.005) # 1344
sum(interactionNO <- foo > 0.9 & huh < 0.005)    # 155

set.seed(333)
graphIt(extractIt(yo <- c(sample(which(interactionYES), size = 4),
                           sample(which(interactionNO), size = 4))),
         layout = c(4, 2))

dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/stripplot-interactionYesNoAnova.pdf"),
          width = 7, height = 7)

##########################################################
## does gType matter? ANOVA
##########################################################
huh <- twSumm$Full[ , "Pval"]
foo <- anova(twAnova, dsOnly)[ , "M1:Pval"]
str(foo)

densityplot(~ foo, plot.points = FALSE, n = 300, xlab = "p-value")

dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/densityplot-gTypeMattersAnova.pdf"),
          width = 7, height = 7)

sum(gTypeYES <- foo < 0.005 & huh < 0.005) # 2629

set.seed(333)
graphIt(extractIt(yo <- sample(which(gTypeYES), size = 8)),
        layout = c(4, 2))

dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/stripplot-gTypeYesAnova.pdf"),
          width = 7, height = 7)

##########################################################
## overall significance of model
##########################################################
tw <- twSumm$Full[ , "Pval"]
qb <- qBothSumm$Full[ , "Pval"]

densityplot(~ tw + qb, plot.points = FALSE, n = 300, xlab = "p-value",
            auto.key = TRUE)

dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/densityplot-overallFAnovaVsQuad.pdf"),
          width = 7, height = 7)

xyplot(tw ~ qb,
       panel = function(x, y, ...) {
           panel.smoothScatter(x, y, nbin = 200, raster = TRUE, ...)
           panel.abline(a = 0, b = 1, col = "orange")
           panel.smooth(x, y, col = "pink")
       })

           type = c('p', 'smooth'))

dev.print(pdf,
          paste0(whereLectureLives, "figs/enMasse/stripplot-gTypeYesAnova.pdf"),
          width = 7, height = 7)





