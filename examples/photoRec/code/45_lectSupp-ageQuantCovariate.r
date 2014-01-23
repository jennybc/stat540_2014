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

## handy for scatterplotting the data for a small set of probesets
graphIt2 <- function(jDat, yFree = FALSE) {
    thePlot <- xyplot(gExp ~ age | gene, jDat,
                      type = c('p', 'a'), grid = TRUE,
                      group = gType, auto.key = TRUE,
                      jitter.data = TRUE)
    if(yFree) {
        thePlot <- update(thePlot,
                          scales = list(y = list(relation = "free")))
    }
    print(thePlot)
}

## handy for scatterplotting the data for a small set of probesets
graphIt3 <- function(jDat, yFree = FALSE) {
    thePlot <- xyplot(gExp ~ age | gene, jDat,
                      type = c('p', 'r'), grid = TRUE,
                      subset = gType == "wt",
                      jitter.data = TRUE)
    if(yFree) {
        thePlot <- update(thePlot,
                          scales = list(y = list(relation = "free")))
    }
    print(thePlot)
}


##########################################################
## finding good exemplars ... line, first try
##########################################################

## I want a gene where linear regression on age will work nicely
rFit <- lm(prMat ~ devStage, subset = gType == "wt")
rfSumm <- summary(rFit)
str(rfSumm$Coef)
## isolate the non-intercept, non-4_weeks param ests
str(rfSumm$Coef[ , 2:4, "Estimate"])
litmus <-
    data.frame(maxPval = apply(rfSumm$Coef[ , 2:4, "Pr(>|t|)"], 1, max),
               magnitude = apply(abs(rfSumm$Coef[ , 2:4, "Estimate"]), 1, min),
               sameSign = abs(rowSums(sign(rfSumm$Coef[ , 2:4,
                                                       "Estimate"]))) == 4,
               range = apply(rfSumm$Coef[ , 2:4, "Estimate"], 1,
                             function(z) max(z) - min(z)))
summary(litmus)
with(litmus,
     sum(getMe <- maxPval < 0.05 & sameSign &
         magnitude > quantile(magnitude, 0.95))) # 1028

set.seed(999)
graphIt2(extractIt(yo <- sample(which(foo), size = 12)))

## interesting
## "1424963_at"
## "1437687_x_at"
## "1439571_at"
## "1423506_a_at"
## "1425222_x_at"
## "1436552_at"
(yo <- c("1424963_at", "1437687_x_at", "1439571_at", "1423506_a_at",
         "1425222_x_at", "1436552_at"))
graphIt2(extractIt(yo))

##########################################################
## finding good exemplars ... line, second try
##########################################################

## I want a gene where linear regression on age will work nicely
rFit <- lm(prMat ~ age, subset = gType == "wt")
rfSumm <- summary(rFit)
str(rfSumm$Coef)
## isolate the non-intercept, non-4_weeks param ests
str(rfSumm$Coef[ , "age", "Estimate"])
litmus <- data.frame(rfSumm$Coef[ , "age", ])
names(litmus) <- c("est","se","t","pval")
summary(litmus)
sum(getMe <-
    with(litmus,
         pval < 0.05 &
         est > quantile(est, 0.95))) # 1464

set.seed(999)
graphIt3(extractIt(yo <- sample(which(foo), size = 12)))

## interesting
## "1441811_x_at"

##########################################################
## finding good exemplars ... parabola
##########################################################

## I want a gene where quadratic regression on age will work nicely
rFit <- lm(prMat ~ age + I(age^2), subset = gType == "wt")
rfSumm <- summary(rFit)
str(rfSumm$Coef)
## isolate the non-intercept param ests
str(rfSumm$Coef[ , -1, "Estimate"])
litmus <-
    data.frame(quadPval = rfSumm$Coef[ , "I(age^2)", "Pr(>|t|)"],
               quadMagnitude = abs(rfSumm$Coef[ , "I(age^2)", "Estimate"]))
summary(litmus)
sum(getMe <- with(litmus,
                  quadPval < 0.05 &
                  quadMagnitude > quantile(quadMagnitude, 0.9))) # 9550
set.seed(222)
graphIt2(extractIt(yo <- sample(which(getMe), size = 12)))

## fairly steady decrease in expression
## practically no difference in wt and NrlKO
## "1429688_at"
## "1415840_at"
## "1419740_at"
## "1457682_at"
## "1427275_at"
## "1456341_a_at"
## "1449519_at"
## "1419864_x_at"
## "1455972_x_at"

(yo <- c("1429688_at", "1415840_at", "1419740_at", "1457682_at",
         "1427275_at", "1456341_a_at", "1449519_at", "1419864_x_at",
         "1455972_x_at"))
graphIt2(extractIt(yo))

##########################################################
## finding good exemplars ... weird
##########################################################

## I want a gene where neither linear nor quadratic regression will
## work well
rFit <- lm(prMat ~ devStage, subset = gType == "wt")
rfSumm <- summary(rFit)
str(rfSumm$Coef)
## isolate the non-intercept param ests
str(rfSumm$Coef[ , 2:4, "Estimate"])
litmus <-
    data.frame(maxPval = apply(rfSumm$Coef[ , 2:4, "Pr(>|t|)"], 1,
               max),
               minPval = apply(rfSumm$Coef[ , 2:4, "Pr(>|t|)"], 1, min),
               magnitude = apply(abs(rfSumm$Coef[ , 2:4, "Estimate"]), 1, min),
               signSum = rowSums(sign(rfSumm$Coef[ , 2:4, "Estimate"])),
               range = apply(rfSumm$Coef[ , 2:4, "Estimate"], 1,
               function(z) max(z) - min(z)))

summary(litmus)
sum(getMe <- with(litmus,
                  minPval < 0.025 &
                  signSum >= -1 &
                  signSum <= 1))        # 413
set.seed(222)
graphIt2(extractIt(yo <- sample(which(getMe), size = 12)))

## interestingly weird
## "1456341_a_at"
## "1444790_at"

(yo <- c("1456341_a_at", "1444790_at"))
graphIt2(extractIt(yo))

##########################################################
## check on my exemplars
##########################################################
(luckyGenes <- c("1438815_at",          # selected from prev lect
                 ## found in the 'linear' search, part 1
                 "1424963_at", "1437687_x_at",
                 ## just need to narrow things down for viewing!
#                 "1439571_at", "1423506_a_at",
                 "1425222_x_at", "1436552_at",
                 "1441811_x_at"
                 ## found in the 'quadratic' search
                 "1429688_at", "1415840_at", "1419740_at", "1457682_at",
                 "1427275_at", "1456341_a_at", "1449519_at", "1419864_x_at",
                 "1455972_x_at",
                 ## found in the 'weird' search
                 "1456341_a_at", "1444790_at"))
#(luckyGene <- "1452702_at")

(luckyGenes <- c("1427275_at", "1456341_a_at", "1441811_x_at"))

jDat <- extractIt(luckyGenes)

xyplot(gExp ~ age | gene, jDat,
       group = gType, grid = TRUE, type = c('p', 'a'))

linFits <- with(jDat,
                by(jDat, gene, function(z) {
                    lm(gExp ~ age, z, subset = gType == 'wt')
                }))
quadFits <- with(jDat,
                 by(jDat, gene, function(z) {
                     lm(gExp ~ age + I(age^2), z, subset = gType == 'wt')
                 }))

nPts <- 100
linFitted <- lapply(linFits, function(theFit) {
    newDat <-
        with(prDes,
             data.frame(age = seq(from = min(age), to = max(age), length = nPts)))
    newDat$fitted <- predict(theFit, newDat)
    newDat
})
quadFitted <- lapply(quadFits, function(theFit) {
    newDat <-
        with(prDes,
             data.frame(age = seq(from = min(age), to = max(age), length = nPts)))
    newDat$fitted <- predict(theFit, newDat)
    newDat
})

xyplot(gExp ~ age | gene, jDat,
       subset = gType == "wt", grid = TRUE,
#       group = gType, grid = TRUE,
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           with(linFitted[[packet.number()]],
                panel.lines(age, fitted))
       })

xyplot(gExp ~ age | gene, jDat,
       subset = gType == "wt", grid = TRUE,
#       group = gType, grid = TRUE,
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           with(quadFitted[[packet.number()]],
                panel.lines(age, fitted))
       })

##########################################################
## simple linear regression on age for the exemplars
## wild type data only
##########################################################

jDat <- extractIt(luckyGenes)
jDat <- droplevels(subset(jDat, gType == "wt", select = -gType))
str(jDat)

print(jStrip <- stripplot(gExp ~ devStage | gene, jDat,
                    type = c('p', 'a'), grid = TRUE, layout = c(3, 1)))

print(jScatter <- xyplot(gExp ~ age | gene, jDat,
                   type = c('p', 'a'), grid = TRUE, layout = c(3, 1)))

print(jStrip, c(0, 0, 1, 0.5), more = TRUE)
print(jScatter, c(0, 0.5, 1, 1), more = FALSE)

dev.print(pdf,
          paste0(whereLectureLives, "figs/ageQuantCovariate/exemplars-wtOnly-devStageVsAge.pdf"),
          width = 10, height = 7)

linFits <- with(jDat, by(jDat, gene, function(z) lm(gExp ~ age, z)))

## there are easier ways to put simple linear regression lines down
## but they don't generalize well, so I'm doing it the 'hard' way
nPts <- 100
linFitted <- lapply(linFits, function(theFit) {
    newDat <-
        with(prDes,
             data.frame(age = seq(from = min(age), to = max(age), length = nPts)))
    newDat$fitted <- predict(theFit, newDat)
    newDat
})

xyplot(gExp ~ age | gene, jDat,
       grid = TRUE, layout = c(3, 1),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           with(linFitted[[packet.number()]],
                panel.lines(age, fitted))
       })

dev.print(pdf,
          paste0(whereLectureLives, "figs/ageQuantCovariate/exemplars-wtOnly-scatterAndLine.pdf"),
          width = 10, height = 4)

summary(linFits[["1441811_x_at"]])

##########################################################
## quadratic model for the exemplars
## wild type data only
##########################################################

jDat <- extractIt(luckyGenes)
jDat <- droplevels(subset(jDat, gType == "wt", select = -gType))
str(jDat)

quadFits <- with(jDat, by(jDat, gene,
                          function(z) lm(gExp ~ age + I(age^2), z)))

## there are easier ways to put simple linear regression lines down
## but they don't generalize well, so I'm doing it the 'hard' way
nPts <- 100
quadFitted <- lapply(quadFits, function(theFit) {
    newDat <-
        with(prDes,
             data.frame(age = seq(from = min(age), to = max(age), length = nPts)))
    newDat$fitted <- predict(theFit, newDat)
    newDat
})

xyplot(gExp ~ age | gene, jDat,
       grid = TRUE, layout = c(3, 1),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           with(quadFitted[[packet.number()]],
                panel.lines(age, fitted))
       })

dev.print(pdf,
          paste0(whereLectureLives, "figs/ageQuantCovariate/exemplars-wtOnly-scatterAndParabola.pdf"),
          width = 10, height = 4)

summary(quadFits[["1427275_at"]])
summary(quadFits[["1441811_x_at"]])
summary(quadFits[["1456341_a_at"]])

##########################################################
## compare models for the exemplars
## wild type data only
##########################################################

factFits <- with(jDat, by(jDat, gene, function(z) lm(gExp ~ devStage, z)))

(jGene <- luckyGenes[1])
(jGene <- luckyGenes[2])
(jGene <- luckyGenes[3])

anova(linFits[[jGene]], quadFits[[jGene]])
anova(quadFits[[jGene]], factFits[[jGene]])

print(jLinear <-
      xyplot(gExp ~ age, jDat,
             subset = gene == jGene,
             grid = TRUE,
             panel = function(x, y, ...) {
                 panel.xyplot(x, y, ...)
                 with(linFitted[[jGene]],
                      panel.lines(age, fitted))
             }))


print(jQuadratic <-
      xyplot(gExp ~ age, jDat,
             subset = gene == jGene,
             grid = TRUE,
             panel = function(x, y, ...) {
                 panel.xyplot(x, y, ...)
                 with(quadFitted[[jGene]],
                      panel.lines(age, fitted))
             }))

print(jAnova <-
      xyplot(gExp ~ age, jDat,
             subset = gene == jGene,
             grid = TRUE, type = c('p','a')))

print(jAnova, c(0, 0, 1, 0.33), more = TRUE)
print(jQuadratic, c(0, 0.33, 1, 0.67), more = TRUE)
print(jLinear, c(0, 0.67, 1, 1), more = FALSE)

dev.print(pdf,
          paste0(whereLectureLives, "figs/ageQuantCovariate/", jGene, "-wtOnly-linearVsQuad.pdf"),
          width = 6, height = 10)

##########################################################
## now put knockouts back in, i.e. use gType
## just the linear term
##########################################################

jDat <- extractIt(luckyGenes)
str(jDat)

xyplot(gExp ~ age | gene, jDat,
       group = gType, grid = TRUE, layout = c(3, 1),
       type = c('p', 'a'))

dev.print(pdf,
          paste0(whereLectureLives, "figs/ageQuantCovariate/exemplars-scatterAndMeans.pdf"),
          width = 10, height = 5)

## let's focus on "1441811_x_at" for linear
jDat <- extractIt("1441811_x_at")
str(jDat)

xyplot(gExp ~ age | gene, jDat,
       group = gType, grid = TRUE,
       type = c('p', 'r'))

dev.print(pdf,
          paste0(whereLectureLives, "figs/ageQuantCovariate/1441811_x_at-scatterAndLines.pdf"),
          width = 7, height = 5)

jFit <- lm(gExp ~ gType * age, jDat)
summary(jFit)

jFitAlt <- lm(gExp ~ gType/age - 1, jDat)
summary(jFitAlt)

coef(jFit)
coef(jFitAlt)

(contMat <- rbind(c(1, 0, 0, 0),
                  c(1, 1, 0, 0),
                  c(0, 0, 1, 0),
                  c(0, 0, 1, 1)))

cbind(coefDefault = coef(jFit),
      coefAlt = coef(jFitAlt),
      matrixResult = as.vector(contMat %*% coef(jFit)))

anova(lm(gExp ~ age, jDat), jFit)

##########################################################
## now put knockouts back in, i.e. use gType
## include quadratic term
##########################################################

## let's focus on "1427275_at" for quadratic
jDat <- extractIt("1427275_at")
str(jDat)

jFit <- lm(gExp ~ gType * (age + I(age^2)), jDat)
summary(jFit)

nPts <- 100
quadFitted <-
    with(jDat,
         data.frame(age = seq(from = min(age), to = max(age),
                    length = nPts),
                    gType = factor(rep(levels(gType), each = nPts),
                    levels(gType))))
quadFitted$fitted <- predict(jFit, quadFitted)

xyplot(gExp ~ age | gene, jDat,
       group = gType, grid = TRUE,
       panel = panel.superpose,
       panel.group = function(x, y, ..., group.number) {
           panel.xyplot(x, y, ...)
           with(quadFitted,
                panel.lines(age[gType == levels(gType)[group.number]],
                            fitted[gType == levels(gType)[group.number]],
                            col =
                            trellis.par.get("superpose.line")$col[group.number]))
       })


dev.print(pdf,
          paste0(whereLectureLives, "figs/ageQuantCovariate/1427275_at-scatterAndParabolas.pdf"),
          width = 7, height = 5)

anova(lm(gExp ~ age + I(age^2), jDat),
      lm(gExp ~ gType * (age + I(age^2)), jDat))


##########################################################
## code below here is rough, not currently used but I'm afraid to delete
##########################################################

rFactFit <- lm(prMat ~ gType * devStage)
rQuadFullFit <- lm(prMat ~ gType * (age + I(age^2)))
rQuadRestrFit <- lm(prMat ~ age + I(age^2))

rffSumm <- summary(rFactFit)

rqfSumm <- summary(rQuadFullFit)
rqrSumm <- summary(rQuadRestrFit)

##########################################################
## testing for interaction
##########################################################

rfAnova <- anova(rQuadFullFit, rQuadRestrFit)
plot(rfAnova)                           # things whizz by ....
rfAnova                                 # first 2 fits, by default

peek(rfAnova)

##########################################################
## make it easy to find various types of 'hits'
##########################################################

## get p-values for the the main effects and interaction
pVals <- cbind(interaction = rfAnova[ , "M1:Pval"],
               rqfSumm$Coef[ , c("gTypeNrlKO:I(age^2)",
                                 "gTypeNrlKO:age",
                                 "gTypeNrlKO",
                                 "I(age^2)",
                                 "age"), "Pr(>|t|)"])
splom(~ pVals, panel = panel.smoothScatter)
pValsTall <- data.frame(pVals = as.vector(pVals),
                        what = rep(colnames(pVals),
                        each = nrow(pVals)))
densityplot(~ pVals, pValsTall,
            group = what, auto.key = TRUE,
            plot.points = FALSE)


##########################################################
## finding some genes to feature in lecture
##########################################################

## boring
getMe <- which(rowSums(pVals < 0.05) == 0)
length(getMe)                           # 11244
set.seed(111)
(boring <- sample(getMe, size = 3))
graphIt(extractIt(boring))
print(rqfSumm, show = which(colnames(prMat) %in% names(boring)))
print(rqfSumm, show = which(colnames(prMat) == names(boring)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(boring)))

#dev.print(pdf,
#          paste0(whereLectureLives, "figs/gTypeByDevStage/stripplot-boring.pdf"),
#          width = 7, height = 7)

## age YES NrlKO NO interaction NO
getMe <- which(pVals[ , "interaction"] > 0.2 &
               (pVals[ , "age"] < 0.05 | pVals[ , "I(age^2)"] < 0.05) &
               (pVals[ , "gTypeNrlKO"] > 0.15 &
                pVals[ , "gTypeNrlKO:age"] > 0.15 &
                pVals[ , "gTypeNrlKO:I(age^2)"] > 0.15))
length(getMe)                           # 3550
set.seed(111)
(ageOnly <- sample(getMe, size = 3))
graphIt(extractIt(ageOnly))
print(rqfSumm, show = which(colnames(prMat) %in% names(ageOnly)))
print(rqfSumm, show = which(colnames(prMat) == names(ageOnly)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(ageOnly)))

## age NO NrlKO YES interaction NO
getMe <- which(pVals[ , "interaction"] > 0.2 &
               pVals[ , "age"] > 0.20 & pVals[ , "I(age^2)"] > 0.2 &
               pVals[ , "gTypeNrlKO"] < 0.05)
length(getMe)                           # 28
set.seed(111)
(koOnly <- sample(getMe, size = 3))
graphIt(extractIt(koOnly))
print(rqfSumm, show = which(colnames(prMat) %in% names(koOnly)))
print(rqfSumm, show = which(colnames(prMat) == names(koOnly)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(koOnly)))

## age YES NrlKO YES interaction NO
getMe <- which(pVals[ , "interaction"] > 0.2 &
               (pVals[ , "age"] < 0.05 | pVals[ , "I(age^2)"] < 0.05) &
               pVals[ , "gTypeNrlKO"] < 0.05)
length(getMe)                           # 90
set.seed(111)
(mainEffOnly <- sample(getMe, size = 3))
graphIt(extractIt(mainEffOnly))
print(rqfSumm, show = which(colnames(prMat) %in% names(mainEffOnly)))
print(rqfSumm, show = which(colnames(prMat) == names(mainEffOnly)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(mainEffOnly)))

## devStage YES NrlKO YES interaction YES
getMe <- which(pVals[ , "interaction"] < 0.05 &
               (pVals[ , "age"] < 0.05 | pVals[ , "I(age^2)"] < 0.05) &
               pVals[ , "gTypeNrlKO"] < 0.05)
length(getMe)                           # 1771
set.seed(111)
(exciting <- sample(getMe, size = 3))
graphIt(extractIt(exciting))
print(rqfSumm, show = which(colnames(prMat) %in% names(exciting)))
print(rqfSumm, show = which(colnames(prMat) == names(exciting)[2]))
print(rfAnova, show = which(colnames(prMat) %in% names(exciting)))

##########################################################
## writing results in case I want to mine this again soon
##########################################################

## not edited for this particular diff exp anal yet!

## deGtypeByDevStage <-
##     data.frame(coef(ebFitGtypeStarDevStage), # est coefs
##                rfSumm$Coef[ , , "Std. Error"], # est se's
##                rfSumm$Coef[ , , "t value"], # t stats
##                rfSumm$Coef[ , , "Pr(>|t|)"], # p vals
##                deviance(rFit)/df.residual(rFit), # resid var
##                rfSumm$FullModelStats[ , "Fstat"], # overall F test
##                rfSumm$FullModelStats[ , "Pval"],  # overall F pval
##                rfAnova[ , c("gType:Pval", "devStage:Pval", # anova
##                             "gType:devStage:Pval")]) # table stuff
## str(deGtypeByDevStage)
## peek(deGtypeByDevStage)

## yo <- c(rep(c("est", "se", "tStat", "pVal"), each = 10))
## yo2 <- c(rep(c("int", "NrlKO", "4_weeks", "NrlKO:4_weeks"), 10))
## colnames(deGtypeByDevStage) <- c(paste(yo, yo2, sep = "."),
##                                  c("sigSq", "overallFstat",
##                                    "overallFp",
##                                    "gTypeMainEffPval",
##                                    "devStageMainEffPval",
##                                    "interactionPval"))
## peek(deGtypeByDevStage)

## ## good news: writes a plain text version to file, easy to read
## ## by both humans and machines (e.g. Excel)
## ## bad news: order of factor levels not captured
## write.table(deGtypeByDevStage,
##             file = file.path(whereAmI,
##             "rmd/data/photoRec/deGtypeDevStage.txt"),
##             quote = FALSE, row.names = FALSE)

## ## good news: writes an easy-to-reload-in-R version that
## ## captures order of factor levels
## ## bad news: the file is binary and totally specific to R
## save(deGtypeByDevStage,
##      file = file.path(whereAmI,
##      "rmd/data/photoRec/deGtypeByDevStage.robj"))





## factFitted <-
##     with(prDes,
##          data.frame(gType = factor(rep(levels(gType), each = nlevels(devStage)),
##                     levels = levels(gType)),
##                     devStage = factor(levels(devStage),
##     levels(devStage))))
## factFitted$fitted <- predict(factFit, factFitted)

## jStrip <-
##     stripplot(gExp ~ devStage, jDat, group = gType, grid = TRUE,
##               panel = panel.superpose,
##               panel.groups = function(x, y , ..., group.number) {
##                   panel.stripplot(x, y, ...)
##                   panel.lines(seq_len(nlevels(prDes$devStage)),
##                               factFitted$fitted[factFitted$gType ==
##                                                 levels(factFitted$gType)[group.number]],
##                               col = trellis.par.get("superpose.line")$col[group.number])
##               })
