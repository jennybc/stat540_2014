library(lattice)

## load photoRec data and design from the interwebs to delay filepath pain
## use the dput/dget workflow since works better with URLs 
jURL <- "http://www.ugrad.stat.ubc.ca/~stat540/examples/photoRec/data/GSE4051_design_DPUT.txt"
str(prDes <- dget(jURL))
jURL <- "http://www.ugrad.stat.ubc.ca/~stat540/examples/photoRec/data/GSE4051_data.tsv"
str(prDat <- read.table(jURL), max.level = 0)

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

## I've already done differential expression analysis for gType * a
## simplified version of devStage, including only E16 and 4_weeks

## I've selected this as our example
(luckyGene <- which(rownames(prDat) == "1455695_at")) # 26861
miniDat <- data.frame(gExp = unlist(prDat[luckyGene, ]))
miniDat <- data.frame(prDes, miniDat)
miniDat <- droplevels(subset(miniDat,
                             devStage %in%
                             levels(devStage)[c(1, nlevels(devStage))]))
miniDat$grp <- with(miniDat, interaction(gType, devStage))
str(miniDat)
with(miniDat, table(gType, devStage))
table(miniDat$grp)

stripplot(gExp ~ devStage, miniDat,
          group = gType, auto.key = TRUE,
          grid = TRUE, type = c('p', 'a'),
          jitter.data = TRUE)

dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-twoFact.pdf",
          width = 7, height = 4)

stripplot(gExp ~ grp, miniDat,
          grid = TRUE, type = c('p', 'a'),
          jitter.data = TRUE)

dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-grp.pdf",
          width = 7, height = 4)

(theAvgs <- with(miniDat, tapply(gExp, grp, mean)))

with(miniDat, tapply(gExp, list(gType, devStage), mean))

twoFactFit <- lm(gExp ~ gType * devStage, miniDat)
grpFit <- lm(gExp ~ grp, miniDat)

summary(grpFit)
cbind(sampleMeans = theAvgs,
      minuRef = theAvgs - theAvgs["wt.E16"],
      grpFit = coef(grpFit))

summary(twoFactFit)
cbind(sampleMeans = theAvgs,
      minuRef = theAvgs - theAvgs["wt.E16"],
      twoFactFit = coef(twoFactFit))

theAvgs["NrlKO.4_weeks"] -
    (theAvgs["wt.E16"] +
     (theAvgs["NrlKO.E16"] - theAvgs["wt.E16"]) +
      (theAvgs["wt.4_weeks"] - theAvgs["wt.E16"]))

## take-home challenges

## another parametrization to try
summary(lm(gExp ~ gType/devStage, miniDat))

xyplot(fitted(twoFactFit) ~ fitted(grpFit))
xyplot(resid(twoFactFit) ~ resid(grpFit))

## using a contrast matrix
(contMat <- rbind(c(1, 0, 0, 0),
                  c(1, 1, 0, 0),
                  c(1, 0, 1, 0),
                  c(1, 0, 0, 1)))

cbind(matrixResult = as.vector(contMat %*% coef(grpFit)),
      cellMeansFit = coef(lm(gExp ~ 0 + grp, miniDat)))

##########################################################
## homegrown mlm code
##########################################################

## code Rick developed for working with mlm objects, ie fitted multivariate
## regression
source("80_anova-mlm.r")

keepMe <- with(prDes,
               devStage %in%
               levels(devStage)[c(1, nlevels(devStage))])
prDesSimple <- droplevels(subset(prDes, keepMe))
str(prDesSimple)
prDatSimple <- subset(prDat, select = prDesSimple$sidChar)

## responses must be a matrix, one *row* per response
prMatSimple <- t(as.matrix(prDatSimple))
gType <- prDesSimple$gType              # lesser of two evils
devStage <- prDesSimple$devStage        # lesser of two evils
rFit <- lm(prMatSimple ~ gType * devStage)

## summary(rFit) hangs unless you've loaded our homegrown code
## source()ed above!
rfSumm <- summary(rFit)
rfSumm                                  # by default, info on first 2
                                        # fitted models

## tour through different classes of genes w.r.t. which terms appear
## to be nonzero

##########################################################
## make it easy to find various types of 'hits'
##########################################################

## get p-values for the the main effects, interaction, overall F test
pVals <-
    cbind("NrlKO" = rfSumm$Coef[ , "gTypeNrlKO", "Pr(>|t|)"],
          "4_weeks" = rfSumm$Coef[ , "devStage4_weeks", "Pr(>|t|)"],
          "interaction" = rfSumm$Coef[ , "gTypeNrlKO:devStage4_weeks",
          "Pr(>|t|)"],
          "overallF" = rfSumm$FullModelStats[ , "Pval"])
splom(~ pVals, panel = panel.smoothScatter)
pValsTall <- data.frame(pVals = as.vector(pVals),
                        what = rep(colnames(pVals),
                        each = nrow(pVals)))
densityplot(~ pVals, pValsTall,
            group = what, auto.key = TRUE,
            plot.points = FALSE)

## NOTE written January 2013 for my future self: I got better at this
## task of finding interesting genes, graphing them, and printing some
## model results as I went along. So the next analysis -- diff exp
## anal by gType * devStage -- has a superior section on this. I think
## using p-values (vs. pre-threshholded ones) is a better way to
## go. Pulling out data for a few probeets and graphing that is a
## repetitive task where a function should probably be defined
## externally and reused in various files like this one and other
## differential expression analyses.

pValCutoff <- 0.15
pHits <- apply(pVals, 2, function(yo) yo <= pValCutoff)
table(NrlKO = pHits[ , "NrlKO"], four_weeks = pHits[ , "4_weeks"],
      interaction = pHits[ , "interaction"],
      overallF = pHits[ , "overallF"])

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
    miniDat <- droplevels(subset(miniDat,
                                 devStage %in%
                                 levels(devStage)[c(1, nlevels(devStage))]))
    miniDat
}

## handy for stripplotting the data for a small set of probesets with
## our 2x2 factorial mentality in mind
graphIt <- function(jDat) {
    print(stripplot(gExp ~ devStage | gene, jDat,
                    type = c('p', 'a'), grid = TRUE,
                    group = gType, auto.key = TRUE, jitter.data = TRUE))
}


## boring
getMe <- which(rowSums(pHits) == 0)
length(getMe)
set.seed(111)
(boring <- sample(getMe, size = 3))
graphIt(extractIt(boring))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-boring.pdf",
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMatSimple) == names(boring)[2]))

## devStage YES NrlKO NO interaction NO
getMe <- which(pHits[ , "overallF"] &
                !pHits[ , "interaction"] &
                pHits[ , "4_weeks"] &
                !pHits[ , "NrlKO"])
length(getMe)                           # 3666
set.seed(111)
(devStageOnly <- sample(getMe, size = 3))
graphIt(extractIt(devStageOnly))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-4_weeks.pdf",
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMatSimple) == names(devStageOnly)[2]))

## devStage NO NrlKO YES interaction NO
getMe <- which(pHits[ , "overallF"] &
               !pHits[ , "interaction"] &
               !pHits[ , "4_weeks"] &
               pHits[ , "NrlKO"])
length(getMe)                           # 1043
set.seed(111)
(NrlKoOnly <- sample(getMe, size = 3))
graphIt(extractIt(NrlKoOnly))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-NrlKO.pdf",
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMatSimple) == names(NrlKoOnly)[2]))

## devStage YES NrlKO YES interaction NO
getMe <- which(pHits[ , "overallF"] &
               !pHits[ , "interaction"] &
               pHits[ , "4_weeks"] &
               pHits[ , "NrlKO"])
length(getMe)                           # 1307
set.seed(222)
(mainEffOnly <- sample(getMe, size = 3))
graphIt(extractIt(mainEffOnly))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-mainYESinterNO.pdf",
          width = 7, height = 7)
print(rfSumm, show = which(colnames(prMatSimple) == names(mainEffOnly)[2]))

## devStage YES NrlKO YES interaction YES
## the randomly selected ones were not as good as some I had found
## earlier with more stringent methods
(exciting <- which(colnames(prMatSimple) %in% c("1434709_at", "1458220_at", "1455695_at")))
graphIt(extractIt(exciting))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-mainYESinterYES.pdf",
          width = 7, height = 7)
## little different here, due to hard-wiring
print(rfSumm, show = exciting[2])






