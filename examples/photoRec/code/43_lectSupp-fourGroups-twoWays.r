library(lattice)

## load design data.frame
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes)

## load gene expression data.frame
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)
## 'data.frame':  29949 obs. of  39 variables:

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

##########################################################
## simplying devStage to first and last timepoints
##########################################################
keepMe <- with(prDes,
               devStage %in%
                 levels(devStage)[c(1, nlevels(devStage))])
prDes <- droplevels(subset(prDes, keepMe))
str(prDes)
prDat <- subset(prDat, select = prDes$sidChar)
str(prDat)

##########################################################
## show an example
##########################################################

## I've already done differential expression analysis for gType * a
## simplified version of devStage, including only E16 and 4_weeks

## I've selected this as our example
(luckyGene <- which(rownames(prDat) == "1455695_at")) # 26861
miniDat <- data.frame(gExp = unlist(prDat[luckyGene, ]))
miniDat <- data.frame(prDes, miniDat)
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
## load DEA results computed in 05_dea-by-gType-oneDevStage.r
##########################################################
deaRes <- read.table("../results/deGtypeBy4weeks.txt", header = TRUE)
str(deaRes)

##########################################################
## marshall p-values to help us find various types of 'hits'
##########################################################

## get p-values for both the main effects, interaction, overall F test
pVals <- subset(deaRes,
                select = c("pVal.NrlKO", "pVal.4_weeks",
                           "pVal.NrlKO.4_weeks", "Fp"))
names(pVals) <- c("NrlKO", "4_weeks", "interaction", "overallF")
splom(~ pVals, panel = panel.smoothScatter)
pValsTall <- data.frame(pVals = unlist(pVals),
                        what = rep(colnames(pVals), each = nrow(pVals)))
densityplot(~ pVals, pValsTall,
            group = what, auto.key = TRUE, plot.points = FALSE, n = 200)

##########################################################
## helper functions for finding and featuring genes
##########################################################

source("98_lectSupp-prepareDataFunction.r")
source("99_lectSupp-stripplotItFunction.r")

##########################################################
## designate "hits"
##########################################################

## NOTE written January 2013 and updated January 2014, for my future self: I got
## better at this task of finding interesting genes, graphing them, and printing
## some model results as I went along. So the next analysis -- DEA by gType *
## devStage -- has a superior approach for using the p-values to find hits. I
## haven't fully switched over here because I amd [1] lazy and [2] want to keep
## the code that produces the figures actually in the lecture slides.

pValCutoff <- 0.15
pHits <- apply(pVals, 2, function(yo) yo <= pValCutoff)
table(NrlKO = pHits[ , "NrlKO"], four_weeks = pHits[ , "4_weeks"],
      interaction = pHits[ , "interaction"],
      overallF = pHits[ , "overallF"])

# NOTE to self January 2014: it appears that the "boring" genes featured in 
# lecture slides are NOT the ones I find below. I do reproduce everything else.
# So I suspect I made that slide before adding set.seed() to this code. I'm not 
# replacing the figure and result because I was also using a different lattice 
# theme at the time and then this new figure would stick out. Sigh.

## NOTE to self January 2014: in the past I re-did the DEA in this script,
## duplicating the efforts in 05_dea-by-gType-oneDevStage.r. This gave me access
## to the full model results in the summary of the mlm object. Now that I am
## just loading stored results, I don't have access to such detail. That's why
## I've commented out the lines that produced typical "summary of a lm". They
## don't work but this is a reminder of how I got the reporting found on current
## lecture slides.

## boring genes
## not a hit w/r/t either main effect OR the interaction OR overall
getMe <- which(rowSums(pHits) == 0)
length(getMe) # 13315
set.seed(111)
(boring <- sample(getMe, size = 3))
stripplotIt(prepareData(boring))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-boring.pdf",
          width = 7, height = 7)
#print(rfSumm, show = which(colnames(prMatSimple) == names(boring)[2]))

## devStage YES NrlKO NO interaction NO
getMe <- which(pHits[ , "overallF"] &
                !pHits[ , "interaction"] &
                pHits[ , "4_weeks"] &
                !pHits[ , "NrlKO"])
length(getMe)                           # 3666
set.seed(111)
(devStageOnly <- sample(getMe, size = 3))
stripplotIt(prepareData(devStageOnly))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-4_weeks.pdf",
          width = 7, height = 7)
#print(rfSumm, show = which(colnames(prMatSimple) == names(devStageOnly)[2]))

## devStage NO NrlKO YES interaction NO
getMe <- which(pHits[ , "overallF"] &
               !pHits[ , "interaction"] &
               !pHits[ , "4_weeks"] &
               pHits[ , "NrlKO"])
length(getMe)                           # 1043
set.seed(111)
(NrlKoOnly <- sample(getMe, size = 3))
stripplotIt(prepareData(NrlKoOnly))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-NrlKO.pdf",
          width = 7, height = 7)
#print(rfSumm, show = which(colnames(prMatSimple) == names(NrlKoOnly)[2]))

## devStage YES NrlKO YES interaction NO
getMe <- which(pHits[ , "overallF"] &
               !pHits[ , "interaction"] &
               pHits[ , "4_weeks"] &
               pHits[ , "NrlKO"])
length(getMe)                           # 1307
set.seed(222)
(mainEffOnly <- sample(getMe, size = 3))
stripplotIt(prepareData(mainEffOnly))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-mainYESinterNO.pdf",
          width = 7, height = 7)
#print(rfSumm, show = which(colnames(prMatSimple) == names(mainEffOnly)[2]))

## devStage YES NrlKO YES interaction YES
## the randomly selected ones were not as good as some I had found
## earlier with more stringent methods
(exciting <-
   which(rownames(prDat) %in% c("1434709_at", "1458220_at", "1455695_at")))
stripplotIt(prepareData(exciting))
dev.print(pdf, "../figs/gTypeBy4weeks/stripplot-mainYESinterYES.pdf",
          width = 7, height = 7)
## little different here, due to hard-wiring
#print(rfSumm, show = exciting[2])