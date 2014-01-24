# see 20_dea-by-devStage for a more chatty, educational tour through
# this analysis; this is stripped down to bare minimum

## load design data.frame
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes)

## load gene expression data.frame
prDat <- read.table("../data/GSE4051_data.tsv")
str(prDat, max.level = 0)
## 'data.frame':  29949 obs. of  39 variables:

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

## code Rick developed for working with mlm objects, ie fitted multivariate
## regression
source("80_anova-mlm.r")

## responses must be a matrix, one *row* per response
prMat <- t(as.matrix(prDat))
devStage <- prDes$devStage        # lesser of two evils
rFit <- lm(prMat ~ devStage)

## summary(rFit) hangs unless you've loaded our homegrown code
## source()ed above!
rfSumm <- summary(rFit)
str(rfSumm, max.level = 1)
str(rfSumm$Coef)
## one row per probeset
## one column per parameter
## for each row * param, 4 pieces of stat inf:
## "Estimate" "Std. Error" "t value" "Pr(>|t|)"

deDevStage <-
  data.frame(F = (rfSumm$FullModelStats)[ , "Fstat"],
             pVal = (rfSumm$FullModelStats)[ , "Pval"],
             sigma = (rfSumm$FullModelStats)[ , "Res. SE"],
             rank = rank((rfSumm$FullModelStats)[ , "Pval"]))

write.table(deDevStage,
            file = "../results/deDevStage.txt",
            quote = FALSE)


##########################################################
## limma: use it to fit linear model
##        gene expression ~ devStage
##        for each probeset
##########################################################

# library(limma)
# 
# (dmDevStage <- model.matrix(~ devStage, prDes))
# colnames(dmDevStage) <- levels(prDes$devStage)
# 
# fitDevStage <- lmFit(prDat, dmDevStage)
# fitDevStage
# 
# ebFitDevStage <- eBayes(fitDevStage)


