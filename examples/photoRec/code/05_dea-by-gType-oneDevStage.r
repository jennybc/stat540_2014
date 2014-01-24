# see 21_dea-by-gType-oneDevStage for a more chatty, educational tour through
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

##########################################################
## simplying devStage to first and last timepoints
##########################################################
keepMe <- with(prDes,
               devStage %in%
                 levels(devStage)[c(1, nlevels(devStage))])
prDesSimple <- droplevels(subset(prDes, keepMe))
str(prDesSimple)
prDatSimple <- subset(prDat, select = prDesSimple$sidChar)
str(prDatSimple)

## code Rick developed for working with mlm objects, ie fitted multivariate
## regression
source("80_anova-mlm.r")

## responses must be a matrix, one *row* per response
prMatSimple <- t(as.matrix(prDatSimple))
gType <- prDesSimple$gType              # lesser of two evils
devStage <- prDesSimple$devStage        # lesser of two evils
rFit <- lm(prMatSimple ~ gType * devStage)

## summary(rFit) hangs unless you've loaded our homegrown code
## source()ed at top!
rfSumm <- summary(rFit)

est_coefs <- t(coef(rFit))

deGtypeBy4weeks <-
  data.frame(est_coefs, # est coefs
             rfSumm$Coef[ , , "Std. Error"], # est se's
             rfSumm$Coef[ , , "t value"], # t stats
             rfSumm$Coef[ , , "Pr(>|t|)"], # p vals
             deviance(rFit)/df.residual(rFit), # resid var
             rfSumm$FullModelStats[ , "Fstat"], # F test
             rfSumm$FullModelStats[ , "Pval"])

yo <- c(rep(c("est", "se", "tStat", "pVal"), each = 4))
yo2 <- c(rep(c("int", "NrlKO", "4_weeks", "NrlKO:4_weeks"), 4))
colnames(deGtypeBy4weeks) <-
  c(paste(yo, yo2, sep = "."), c("sigSq", "Fstat", "Fp"))

head(deGtypeBy4weeks)
str(deGtypeBy4weeks)

write.table(deGtypeBy4weeks,
            file = "../results/deGtypeBy4weeks.txt",
            quote = FALSE, row.names = FALSE)


##########################################################
## limma
##########################################################
# (dmGtypeStarDevStage <- model.matrix(~ gType * devStage, prDesSimple))
# fitGtypeStarDevStage <- lmFit(prDatSimple, dmGtypeStarDevStage)
# ebFitGtypeStarDevStage <- eBayes(fitGtypeStarDevStage)
