## Note: certain tasks done here, such as pulling out data for a few
## probests and graphing or selecting probesets via logical operations
## on p-values, are done better in differential expression analyses
## that came later, such as the full-blown gType * devStage
## analysis. Ideally, some of that would propagate back to here via
## improved code or functions defined externally and used both places.

library(lattice)
library(limma)

## code Rick developed for working with mlm objects, ie fitted multivariate
## regression
source("80_anova-mlm.r")

## load photoRec data and design from the interwebs to delay filepath pain
## use the dput/dget workflow since works better with URLs 
jURL <- "http://www.ugrad.stat.ubc.ca/~stat540/examples/photoRec/data/GSE4051_design_DPUT.txt"
str(prDes <- dget(jURL))
jURL <- "http://www.ugrad.stat.ubc.ca/~stat540/examples/photoRec/data/GSE4051_data.tsv"
str(prDat <- read.table(jURL), max.level = 0)

with(prDes, table(devStage, gType))
head(subset(prDat, select = 1:5))

##########################################################
## limma
##########################################################

(dmDevStage <- model.matrix(~ devStage, prDes))
colnames(dmDevStage) <- levels(prDes$devStage)

fitDevStage <- lmFit(prDat, dmDevStage)
fitDevStage

## matrix, one row per probeset, five columns
## intercept = E16, then 4 columns for the effect of subsequent time
## points / developmental stages
head(fitDevStage$coef)
str(fitDevStage$coef)

## square root of diagonal elements of (t(X) %*% X)^{-1}
## same for each probeset, of course, in our case
head(fitDevStage$stdev.unscaled)
sqrt(diag(solve(t(dmDevStage) %*% dmDevStage))) # yep

## estimate residual std dev'n
head(fitDevStage$sigma)

## residual degress of freedom
## n - #params = 39 - 5 = 34
## same for each probeset, of course, in our case
head(fitDevStage$df.residual)

## top of page 62 in 2013 Jan 04 limma User's Guide
## "The ordinary t-statistics can be recovered by ..."
head(fitDevStage$coef/fitDevStage$stdev.unscaled/fitDevStage$sigma)
## one row per probeset, one col for wt or intercept, one col for the
## knockout effect


##########################################################
## one-way ANOVA "by hand" for one probeset
##########################################################
## working with first probeset
jDat <- data.frame(devStage = prDes$devStage,
                   gExp = unlist(prDat[1, ]))
lmRes <- lm(gExp ~ devStage, jDat)
summary(lmRes)

(foo <- with(jDat, tapply(gExp, devStage, mean)))
foo - foo[1]
coef(lmRes)
fitDevStage$coef[1, ]
## agree agree agree agree!!!!

summary(lmRes)
fitDevStage$sigma[1] * fitDevStage$stdev.unscaled[1, ] # std errors
round((fitDevStage$coef/fitDevStage$stdev.unscaled/
       fitDevStage$sigma)[1, ], 5) # t stats
## agree agree agree agree!!!!

##########################################################
## homegrown mlm code
##########################################################

## responses must be a matrix, one *row* per response
prMat <- t(as.matrix(prDat))
devStage <- prDes$devStage        # lesser of two evils
rFit <- lm(prMat ~ devStage)

## summary(rFit) hangs unless you've loaded our homegrown code
## source()ed at top!
rfSumm <- summary(rFit)
rfSumm                                  # by default, info on first 2
                                        # fitted models
print(rfSumm, show = c(2, 4555, 29403)) # show gives more flexibility

str(rfSumm, max.level = 1)
str(rfSumm$Coef)

## one row per probeset
## one column per parameter
## for each row * param, 4 pieces of stat inf:
## "Estimate" "Std. Error" "t value" "Pr(>|t|)"

## YES getting same wild type or intercept estimates
head(cbind(coef(rFit)["(Intercept)", ],
           coef(fitDevStage)[ , "E16"]))

## YES getting same devStage effect estimates
head(cbind(coef(rFit)["devStageP2", ],
           coef(fitDevStage)[ , "P2"]))

## YES getting same estimated residual variance
head(cbind(rFit = deviance(rFit)/df.residual(rFit),
           limma = fitDevStage$sigma ^ 2))

## estimated std err of estimated intercept = wt value
head(cbind(rFit = rfSumm$Coef[ , "(Intercept)", "Std. Error"],
           limma = fitDevStage$sigma * fitDevStage$stdev.unscaled[ , "E16"]))

## estimated std err of estimated P6 effect
head(cbind(rFit = rfSumm$Coef[ , "devStageP6", "Std. Error"],
           limma = fitDevStage$sigma * fitDevStage$stdev.unscaled[ , "P6"]))

## t statistic for 4_weeks effect
head(cbind(rFit = rfSumm$Coef[ , "devStage4_weeks", "t value"],
           limma = coef(fitDevStage)[ , "4_weeks"] /
           (fitDevStage$sigma * fitDevStage$stdev.unscaled[ , "4_weeks"])))

## matrix with one row per probeset
## columns/variables are est, se, t stat, p val
str(rfSumm$Coefficients[ , "devStageP10", ])
## num [1:29949, 1:4] 0.000161 0.143103 -0.045129 -0.041076 0.015492 ...
## - attr(*, "dimnames")=List of 2
##  ..$ : chr [1:29949] "1415670_at" "1415671_at" "1415672_at" "1415673_at" ...
##  ..$ : chr [1:4] "Estimate" "Std. Error" "t value" "Pr(>|t|)"

## matrix with one row per probeset
## columns are "Res. SE"   "Fstat"     "Pval"      "R-Squared" "Adj R-Sqr"
## for our model Fstat is t^2 and Pval is the t-test p-value.
str(rfSumm$FullModelStats)

## YES getting same estimated residual variance
head(cbind(rFit1 = deviance(rFit)/df.residual(rFit),
           rFit2 = rfSumm$Full[ , "Res. SE"] ^ 2,
           limma = fitDevStage$sigma ^ 2))

##########################################################
## limma continued ... Empirical Bayes
##########################################################

## not populated yet, i.e. not before eBayes()
head(fitDevStage$F)
fitDevStage$s2.prior
fitDevStage$df.prior
head(fitDevStage$var.prior)
head(fitDevStage$t)
## NULL

ebFitDevStage <- eBayes(fitDevStage)

## Yes these are the same
head(cbind(coef(ebFitDevStage),
           coef(fitDevStage)))

## this is the mean of the inverse Chisquare prior for the
## gene-specific variances
ebFitDevStage$s2.prior
## [1] 0.06624542

## this is the degrees of freedom associated with the inverse
## ChiSquare prior for gene-specific variances
ebFitDevStage$df.prior
## [1] 2.881327

head(cbind(rFit1 = deviance(rFit)/df.residual(rFit),
           rFit2 = rfSumm$Full[ , "Res. SE"] ^ 2,
           limma = fitDevStage$sigma ^ 2,
           ebLimma = ebFitDevStage$s2.post))

## t statistics for P2
head(cbind(rFit1 = rfSumm$Coef[ , "devStageP2", "t value"],
           limma = coef(fitDevStage)[ , "P2"] /
           (fitDevStage$sigma * fitDevStage$stdev.unscaled[ , "P2"]),
           ebLimma = ebFitDevStage$t[ , "P2"]))

## pvalues for P10
head(cbind(rFit = rfSumm$Coef[ , "devStageP10", "Pr(>|t|)"],
           ebLimma = ebFitDevStage$p.value[ , "P10"]))

## look at some hits
(foo <- topTable(ebFitDevStage))
(getMe <- rownames(foo)[1:6])
cbind(getMe)
## [1,] "1438940_x_at"
## [2,] "1450084_s_at"
## [3,] "1455897_x_at"
## [4,] "1456736_x_at" # mitochondrial fission factor Mff
## [5,] "1426384_a_at"
## [6,] "1460547_a_at"
jDat <- prDat[getMe, ]
jDat <- t(jDat)
jDat <- data.frame(gExp = as.vector(jDat),
                   probeset = rep(colnames(jDat),
                   each = nrow(jDat)))
kDat <- data.frame(prDes, jDat)
str(kDat)

stripplot(gExp ~ devStage | probeset, kDat,
          type = c('p', 'a'), grid = TRUE)

stripplot(gExp ~ devStage | probeset, kDat,
          groups = devStage, auto.key = TRUE)

## look at some non-hits
n <- nrow(prDat)
foo <- topTable(ebFitDevStage, n = Inf)
foo <- topTable(ebFitDevStage, n = Inf)[(n - 5):n, ]
(getMe <- rownames(foo)[1:6])
cbind(getMe)
##      getMe
## [1,] "1458418_at"
## [2,] "1426288_at"
## [3,] "1425171_at"
## [4,] "1436646_at"
## [5,] "1451617_at"
## [6,] "1450946_at"
jDat <- prDat[getMe, ]
jDat <- t(jDat)
jDat <- data.frame(gExp = as.vector(jDat),
                   probeset = rep(colnames(jDat),
                   each = nrow(jDat)))
kDat <- data.frame(prDes, jDat)
str(kDat)

## hello gType * devStage interaction effects!
stripplot(gExp ~ devStage | probeset, kDat,
          type = c('p', 'a'), grid = TRUE,
          group = gType)
## not what I really wanted, though

## another try to find some boring non-hits
## matrix with one row per probeset
## columns are "Res. SE"   "Fstat"     "Pval"      "R-Squared" "Adj R-Sqr"
## for our model Fstat is t^2 and Pval is the t-test p-value.
str(rfSumm$FullModelStats)
head(rfSumm$FullModelStats)
head((rfSumm$FullModelStats)[ , "Fstat"])
(getMe <- which(rank( (rfSumm$FullModelStats)[ , "Fstat"]) <= 10))
jDat <- prDat[getMe, ]
jDat <- t(jDat)
jDat <- data.frame(gExp = as.vector(jDat),
                   probeset = rep(colnames(jDat),
                   each = nrow(jDat)))
kDat <- data.frame(prDes, jDat)
str(kDat)

## hello gType * devStage interaction effects!
stripplot(gExp ~ devStage | probeset, kDat,
                         type = c('p', 'a'), grid = TRUE,
          group = gType)

## 1447281_at is super boring, row 21057
## 1443184_at also nice, row 18898

deDevStage <-
    data.frame(F = (rfSumm$FullModelStats)[ , "Fstat"],
               pVal = (rfSumm$FullModelStats)[ , "Pval"],
               sigma = (rfSumm$FullModelStats)[ , "Res. SE"],
               rank = rank((rfSumm$FullModelStats)[ , "Pval"]))
peek(deDevStage)

str(deDevStage)

write.table(deDevStage,
            file = "../results/deDevStage.txt",
            quote = FALSE)


