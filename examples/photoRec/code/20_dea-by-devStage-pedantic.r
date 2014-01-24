##########################################################
## homegrown mlm code
##########################################################


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


