library(qvalue)
library(lattice)

nStats <- 10000

set.seed(111503)

pi0 <- 0.85                             # proportion of truly null

jDat <- data.frame(yNull = rnorm(nStats),
                   yExciting = rnorm(nStats, mean = 0.75),
                   stateOfNature = rbinom(n = nStats, size = 1, prob = pi0))
jDat$y <- with(jDat, ifelse(stateOfNature, yNull, yExciting))

head(jDat)

densityplot(~ yNull + yExciting + y, jDat,
            allow.multiple = TRUE, auto.key = TRUE,
            type = c("p", "r", "g"))

# dev.print(pdf,
#           jPaste(whereAmI, "figs/densityplot-variousPvalues.pdf"),
#           width = 8, height = 8)

foo <- lapply(jDat[c("yNull", "yExciting", "y")],
              function(z) {
                return(cbind(pnorm(z, lower.tail = FALSE),
                             pnorm(abs(z), lower.tail = FALSE) * 2))})
foo <- do.call(cbind, foo)
colnames(foo) <- paste0(c("rtProbs", "tsProbs"),
                        rep(c("Null", "Exc", "Obs"), each = 2))

jDat <- data.frame(jDat, foo)

histogram(~ rtProbsNull + tsProbsNull, jDat,
          endpoints = c(0, 1), allow.multiple = TRUE, outer = TRUE,
          xlab = "p-values for null data")

histogram(~ rtProbsExc + tsProbsExc, jDat,
          endpoints = c(0, 1), allow.multiple = TRUE, outer = TRUE,
          xlab = "p-values for exciting data!")

histogram(~ rtProbsObs + tsProbsObs, jDat,
          endpoints = c(0, 1), allow.multiple = TRUE, outer = TRUE,
          xlab =
          'p-values for realistic data (mostly null + bit of exciting)')

histogram(~ rtProbsNull + tsProbsNull +
          rtProbsExc + tsProbsExc +
          rtProbsObs + tsProbsObs, jDat,
          endpoints = c(0, 1), allow.multiple = TRUE, outer = TRUE,
          xlab = "p-values", layout = c(2, 3),
          scale = "free")

# dev.print(pdf,
#           jPaste(whereAmI, "figs/histogram-variousPvalues.pdf"),
#           width = 8, height = 6.5)


ks.test(x = jDat$rtProbsNull, y = "punif")
ks.test(x = jDat$rtProbsExc, y = "punif")
ks.test(x = jDat$rtProbsObs, y = "punif")

ks.test(x = jDat$tsProbsNull, y = "punif")
ks.test(x = jDat$tsProbsExc, y = "punif")
ks.test(x = jDat$tsProbsObs, y = "punif")

foo <- qvalue(jDat$rtProbsNull)
foo$pi0
plot(foo)

# dev.print(pdf,
#           jPaste(whereAmI, "figs/qvaluePlot-rtProbsNull.pdf"),
#           width = 8, height = 8)

foo <- qvalue(jDat$rtProbsExc)
foo$pi0
plot(foo)

# dev.print(pdf,
#           jPaste(whereAmI, "figs/qvaluePlot-rtProbsExc.pdf"),
#           width = 8, height = 8)

foo <- qvalue(jDat$rtProbsObs)
foo$pi0
plot(foo)

# dev.print(pdf,
#           jPaste(whereAmI, "figs/qvaluePlot-rtProbsObs.pdf"),
#           width = 8, height = 8)

addmargins(table(foo$qvalues < 0.40, jDat$stateOfNature), 2)
17/30

max(jDat$rtProbsObs[foo$qvalues < 0.40])

xyplot(foo$qvalues ~ rtProbsObs | stateOfNature, jDat,
       groups = yo, auto.key = TRUE)
##Error in eval(expr, envir, enclos) : object 'yo' not found

# dev.print(pdf,
#           jPaste(whereAmI, "figs/qvalue-vs-pvalue--rtProbsObs.pdf"),
#           width = 8, height = 8)

densityplot(~ foo$qvalues | stateOfNature, jDat)

foo <- qvalue(jDat$tsProbsNull)
foo$pi0
plot(foo)

foo <- qvalue(jDat$tsProbsExc)
foo$pi0
plot(foo)

foo <- qvalue(jDat$tsProbsObs)
foo$pi0
plot(foo)

addmargins(table(foo$qvalues < 0.70, jDat$stateOfNature), 2)
17/40

max(jDat$rtProbsObs[foo$qvalues < 0.40])

xyplot(foo$qvalues ~ tsProbsObs | stateOfNature, jDat,
       groups = yo, auto.key = TRUE)
## Error in eval(expr, envir, enclos) : object 'yo' not found

densityplot(~ foo$qvalues | stateOfNature, jDat)



## draw the null dist'n
tstStat <- 1.3
nullDistn <- data.frame(testStats = seq(-3.5, 3.5, length = 300))
nullDistn$densNorm <- dnorm(nullDistn$testStats)

with(nullDistn,
     plot(x = nullDistn$testStats,
          y = nullDistn$densNorm, type = "l",
          xlab = "test statistics", ylab = ""))

abline(v = tstStat, lty = "dashed")

foo <- subset(nullDistn, testStats > tstStat)
polygon(x = c(foo$testStats, rev(foo$testStats)),
        y = c(rep(0, nrow(foo)), rev(foo$densNorm)),
        col = "blue")

# dev.print(pdf,
#           jPaste(whereAmI, "figs/densityplot-rightTailPvalue.pdf"),
#           width = 8, height = 8)

polygon(x = -1 * c(foo$testStats, rev(foo$testStats)),
        y = c(rep(0, nrow(foo)), rev(foo$densNorm)),
        col = "blue")

# dev.print(pdf,
#           jPaste(whereAmI, "figs/densityplot-twoSidedPvalue.pdf"),
#           width = 8, height = 8)


