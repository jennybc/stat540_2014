library(qvalue)
library(ggplot2)
library(plyr)

nStats <- 10000

set.seed(111503)

pi0 <- 0.85                             # proportion of truly null

jDat <- data.frame(tNull = rnorm(nStats),
                   tExciting = rnorm(nStats, mean = 0.75),
                   isNull = rbinom(n = nStats, size = 1, prob = pi0))
jDat$tObs <- with(jDat, ifelse(isNull, tNull, tExciting))
table(jDat$isNull)
#    0    1 
# 1521 8479 

head(jDat)

jDatTall <- stack(jDat, select = c(tNull, tObs, tExciting))
names(jDatTall) <- c("testStat", "source")
jDatTall$source <-
  factor(jDatTall$source, levels = c("tNull", "tObs", "tExciting"))
str(jDatTall)

p <- ggplot(jDatTall, aes(x = testStat))

p + geom_density(aes(col = source)) +
  geom_point(aes(x = testStat,
                 y = rep(c(0.00, -0.035, -0.07), each = nStats),
                 col = source), alpha = I(1/60))

dev.print(pdf,
          "../figs/densityplot-various-test-stats.pdf",
          width = 8, height = 8)

jDatTall$oneSided <- pnorm(jDatTall$testStat, lower.tail = FALSE)
jDatTall$twoSided <- pnorm(abs(jDatTall$testStat), lower.tail = FALSE) * 2

pDat <- data.frame(source = factor(jDatTall$source,
                                   levels = levels(jDatTall$source)),
                   isNull = jDat$isNull,
                   stack(jDatTall, select = c(oneSided, twoSided)))
str(pDat)
names(pDat) <- c("source", "isNull", "pVal", "test")
str(pDat)

p <- ggplot(pDat, aes(x = pVal)) + xlab("p-value")

p + geom_histogram(binwidth = 1/20, fill = I("grey40")) +
  facet_grid(test ~ source) + aes(y = ..density..)

dev.print(pdf,
          "../figs/histogram-variousPvalues.pdf",
          width = 8, height = 6.5)

ksRes <-
  aggregate(pVal ~ source * test, pDat, function(x) ks.test(x, y = "punif")$p)

with(ksRes,
     data.frame(source, test, ks.pval = round(pVal, 8)))
#      source     test    ks.pval
# 1     tNull oneSided 0.80181209
# 2      tObs oneSided 0.00000000
# 3 tExciting oneSided 0.00000000
# 4     tNull twoSided 0.83212279
# 5      tObs twoSided 0.00283168
# 6 tExciting twoSided 0.00000000

## get q-values
qVals <- dlply(pDat, ~ test * source, function(x) qvalue(x$pVal))
str(qVals, max.level = 0)
names(qVals)

## extract pi0
cbind(round(sapply(qVals, "[[", "pi0"), 2))
# oneSided.tNull     0.96
# oneSided.tObs      0.85
# oneSided.tExciting 0.20
# twoSided.tNull     1.00
# twoSided.tObs      0.98
# twoSided.tExciting 0.79

pDat$qVal <- unlist(lapply(qVals, "[[", "qvalues"))

p <- ggplot(pDat, aes(x = pVal, y = qVal))

p + geom_point() + facet_grid(test ~ source) + coord_equal()

dev.print(pdf,
          "../figs/scatterplot-qValues-vs-pValues.pdf",
          width = 8, height = 6.5)

## if we threshholded based on q-values from the "observed" test statistics, is
## the empirical FDR less than or equal to the threshhold?
realityCheck <-
  ddply(droplevels(subset(pDat, source == "tObs")),
        ~ test, function(x) {
          x <- arrange(x, qVal)
          x$nD <- seq_len(nrow(x))
          x$nFD <- cumsum(x$isNull)
          x$fdr <- with(x, nFD / nD)
          x
        })
str(realityCheck)

p <- ggplot(realityCheck, aes(x = qVal, y = fdr)) +
  xlab("q-value cutoff = putative FDR") + ylab("achieved FDR")

p + facet_wrap(~ test) + coord_equal(xlim = 0:1, ylim = 0:1) +
  geom_abline(intercept = 0, slope = 1, col = "orange", lwd = 1.5) +
  geom_point()

dev.print(pdf,
          "../figs/achieved-fdr-vs-putative-fdr.pdf",
           width = 8, height = 8)

## old code, using base graphics, to make the pretty pictures for depicting the
## p-values; not worth redoing in ggplot2; keeping for posterity, though

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


