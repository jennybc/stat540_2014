## makes stripplot for a handful of probesets ... used in all my lecture support
## scripts that show what different kinds of DEA "hits" look like

stripplotIt <- function(myData, ...) {
    stripplot(gExp ~ devStage | gene, myData,
              group = gType, jitter.data = TRUE,
              auto.key = TRUE, type = c('p', 'a'), grid = TRUE, ...)
}
