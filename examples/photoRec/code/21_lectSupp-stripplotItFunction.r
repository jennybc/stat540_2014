stripplotIt <- function(myData, ...) {
    stripplot(gExp ~ devStage | gene, myData,
              group = gType, jitter.data = TRUE,
              auto.key = TRUE, type = c('p', 'a'), grid = TRUE, ...)
}
