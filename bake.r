library(fileherd)

setwd("~/teaching/STAT540/2014_STAT540")

findPreview()

## be careful! do not want to delete photoRec/README.md!
findPreview(remove = TRUE, verbose = TRUE)

herd(pub.dir = "/Users/jenny/web/STAT540", backend = "pandoc")

