library(meta)
library(metafor)

getsd <- function() {
  path <- try(sys.frame(1)$ofile, silent=T)
  if (is.null(path)) {
    # Rscript
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    path <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  }
  dirname(path)
}

suffix <- "combin"
csv.file <- sprintf("means-%s.txt", suffix)

script.dir <- getsd()

dataCbl <- read.table(paste0(script.dir, "/means-Cbl.txt"), h=T, sep="\t")
dataAbide <- read.table(paste0(script.dir, "/means-abide.txt"), h=T, sep="\t")

dataComb <- rbind(dataCbl, dataAbide)
write.table(dataComb, paste(script.dir, csv.file, sep="/"), row.names=F, quote=F, sep="\t")

source(paste(script.dir, "meta-analysis.R", sep="/"))

funnel(meta)
funnel(meta, pch=c(rep(21, nrow(dataCbl)), rep(24, nrow(dataAbide))))
legend("topright", inset=.02, legend = c("Literature", "ABIDE"),
       pch = c(21, 24), pt.bg = "grey")
