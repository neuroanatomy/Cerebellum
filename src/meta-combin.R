library(meta)
library(metafor)

get.script.dir <- function(){
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  sourceDir <- getSrcDirectory(function(dummy) {dummy})
  if (length(script.name)) { # called from command
    (dirname(script.name))
  } else if (nchar(sourceDir)) { # called with source
    sourceDir
  } else if (rstudioapi::isAvailable()) { # called from RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else getwd()
}

suffix <- "combin"
csv.file <- sprintf("means-%s.txt", suffix)

script.dir <-get.script.dir()
base.dir <- system(paste("cd", script.dir, "&& git rev-parse --show-toplevel"), intern=T)
meta.dir <- file.path(base.dir, "data", "meta-analysis")

dataCbl <- read.table(file.path(meta.dir, "means-Cbl.txt"), h=T, sep="\t")
dataAbide <- read.table(file.path(meta.dir, "means-abide.txt"), h=T, sep="\t")

dataComb <- rbind(dataCbl, dataAbide)
write.table(dataComb, paste(script.dir, csv.file, sep="/"), row.names=F, quote=F, sep="\t")

source(file.path(script.dir, "meta-analysis.R"))

funnel(meta)
funnel(meta, pch=c(rep(21, nrow(dataCbl)), rep(24, nrow(dataAbide))))
legend("topright", inset=.02, legend = c("Literature", "ABIDE"),
       pch = c(21, 24), pt.bg = "grey")
