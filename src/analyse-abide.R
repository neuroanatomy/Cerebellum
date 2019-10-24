# Fit linear model to evaluate the effect of group and covariables on cerebellum volume

library(foreign)

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


script.dir <- get.script.dir()
base.dir <- system(paste("cd", script.dir, "&& git rev-parse --show-toplevel"), intern=T)
abide.dir <- file.path(base.dir, "data", "abide")


# table import
data.jmp <- read.xport(file.path(abide.dir, "cerebellum.stx"))

# quality check filter
df <- subset(data.jmp, CBANALYS == "Include")
df$DX_GROUP <- relevel(df$DX_GROUP, ref="Control")
df$SITE_ID2 <- relevel(df$SITE_ID2, ref="YALE")
df$SEX <- factor(df$SEX, labels=c("Male", "Female"))

# linear regressions to test the effect of group, site, age, IQ, sex and BV
fit.cb <- lm(CB~DX_GROUP+SITE_ID2+AGE_AT_S+FIQ2+SEX+BV, data=df)
summary(fit.cb)
confint(fit.cb)
fit.cbw <- lm(CBWM~DX_GROUP+SITE_ID2+AGE_AT_S+FIQ2+SEX+BV, data=df)
summary(fit.cbw)
confint(fit.cbw)
fit.cbg <- lm(CBGM~DX_GROUP+SITE_ID2+AGE_AT_S+FIQ2+SEX+BV, data=df)
summary(fit.cbg)
confint(fit.cbg)

# linear regressions to test the effect of group combined with site, age, IQ, sex and BV
fit.cb.int <- lm(CB~DX_GROUP*(SITE_ID2+AGE_AT_S+FIQ2+SEX+BV), data=df)
summary(fit.cb.int)
confint(fit.cb.int)
fit.cbw <- lm(CBWM~DX_GROUP*(SITE_ID2+AGE_AT_S+FIQ2+SEX+BV), data=df)
summary(fit.cbw.int)
confint(fit.cbw.int)
fit.cbg <- lm(CBGM~DX_GROUP*(SITE_ID2+AGE_AT_S+FIQ2+SEX+BV), data=df)
summary(fit.cbg.int)
confint(fit.cbg.int)
