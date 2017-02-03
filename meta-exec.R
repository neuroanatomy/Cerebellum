library(meta)
library(metafor)

getsd <- function() {
  path <- try(sys.frame(1)$ofile, silent=T)
  if (is.null(path)) {
    path <- paste(getSrcDirectory(function(dummy) {dummy}), "dummy", sep="/")
  } else if (is.null(path)) {
    # Rscript
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    path <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  }
  dirname(path)
}

csv.file <- sprintf("means-%s.txt", suffix)
pcurve.file <- sprintf("pcurve-%s.txt", suffix)

script.dir <- getsd()
alpha <- 0.05

data <- read.table(paste(script.dir, csv.file, sep="/"), h=T, sep="\t")
data[data == 0] <- NA

# Meta-analysis on standard mean differences

meta_SMD <- metacont(data$n.asd, data$mean.asd, data$sd.asd, data$n.ctrl, data$mean.ctrl, data$sd.ctrl,
		 data=data, sm="SMD", studlab=label, comb.fixed=F, method.tau="REML")

# Independent two groups Student's test

pooledsd <- sqrt(((data$n.asd-1)*data$sd.asd^2+(data$n.ctrl-1)*data$sd.ctrl^2)/(data$n.asd+data$n.ctrl-2))
df <- data$n.asd+data$n.ctrl-2
tvalue <- (data$mean.asd-data$mean.ctrl)/(pooledsd * sqrt(1/data$n.asd+1/data$n.ctrl))

ttext <- paste("t(", df, ")=", tvalue, sep="")
writeLines(ttext, paste(script.dir, pcurve.file, sep="/"))

p_SMD <- 2*pt(-abs(tvalue), df=df)

# Meta-regression

reg_SMD <- metareg(meta_SMD, ~age.asd+iq.asd)
reg_SMD$slab <- data$label

# Power for effect size d estimated by mixed effect model

ct <- qt(1-alpha/2, df)

TE.reg_SMD <- colSums(rbind(1, data$age.asd, data$iq.asd) * c(reg_SMD$b))
ncp.reg_SMD <- TE.reg_SMD/sqrt(1/data$n.asd+1/data$n.ctrl)
power.reg_SMD <- 1 - pt(ct, df, ncp.reg_SMD) + pt(-ct, df, ncp.reg_SMD)

# Publication bias

bias_SMD <- metabias(meta_SMD, k.min=5)

# p-hacking

source(paste(script.dir, "p-curve.R", sep="/"))
pcurve_SMD <- pcurve(ttext)

# Meta-analysis on ratio of variances

# Remove studies with strictly equal variances
nequal = sum(data$sd.asd == data$sd.ctrl)
if (nequal)
	print(sprintf(paste("warning: found %d studies with exactly equal standard deviation. ",
				"Removing them from meta-analysis on ratio of variances",
				" because it seems to be artificial.", sep=""), nequal))
data_F <- data[data$sd.asd != data$sd.ctrl, ]

# Independent groups F test on ratio of variance
fvalue <- data_F$sd.asd^2 / data_F$sd.ctrl^2
df1 <- data_F$n.asd - 1
df2 <- data_F$n.ctrl - 1

ftext <- paste("F(", df1, ",", df2, ")=", fvalue, sep="")

p_F <- 1 - abs(pf(fvalue, df1, df2) - 0.5) * 2
ci_F <- cbind(fvalue/qf(1-alpha/2, df1, df2), fvalue/qf(alpha/2, df1, df2))

TE_F <- log(fvalue) - digamma(df1/2) + digamma(df2/2) + log(df1/df2)
TEse_F <- sqrt(trigamma(df1/2) + trigamma(df2/2))
meta_F <- metagen(TE_F, TEse_F, studlab=data_F$label, data=data_F, comb.fixed=F)
reg_F <- metareg(meta_F, ~age.asd+iq.asd)
reg_F$slab <- data_F$label


# Power for effect size d estimated by mixed effect model

# critical F
cf <- cbind(qf(alpha/2, df1, df2), qf(1-alpha/2, df1, df2))

TE.reg_F <- colSums(rbind(1, data_F$age.asd, data_F$iq.asd) * c(reg_F$b))
Fr.reg <- exp(TE.reg_F)
power.reg_F <- pf(cf[,1]/Fr.reg, df1, df2, lower.tail=T) + pf(cf[,2]/Fr.reg, df1, df2, lower.tail=F)

# Publication bias
bias_F <- metabias(meta_F, k.min=5)

# p-hacking
pcurve_F <- pcurve(ftext)
