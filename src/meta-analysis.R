library(meta)
library(metafor)
library(reshape)
library(ggplot2)

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

combinedsd <- function(n, mean, sd, na.rm=F) {
  data <- data.frame(n, mean, sd)
  if (na.rm)
    data <- na.omit(data)
  with(data, {
    ncomb <- sum(n)
    meancomb <- weighted.mean(mean, n)
    sqrt(sum((n-1)*sd^2+n*(meancomb-mean)^2)/(ncomb-1))
  })
}

# Set suffix according to region of interest
suffix <- "Cbl"

# input files
csv.file <- sprintf("means-%s.txt", suffix)
# output files
pcurve.file <- sprintf("pcurve-%s.txt", suffix)
data.file <- sprintf("data-%s.txt", suffix)

script.dir <- get.script.dir()
base.dir <- system(paste("cd", script.dir, "&& git rev-parse --show-toplevel"), intern=T)
meta.dir <- file.path(base.dir, "data", "meta-analysis")

data <- read.table(file.path(meta.dir, csv.file), h=T, sep="\t")
data[data == 0] <- NA

total <- data[0,]
total[1,] <- NA
total$label <- "Total"
total$n.asd <- sum(data$n.asd)
total$n.male.asd <- sum(data$n.male.asd)
total$n.ctrl <- sum(data$n.ctrl)
total$n.male.ctrl <- sum(data$n.male.ctrl)
total$age.asd <- weighted.mean(data$age.asd, data$n.asd)
total$age.ctrl <- weighted.mean(data$age.ctrl, data$n.ctrl)
total$age.sd.asd <- combinedsd(data$n.asd, data$age.asd, data$age.sd.asd)
total$age.sd.ctrl <- combinedsd(data$n.ctrl, data$age.ctrl, data$age.sd.ctrl)
total$iq.asd <- weighted.mean(data$iq.asd, data$n.asd, na.rm=T)
total$iq.ctrl <- weighted.mean(data$iq.ctrl, data$n.ctrl, na.rm=T)
total$iq.sd.asd <- combinedsd(data$n.asd, data$iq.asd, data$iq.sd.asd, na.rm=T)
total$iq.sd.ctrl <- combinedsd(data$n.ctrl, data$iq.ctrl, data$iq.sd.ctrl, na.rm=T)
total$mean.asd <- weighted.mean(data$mean.asd, data$n.asd)
total$mean.ctrl <- weighted.mean(data$mean.ctrl, data$n.ctrl)
total$sd.asd <- combinedsd(data$n.asd, data$mean.asd, data$sd.asd)
total$sd.ctrl <- combinedsd(data$n.ctrl, data$mean.ctrl, data$sd.ctrl)

# Meta-analysis on standard mean differences

meta <- metacont(data$n.asd, data$mean.asd, data$sd.asd, data$n.ctrl, data$mean.ctrl, data$sd.ctrl,
		data=data, sm="SMD", studlab=label, comb.fixed=F, method.tau="REML")

meta$number <- paste(data$n.asd, data$n.ctrl, sep="/")
meta$age <- paste(format(data$age.asd, digits=2), format(data$age.ctrl, digits=2), sep="/")
meta$iq <- paste(format(data$iq.asd, digits=2), format(data$iq.ctrl, digits=2), sep="/")
print(meta)
funnel(meta)
metabias(meta)
# forest(meta, lab.e="ASD/Control", leftcols=c("studlab", "number", "age", "iq"), leftlabs=c("Study", "Number", "Age", "IQ"), lab.e.attach.to.col=c("age"))
forest(meta, leftcols=c("studlab"), leftlabs=c("Study"))

# Independent two groups Student's test

pooledsd <- sqrt(((data$n.asd-1)*data$sd.asd^2+(data$n.ctrl-1)*data$sd.ctrl^2)/(data$n.asd+data$n.ctrl-2))
df <- data$n.asd+data$n.ctrl-2
tvalue <- (data$mean.asd-data$mean.ctrl)/(pooledsd * sqrt(1/data$n.asd+1/data$n.ctrl))

ttext <- paste("t(", df, ")=", tvalue, sep="")
writeLines(ttext, file.path(meta.dir, pcurve.file))

pvalue <- 2*pt(-abs(tvalue), df=df)

# Meta-regressions

reg1 <- metareg(meta, ~age.asd)
reg2 <- metareg(meta, ~iq.asd)
reg9 <- metareg(meta, ~age.asd+iq.asd)

# bubble(reg1)
reg9$slab <- data$label
forest.rma(reg9)
regtest(reg9)

plot2Dreg <- function(reg) {
	res <- 64
	min.age <- 0
	max.age <- 40
	min.iq <- 50
	max.iq <- 130
	
	a <- list(x=seq(min.age, max.age, len=res), y=seq(min.iq, max.iq, len=res))
	a$z <- outer(a$x, a$y, function(x, y)(reg$b[1] + reg$b[2] * x + reg$b[3] * y))
	step <- 0.25
	min.z <- floor(min(a$z)/step)*step
	max.z <- ceiling(max(a$z)/step)*step
	levels <- seq(min.z, max.z, step)
	# max.z <- ceiling(max(abs(a$z))*10)/10
	# levels <- seq(-max.z, max.z, 0.25)
	# col = rainbow(length(levels)-1, start=0, end=2/6)
	col <- grey.colors(length(levels)-1, start=0, end=1)
	
	filled.contour(a, col=col, levels=levels,
				   plot.title = title(xlab = "Mean age ASD", ylab = "Mean IQ ASD"),
				   key.title = title(main = "SMD"), plot.axes = { axis(1); axis(2);
				   points(reg$X[,2], reg$X[,3], bg=col[findInterval(reg$yi, levels)],
				   pch=21, cex=2/sqrt(reg$tau2+reg$vi)) })
}

plot2Dreg2 <- function(reg) {
  res <- 64
  min.age <- 0
  max.age <- 40
  min.iq <- 50
  max.iq <- 130
  
  a <- list(x=seq(min.age, max.age, len=res), y=seq(min.iq, max.iq, len=res))
  a$z <- outer(a$x, a$y, function(x, y)(reg$b[1] + reg$b[2] * x + reg$b[3] * y))
  rownames(a$z) <- a$x
  colnames(a$z) <- a$y
  max.z <- ceiling(max(abs(a$z))*10)/10
  levels <- seq(-max.z, max.z, 0.25)
  # col = rainbow(length(levels)-1, start=0, end=2/6)
  # col <- viridis(length(levels)-1)
  col <- grey.colors(length(levels)-1, start=0, end=1)
  
  az.melt <- melt(a$z)
  names(az.melt) <- c("age", "iq", "SMD")
  
  w = ggplot(az.melt, aes(x = age, y = iq, z = SMD))
  w2 <- w + geom_tile(aes(fill = SMD)) +
    xlab("Mean age ASD") + ylab("Mean IQ ASD") +
    guides(fill = guide_colorbar(barwidth = 2, barheight = 10))
  w3 <- w2 + coord_cartesian(xlim=c(min(az.melt$age),max(az.melt$age)), ylim=c(min(az.melt$iq),max(az.melt$iq)), expand=F)
  w4 <- w3 + geom_point(data=as.data.frame(cbind(reg$X, SMD=reg$yi)), aes(x = age.asd, y = iq.asd, fill = SMD),
                        size = 4/sqrt(reg$tau2+reg$vi), pch = 21, col="black")
  w5 <- w4 + scale_fill_gradientn(colours = col, limits = c(-max.z, max.z))
  wf <- w5 + theme_light()
  print(wf)
}

plot2Dreg(reg9)

# Meta-analysis on ratio of standard deviations

Fr <- data$sd.asd^2 / data$sd.ctrl^2
df1 <- data$n.asd - 1
df2 <- data$n.ctrl - 1
p_F <- 1 - abs(pf(Fr, df1, df2) - 0.5) * 2
alpha <- .05
ci_F <- cbind(Fr/qf(1-alpha/2, df1, df2), Fr/qf(alpha/2, df1, df2))

TE <- log(Fr) - digamma(df1/2) + digamma(df2/2) + log(df1/df2)
TEse <- sqrt(trigamma(df1/2) + trigamma(df2/2))
meta_F <- metagen(TE, TEse, studlab=data$label, data=data, comb.fixed=F, method.tau="REML")
forest(meta_F)

# Meta-regression

reg_F <- metareg(meta_F, ~age.asd+iq.asd)

# Comparison of infered volume distributions from Cohen's d

if (F) {
  x <- seq(-4, 4, .01)
  plot(x, dnorm(x), col="blue", type='l', lwd=3,
       main=sprintf("Normalized distributions for Cohen's d = %.2f", meta$TE.random), xlab = "Volume", ylab = "Density")
  lines(x, dnorm(x-meta$TE.random), col="red", lwd=3)
  legend("topright", # places a legend at the appropriate place
         c("Controls","ASD"), # puts text in the legend
         lty=c(1,1), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("blue","red")) # gives the legend lines the correct color and width
}

# Power for effect size d estimated by random model

alpha <- 0.05
ct <- qt(1-alpha/2, df)
ncp.rand <- meta$TE.rand/sqrt(1/data$n.asd+1/data$n.ctrl)
power.rand <- 1 - pt(ct, df, ncp.rand) + pt(-ct, df, ncp.rand)

# Power for effect size d estimated by regression model

TE.reg <- colSums(rbind(1, data$age.asd, data$iq.asd) * c(reg9$b))
ncp.reg <- TE.reg/sqrt(1/data$n.asd+1/data$n.ctrl)
power.reg <- 1 - pt(ct, df, ncp.reg) + pt(-ct, df, ncp.reg)

# Power of the meta-analysis given an expected heterogeneity and effect size

es <- 0.3 # Expected effect size
hg <- 1   # Expected heterogeneity (".33" for small, "1" for moderate, & "3" for large)
# hg <- 1 / (1 - meta$I2) - 1

V <- 1/data$n.asd+1/data$n.ctrl
T2 <- hg / (mean(1 / V))
Vm <- 1 / (sum(1 / (V + T2)))
ncp.z = es/sqrt(Vm)
cz <- qnorm(1-alpha/2)
power.meta = pnorm(-cz, ncp.z) + pnorm(cz, ncp.z, lower.tail = F)

# Data export

dataf <- rbind(data, total)
dataf$n.female.asd <- dataf$n.asd - dataf$n.male.asd
dataf$n.female.ctrl <- dataf$n.ctrl - dataf$n.male.ctrl
dataf <- format(dataf, digits=1, nsmall=1, trim=T)
dataf[dataf == "NA"] <- "-"

table <- data.frame(dataf$label, sprintf("%s (%s)", dataf$n.asd, dataf$n.female.asd), sprintf("%s (%s)", dataf$n.ctrl, dataf$n.female.ctrl),
                     sprintf("%s ± %s", dataf$age.asd, dataf$age.sd.asd), sprintf("%s ± %s", dataf$age.ctrl, dataf$age.sd.ctrl),
                     sprintf("%s ± %s", dataf$iq.asd, dataf$iq.sd.asd), sprintf("%s ± %s", dataf$iq.ctrl, dataf$iq.sd.ctrl),
                     sprintf("%s ± %s", dataf$mean.asd, dataf$sd.asd), sprintf("%s ± %s", dataf$mean.ctrl, dataf$sd.ctrl))

names(table) <- c("Study", "N ASD (F)", "N Ctrl (F)", "Age ASD", "Age Ctrl", "IQ ASD", "IQ Ctrl", "Measure ASD", "Measure Ctrl")
write.table(table, file.path(meta.dir, data.file), row.names=F, quote=F, sep="\t")
