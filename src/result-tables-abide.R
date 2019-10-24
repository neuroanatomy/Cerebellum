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

script.dir <- get.script.dir()
meta.dir <- file.path(base.dir, "data", "meta-analysis")

regions <- c("Cerebellum", "Cerebellum WM", "Cerebellum GM")

suffixes  <- c("abide-Cbl", "abide-Cbl_WM", "abide-Cbl_GM")
names(suffixes) <- regions

values_es <- c("Nstud", "TE", "lTE", "uTE", "p.rand", "I2", "lI2", "uI2", "bAge", "lbAge", "ubAge", "bIQ", "lbIQ", "ubIQ", "rI2", "rI2p")
SMD_es <- data.frame(matrix(nrow=length(regions), ncol=length(values_es)))
names(SMD_es) <- values_es
row.names(SMD_es) <- regions
LVR_es <- data.frame(matrix(nrow=length(regions), ncol=length(values_es)))
names(LVR_es) <- values_es
row.names(LVR_es) <- regions

values_bias <- c("Nstud", "fsign", "apower", "pEgger", "p_right", "p_power", "power_est", "power_est.lb", "power_est.ub")
SMD_bias <- data.frame(matrix(nrow=length(regions), ncol=length(values_bias)))
names(SMD_bias) <- values_bias
row.names(SMD_bias) <- regions
LVR_bias <- data.frame(matrix(nrow=length(regions), ncol=length(values_bias)))
names(LVR_bias) <- values_bias
row.names(LVR_bias) <- regions

results <- list()

for (reg in regions) {
	print(reg)
	suffix <- suffixes[reg]
	
	source(paste(script.dir, "meta-exec.R", sep="/"))
	
	SMD_es[reg,] <- c(meta_SMD$k, meta_SMD$TE.random, meta_SMD$lower.random, meta_SMD$upper.random, meta_SMD$pval.random,
					 meta_SMD$I2, meta_SMD$lower.I2, meta_SMD$upper.I2, reg_SMD$b[2], reg_SMD$ci.lb[2], reg_SMD$ci.ub[2],
					 reg_SMD$b[3], reg_SMD$ci.lb[3], reg_SMD$ci.ub[3], reg_SMD$I2, reg_SMD$QEp)
	SMD_bias[reg,] <- c(meta_SMD$k, mean(p_SMD<alpha), mean(power.reg_SMD), bias_SMD$p.value, pcurve_SMD$p_right, pcurve_SMD$p_power, pcurve_SMD$power_est, pcurve_SMD$power.ci.lb, pcurve_SMD$power.ci.ub)
	
	LVR_es[reg,] <- c(meta_F$k, meta_F$TE.random, meta_F$lower.random, meta_F$upper.random, meta_F$pval.random,
					  meta_F$I2, meta_F$lower.I2, meta_F$upper.I2, reg_F$b[2], reg_F$ci.lb[2], reg_F$ci.ub[2],
					  reg_F$b[3], reg_F$ci.lb[3], reg_F$ci.ub[3], reg_F$I2, reg_F$QEp)
	LVR_bias[reg,] <- c(meta_F$k, mean(p_F<alpha), mean(power.reg_F), bias_F$p.value, pcurve_F$p_right, pcurve_F$p_power, pcurve_F$power_est, pcurve_F$power.ci.lb, pcurve_F$power.ci.ub)
	
	results[[reg]] <- list(meta_SMD=meta_SMD, reg_SMD=reg_SMD, meta_F=meta_F, reg_F=reg_F)
}

# Results format

format_es <- function(results) {
	rf <- results
	rf[, c("Nstud", "rI2")] <- format(rf[, c("Nstud", "rI2")], trim=T, nsmall=0, digits=0)
	rf[, c("TE", "lTE", "uTE")] <- format(rf[, c("TE", "lTE", "uTE")], trim=T, nsmall=2, digits=0)
	rf[, c("I2", "lI2", "uI2")] <- format(rf[, c("I2", "lI2", "uI2")] * 100, trim=T, nsmall=0, digits=0)
	rf[, c("p.rand", "bAge", "lbAge", "ubAge", "bIQ", "lbIQ", "ubIQ", "rI2p")] <- format(rf[, c("p.rand", "bAge", "lbAge", "ubAge", "bIQ", "lbIQ", "ubIQ", "rI2p")], trim=T, nsmall=3, digits=0)
	tab <- data.frame(rf$Nstud, sprintf("%s (%s, %s)", rf$TE, rf$lTE, rf$uTE), rf$p.rand, sprintf("%s%% (%s%%, %s%%)", rf$I2, rf$lI2, rf$uI2),
					  sprintf("%s (%s, %s)", rf$bAge, rf$lbAge, rf$ubAge), sprintf("%s (%s, %s)", rf$bIQ, rf$lbIQ, rf$ubIQ),
					  sprintf("%s%%\np = %s", rf$rI2, rf$rI2p))
	
	names(tab) <- c("N Studies", "Effect Size", "p-value", "I2", "Age impact", "IQ impact", "residual I2")
	row.names(tab) <- regions
	tab
}

format_bias <- function(results) {
	rf <- results
	rf$Nsign <- format(rf$fsign*rf$Nstud, trim=T, nsmall=0, digits=0)
	rf$Nstud <- format(rf$Nstud, trim=T, nsmall=0, digits=0)
	rf[, c("fsign", "apower", "power_est", "power_est.lb", "power_est.ub")] <- format(rf[, c("fsign", "apower", "power_est", "power_est.lb", "power_est.ub")] * 100, trim=T, nsmall=0, digits=0)
	rf[, c("pEgger", "p_right", "p_power")] <- format(rf[, c("pEgger", "p_right", "p_power")], trim=T, nsmall=3, digits=0)
	tab <- data.frame(sprintf("%s%% (%s/%s)", rf$fsign, rf$Nsign, rf$Nstud), sprintf("%s%%", rf$apower),
					  sprintf("p = %s", rf$pEgger), sprintf("p = %s", rf$p_right), sprintf("p = %s", rf$p_power),
					  sprintf("%s%% (%s%%, %s%%)", rf$power_est, rf$power_est.lb, rf$power_est.ub))
	
	names(tab) <- c("Significative studies rate", "Mean achieved power", "Egger's test", "probing value", "inadequate probing value", "infered power")
	row.names(tab) <- regions
	tab
}

tab_se <- format_es(SMD_es)
tab_sb <- format_bias(SMD_bias)
tab_le <- format_es(LVR_es)
tab_lb <- format_bias(LVR_bias)

op <- par(mfrow=c(2,2))
for (reg in regions) {
  funnel(results[[reg]]$meta_SMD, main=reg)
}
par(op)