# Generate two tables by site from ABIDE data:
# 1) Summarize data for subjects used in linear models
# 2) Summarize data for subjects used in meta-analysis

library(foreign)

datatable <- function(dg, roi="CB") {
	label <- names(dg)
	n.asd <- sapply(dg, function(x) nrow(x$ASD))
	n.male.asd <- sapply(dg, function(x) sum(x$ASD$SEX == 1))
	n.ctrl <- sapply(dg, function(x) nrow(x$Control))
	n.male.ctrl <- sapply(dg, function(x) sum(x$Control$SEX == 1))
	age.asd <- sapply(dg, function(x) mean(x$ASD$AGE_AT_S))
	age.ctrl <- sapply(dg, function(x) mean(x$Control$AGE_AT_S))
	age.sd.asd <- sapply(dg, function(x) sd(x$ASD$AGE_AT_S))
	age.sd.ctrl <- sapply(dg, function(x) sd(x$Control$AGE_AT_S))
	iq.asd <- sapply(dg, function(x) mean(x$ASD$FIQ2, na.rm=T))
	iq.ctrl <- sapply(dg, function(x) mean(x$Control$FIQ2, na.rm=T))
	iq.sd.asd <- sapply(dg, function(x) sd(x$ASD$FIQ2, na.rm=T))
	iq.sd.ctrl <- sapply(dg, function(x) sd(x$Control$FIQ2, na.rm=T))
	mean.asd <- sapply(dg, function(x) mean(x$ASD[[roi]]))
	mean.ctrl <- sapply(dg, function(x) mean(x$Control[[roi]]))
	sd.asd <- sapply(dg, function(x) sd(x$ASD[[roi]]))
	sd.ctrl <- sapply(dg, function(x) sd(x$Control[[roi]]))
	
	data=data.frame(label, n.asd, n.male.asd, n.ctrl, n.male.ctrl, 
			age.asd, age.ctrl, age.sd.asd, age.sd.ctrl,
			iq.asd, iq.ctrl, iq.sd.asd, iq.sd.ctrl,
			mean.asd, mean.ctrl, sd.asd, sd.ctrl)
}

# table import
data.jmp <- read.xport("cerebellum.stx")

# quality check filter
df <- subset(data.jmp, CBANALYS == "Include")
ds <- split(df, df$SITE_ID2)
dg <- lapply(ds, function(x) split(x, x$DX_GROUP))

# no quality check filter
ds.tot <- split(data.jmp, data.jmp$SITE_ID2)
dg.tot <- lapply(ds.tot, function(x) split(x, x$DX_GROUP))
data.tot <- datatable(dg.tot)

data <- datatable(dg)
data$n.asd.tot <- data.tot[rownames(data),]$n.asd
data$n.male.asd.tot <- data.tot[rownames(data),]$n.male.asd
data$n.ctrl.tot <- data.tot[rownames(data),]$n.ctrl
data$n.male.ctrl.tot <- data.tot[rownames(data),]$n.male.ctrl

write.table(data, "data-abide.txt", row.names=F, quote=F, sep="\t")


# Equalization of the mean age and the sd of age of the two groups of each cohort

# Stouffer method for combining p-values
stouffer=function(p.values, weights) {
  pnorm(sum(weights*qnorm(p.values),na.rm=TRUE)/sqrt(sum(weights[!is.na(p.values)]^2)))
}

score <- function(s2) {
  # Fisher's exact test
	ni <- sapply(s2, sapply, nrow)
	psr <- fisher.test(ni)$p.value
	
	# Student test
	pmd <- sapply(c(m=1, f=2), function(i)
	    if (rowSums(ni)[i]==0) 1 else
	    t.test(s2$ASD[[i]]$AGE_AT_S, s2$Control[[i]]$AGE_AT_S)$p.value)
	pmdc <- stouffer(pmd, rowSums(ni))
	
	# F-test of equality of variances
	pvr <- sapply(c(m=1, f=2), function(i)
	  if (rowSums(ni)[i]==0) 1 else
	  var.test(s2$ASD[[i]]$AGE_AT_S, s2$Control[[i]]$AGE_AT_S)$p.value)
	pvrc <- stouffer(pvr, rowSums(ni))
	
	min(psr, pmdc, pvrc)
}

dg2 <- lapply(dg, function(s) {
	s2 <- lapply(s, function(x) split(x, x$SEX))
	while (T) {
	  # Contingency table in function of diagnosis and sex
		ni <- sapply(s2, sapply, nrow)
		# effecive of smaller group by sex
		ns <- apply(ni, 1, min)
		# if one group by sex is too small, restore individuals for the other sex and remove individuals for that sex
		sup <- ns < 3 & ns > 0
		if (sum(sup) == 1)
		  s2 <- lapply(s, function(x) split(x, x$SEX))
		s2 <- lapply(s2, function(x) {x[sup] <- lapply(x[sup], function(y) y[0,]); x})
		ni <- sapply(s2, sapply, nrow)
		
		s <- lapply(s2, do.call, what=rbind)
		if (nrow(s$ASD) < 2 | nrow(s$Control) < 2)
			return(lapply(s, function(x) x[0,]))
		
		if (score(s2) >= .2)
			return(s)
		
		scores = lapply(c(asd=1, ctrl=2), function(i) lapply(c(m=1, f=2), function(j)
			sapply(seq(nrow(s2[[i]][[j]])), function(k) {
				s3 <- s2
				s3[[i]][[j]] <- s2[[i]][[j]][-k,]
				score(s3)
			})
		))
		
		score.max <- max(unlist(scores))
		for (i in 1:2) {
			for (j in 1:2) {
				k <- which(scores[[i]][[j]]==score.max)[1]
				if (!is.na(k)) break
			}
			if (!is.na(k)) break
		}
		
		s2[[i]][[j]] <- s2[[i]][[j]][-k,]
	}
})

data.Cbl <- na.omit(datatable(dg2, "CB"))
write.table(data.Cbl, "means-abide-Cbl.txt", row.names=F, quote=F, sep="\t")

data.Cbl_WM <- na.omit(datatable(dg2, "CBWM"))
write.table(data.Cbl_WM, "means-abide-Cbl_WM.txt", row.names=F, quote=F, sep="\t")

data.Cbl_GM <- na.omit(datatable(dg2, "CBGM"))
write.table(data.Cbl_GM, "means-abide-Cbl_GM.txt", row.names=F, quote=F, sep="\t")


