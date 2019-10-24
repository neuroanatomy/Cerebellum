library(stringr)  #Library to process string variables (text of the entered tests)

#Function 5 Stouffer test for a vector of pp-values
stouffer=function(pp) sum(qnorm(pp),na.rm=TRUE)/sqrt(sum(!is.na(pp)))

# Vectorizable version of uniroot
multiroot <- function(f, interval, ..., lower=min(interval), upper=max(interval), maxiter=1000) {
	xl <- lower
	yl <- f(xl, ...)
	xu <- upper
	yu <- f(xu, ...)
	if (any(sign(yl)==sign(yu)))
		stop("interval is not cut by zero")
	
	dir <- yl < yu
	x1 <- dir * xl + (!dir) * xu
	x2 <- dir * xu + (!dir) * xl
	y1 <- dir * yl + (!dir) * yu
	y2 <- dir * yu + (!dir) * yl
	
	epsilon <- 1e-1
	error_tol <- 1e-5
	
	for (i in 1:maxiter) {
		x <- (-y1 / ((y2 - y1) * (1+2*epsilon)) + epsilon / (1+2*epsilon)) * (x2 - x1) + x1
		x[is.nan(x) | x==Inf] <- x1[is.nan(x) | x==Inf]
		y <- f(x, ...)
		if (all((y == y1 | y == y2) & abs(y) < error_tol))
			break
		x1[y<=0] <- x[y<=0]
		x2[y>0] <- x[y>0]
		y1[y<=0] <- y[y<=0]
		y2[y>0] <- y[y>0]
	}
	x
}

#t-test ncp for given power
get.ncp.t =function(df, power)   {
	ct <- qt(p=.975, df=df) # critical t-value
	#Find noncentrality parameter (ncp) that leads requested power to obtain ct
	f <- function(delta) pt(ct, df = df, ncp = delta, lower.tail=F) + pt(-ct, df = df, ncp = delta) - power
	multiroot(f, lower=rep(0, length(df)), upper=rep(37.62, length(df)))
}

#F-test ratio for given power
get.ratio.f =function(df1,df2, power)   {
	cf <- qf(p=alpha/2, df1, df2, lower.tail=F) # critical F-value
	fpower <- qf(p=power, df1, df2, lower.tail=F)
	cf/fpower
}

pcurve <- function(tests) {
	
	par(mfrow=c(1,1))
	alpha=0.05
	
	# Set up for two sided t-test and two sided F-test of equality of variances
	# Note: power calculation are one sided, but with alpha set up for two sided case
	
	#1. create empty vectors for
	#1.1 pp-values
	t.ppr=f.ppr=c();      	#right skew
	t.ppl=f.ppl=c();       	#left
	t.pp33=f.pp33=c();   	#33%
	
	# Parse tests vector
	#2 Split tests into t and F
	
	#2.1 Turn everything to lower case
	tests=tolower(tests)
	
	#2.2 Extract the type of test (stat={t,F})
	stat=substring(tests,1,1)
	
	#2.3 Split vector of tests into these
	#get the t-tests
	t.text=subset(tests,stat=="t")
	#get the f-tests
	f.text=subset(tests,stat=="f")
	
	#3 Get d.f. for the tests 
	#3.1 t-test
	#find the 2nd parenthesis
	t.par=str_locate(t.text,"\\)")[,1]
	#Get the d.f. between both parenthesis
	t.df=as.numeric(substring(t.text,3,t.par -1))
	
	#3.2 f-test
	#find the comma
	f.comma=str_locate(f.text,",")[,1]
	#find the 2nd parenthesis
	f.par=str_locate(f.text,"\\)")[,1]
	#Get the df1  (between "(" and ","
	f.df1=as.numeric(substring(f.text,3,f.comma -1))
	#Get the df2  (between "," and ")"
	f.df2=as.numeric(substring(f.text,f.comma +1,f.par -1))
	
	#4 Get the test values
	#4.1 Find the "=" sign
	t.eq=str_locate(t.text,"=")[,1]
	f.eq=str_locate(f.text,"=")[,1]
	
	t.value=as.numeric(substring(t.text,t.eq+1))
	f.value=as.numeric(substring(f.text,f.eq+1))
	
	#5 Keep significant p-values
	#Compute p-values
	t.p=2*(1-pt(abs(t.value),df=t.df))
	f.p=1-abs(pf(f.value,df1=f.df1,df2=f.df2)-0.5)*2
	
	#Subset statistics and d.f.
	#ts
	t.value.sig=subset(t.value,t.p<.05)
	t.df.sig   =subset(t.df,   t.p<.05)
	t.p.sig    =subset(t.p,    t.p<.05)
	#fs
	# inverse f-values inferior to one
	f.value.sig=c(subset(f.value,f.p<.05 & f.value>=1), subset(1/f.value,f.p<.05 & f.value<1))
	f.df1.sig  =c(subset(f.df1,  f.p<.05 & f.value>=1), subset(f.df2, f.p<.05 & f.value<1))
	f.df2.sig  =c(subset(f.df2,  f.p<.05 & f.value>=1), subset(f.df1, f.p<.05 & f.value<1))
	f.p.sig    =subset(f.p,    f.p<.05)
	
	#All significant p-values
	all.p.sig=c(t.p.sig, f.p.sig)
	#Number of significant results
	ktot=length(all.p.sig)
	if(ktot==0) return(list(p_right=NA, p_left=NA, p_power=NA, power_est=NA, power.ci.lb=NA, power.ci.ub=NA))
	#Number of non-signifcant results in p-curve
	kns=length(tests)-ktot
	
	#6 Compute pp-values
	#6.1 For t-values
	if (length(t.value.sig)>0)  #if nonempty compute pp-values
	{
		#skew
		t.ppr=t.p.sig*(1/.05)               #pp-value for right-skew
		t.ppl=1-t.ppr                         #pp-value for left-skew
		#33%power
		#Find the ncp (uses function from top)
		t.ncp33=get.ncp.t(t.df.sig, 1/3)
		#Using the ncp33 compute pp33.
		t.pp33 <- 1 - 3 * (pt(abs(t.value.sig), df=t.df.sig, ncp=t.ncp33, lower.tail=F) + pt(-abs(t.value.sig), df=t.df.sig, ncp=t.ncp33))
	}
	
	#6.2 For F-values
	if (length(f.value.sig)>0)  #if nonempty compute pp-values
	{
		f.ppr=f.p.sig*(1/.05)             #pp-value for right-skew 
		f.ppl=1-f.ppr                     #pp-value for left-skew
		f.r33 <- get.ratio.f(df1=f.df1.sig, df2=f.df2.sig, power=1/3)
		f.pp33 <- 1 - 3 * pf(f.value.sig/f.r33, df1=f.df1.sig, df2=f.df2.sig, lower.tail=F)
	}
	
	#7 STOUFFER: Overall tests aggregating pp-values 
	#7.1 Convert pp-values to Z scores, aggregate them and divide by sqrt(ktot)
	Zppr =sum(qnorm(c(t.ppr,  f.ppr )))/sqrt(ktot)          #right skew
	Zppl =sum(qnorm(c(t.ppl,  f.ppl )))/sqrt(ktot)          #left skew
	Zpp33=sum(qnorm(c(t.pp33, f.pp33)))/sqrt(ktot)          #33%
	
	#7.2 Compute overall p-values
	p.Zppr =pnorm(Zppr)
	p.Zppl =pnorm(Zppl)
	p.Zpp33=pnorm(Zpp33)
	
	#13.POWER ESTIMATION
	
	# 13.3 CREATE PP-VALUES FOR EACH OF THE TWO DISTRIBUTIONS FOR HOW WELL A GIVEN POWER_EST FITS 
	powerfit.t=function(t_obs, df_obs, power_est) {
		ncp_est=get.ncp.t(df=df_obs,power=power_est)               #find ncp for each  that gives each test power.k
		p_larger=pt(abs(t_obs),df=df_obs,ncp=ncp_est, lower.tail=F) +
		  pt(-abs(t_obs),df=df_obs,ncp=ncp_est)     #prob t>tobs given ncp_est
		ppr=1-p_larger/power_est                                   #condition on p<.05
		return(ppr)
	}
	
	powerfit.f=function(f_obs, df1_obs, df2_obs, power_est) {
		ratio_est=get.ratio.f(df1=df1_obs, df2=df2_obs,power=power_est)      #find ncp for each  that gives each test power.k
		p_larger=pf(f_obs/ratio_est,df1=df1_obs,df2=df2_obs, lower.tail=F)   #prob t>tobs given ncp_est
		ppr=1-p_larger/power_est                                             #condition on p<.05
		return(ppr)
	}
	
	#13.4  STACK-UP ALL THE PP-VALUES INTO A VECTOR AND COMPARE THEM TO UNIFORM DISTRIBUTION USING KOLMOGOROV-SMIRNOV TEST
	
	powerfit=function(power_est)
	{
		ppr.all=c()
		#for each kind of test, check if there are any significant values, if there are, add ppr to overall ppr
		if (length(t.value.sig)>0) ppr.all=c(ppr.all, powerfit.t(t_obs=t.value.sig, df_obs=t.df.sig, power_est=power_est))
		if (length(f.value.sig)>0) ppr.all=c(ppr.all, powerfit.f(f_obs=f.value.sig, df1_obs=f.df1.sig, df2_obs=f.df2.sig, power_est=power_est))
		# KSD=ks.test(ppr.all,punif)$statistic                #KS test on the resulting pprs
		# return(KSD)
		return(stouffer(ppr.all))
	}
	
	#13.5 COMPUTE FIT FOR EACH LEVEL OF POWER, AND PLOT IT

	# Fit will be evaluated at every possible value of power between 5.1% and 99% in steps of 1%, stored in fit()
	fit=c()                                        #Create empty vector
	fit=abs(powerfit(.051))                        #First evaluate fit for power of 5.1%, the lowest one can get for non-directional tests like x2 and F
	for (i in 6:99)   fit=c(fit,abs(powerfit(i/100))) #Now do 6% to 99%
	# Find the minimum
	mini=match(min(fit),fit)       #which ith power level considered leads to best estimate
	hat=(mini+4)/100               #convert that into the power level, the ith value considered is (5+ith)/100
	#Plot results
	#create the x-axis
	x.power=seq(from=5,to=99)/100
	#Draw the line
	par(mar=c(5.1,8.1,4.1,2.1))  #Margins
	plot(x.power,fit,xlab="Underlying Power", ylab="",ylim=c(-.15,max(fit)), main="")
	#Figure title
	mtext("Estimating underlying statistical power",side=3,line=1.75,cex=1.5,at=0.4)
	mtext("(Plot should be V shaped, or a smooth line to 99%; else don't trust estimate)",col='red',side=3,line=.5,cex=1,at=0.4)
	#Make red dot at the estimate
	points(hat,min(fit,na.rm=TRUE),pch=19,col="red",cex=2)
	#Put a label with the estimate value
	sign="="
	if (hat<.06) sign="<"
	text(min(.5,max(.28,hat)),min(fit,na.rm=TRUE)-.15,paste0("Estimated Power ",sign," ",hat*100,"%"))
	#Label the y-axis
	mtext(c("Good","Bad"),side=2,line=3,at=c(0,max(fit)),las=1,cex=1.25,col=c("blue","red"))
	mtext("Fit for observed p-curve",side=2,line=6.5,cex=1.5)
	mtext("(Stouffer test for null of power in x-axis)\n|Z-score|",side=2,line=4.5,col="gray")
	
	
	#4.3 Confidence interval for power estimate
	#4.3.1 Function get.power_pct(pct) 
	get.power_pct =function(pct)   {
	  #Function that finds power that gives p-value=pct for the Stouffer test 
	  #for example, get.power_pct(.5) returns the level of power that leads to p=.5  for the stouffer test.
	  #half the time we would see p-curves more right skewed than the one we see, and half the time
	  #less right-skewed, if the true power were that get.power_pct(.5). So it is the median estimate of power
	  #similarliy, get.power_pct(.1) gives the 10th percentile estimate of power...
	  #Obtain the normalized equivalent of pct, e.g., for 5% it is -1.64, for 95% it is 1.64
	  z=qnorm(pct)  #convert to z because powerfit() outputs a z-score. 
	  #Quantify gap between computed p-value and desired pct
	  error = function(power_est, z)  powerfit(power_est) - z
	  #Find the value of power that makes that gap zero, (root)
	  return(uniroot(error, c(.0501, .99),z)$root)   }
	
	#4.3.2 Boundary conditions (when the end of the ci=5% or 99% we cannot use root to find it, 
	#use boundary value instead)
	
	#Boundary conditions
	p.power.05=pnorm(powerfit(.051)) #Proability p-curve would be at least at right-skewed if power=.051
	p.power.99=pnorm(powerfit(.99))  #Proability p-curve would be at least at right-skewed if power=.99
	
	#4.3.3 Find lower end of ci
	#Low boundary condition? If cannot reject 5% power, don't look for lower levels, use 5% as the end 
	if (p.power.05<=.95) power.ci.lb=.05   
	#High boundary condition? If we reject 99%, from below dont look for higher power, use 99% as the low end
	if (p.power.99>=.95) power.ci.lb=.99   
	#If low bound is higher than 5.1% power and lower than 99% power, estimate it, find interior solution
	if (p.power.05>.95 && p.power.99<.95)  power.ci.lb=get.power_pct(.95)
	
	
	#4.3.4 Higher end of CI
	#If we reject 5% power from below, 5% is above the confidence interval, use 5% as the upper end of the confidence interval
	if (p.power.05<=.05) power.ci.ub=.05
	#If we do not reject that 99% power, don't look higher, use 99% as the higher end 
	if (p.power.99>=.05) power.ci.ub=.99
	#If the the upper bound is between 5% and 99%, find it
	if (p.power.05>.05 && p.power.99<.05) power.ci.ub=get.power_pct(.05)
	
	hist <- hist(all.p.sig, breaks=seq(0, 0.05, 0.01), plot=F)
	plot(hist$breaks[-1], hist$density, ylim=c(0,105), type="o", col="blue", pch=19,
	     xlab="p-value", lwd=2, ylab="frequency", yaxt="n")
	title(main="P-curve")
	#y-axis value labels
	y_=c("0%","25%","50%","75%","100%")
	y=c(0,25,50,75,100)
	axis(2,at=y,labels=y_)
	
	list(p_right=p.Zppr, p_left=p.Zppl, p_power=p.Zpp33, power_est=hat, power.ci.lb=power.ci.lb, power.ci.ub=power.ci.ub)
}