####################################################################
# Updated version of the code for the analysis in:
#
#   "The Impact of Heat Waves on Mortality"
#   Gasparrini A, Armstrong B
#   Epidemiology 2011
#   http://www.ag-myresearch.com/2011_gasparrini_epidem.html
#
# Update: 11 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2011_gasparrini_Epidem_Rcode
####################################################################

####################################################################
# NB: THE CODE HAS BEEN ADAPTED TO THE NEW VERSION OF THE R PACKAGE dlnm
#     ALSO, THE MULTIVARIATE META-ANALYSIS IS NOW CARRIED OUT INTERNALLY
#     IN R, WITHOUT THE NEED TO CALL STATA
####################################################################

require(dlnm);require(mvmeta)
require(Epi);require(tsModel)
require(NMMAPSlite);require(metafor);require(foreign)
# ADD THE PACKAGE splines, NOT LOADED ANYMORE BY dlnm
require(splines)

# CHECK VERSION OF THE PACKAGE
  if(packageVersion("dlnm")<"2.2.0")
    stop("update dlnm package to version >= 2.2.0")

# FUNCTION TO CREATE AN HEAT WAVE INDICATOR FOR A TEMPERATURE SERIES
# 	BASED ON THE THRESHOLD AND THE DURATION, BY GROUPS
fun.hw.thr <- function(x,thr,dur,group=NULL) {
	as.numeric(apply(Lag(x>=thr,0:(dur-1),group=group),
		1,sum,na.rm=T)>(dur-1))
}

# INITIALIZE THE DATASET
initDB()
cities <- listCities()

# CREATE THE MATRICES TO STORE THE RESULTS
# 	DESCRIPTIVE STATS
descr.tmean <- matrix(NA,length(cities),7, dimnames=list(cities,
	names(summary(c(1:10,NA)))))
hw.N <- matrix(NA,length(cities),6, dimnames=list(cities,
	paste("hw",rep(c(2,4),each=3),rep(c(97,98,99),2),sep=".")))
hw.cons <- matrix(NA,length(cities),4, dimnames=list(cities,
	c("N","Max",">3",">7")))
# 	REGRESSION MODELS
main.eff <- added.eff <- matrix(NA,length(cities),12, 
	dimnames=list(cities,paste("hw",rep(c(2,4),each=6),rep(c(97,98,99),
	each=2),c("est","sd"),sep=".")))
strata.eff <- matrix(NA,length(cities),5,dimnames=list(cities,1:5))
strata.vcov <- vector("list",length(cities)) ; names(strata.vcov) <- cities
quad.eff <- strata.eff
quad.vcov <- strata.vcov
# 	MEAN SUMMER TEMPERATURE
meantemp <- 0

###########################################################################

# START THE LOOP FOR CITIES
# COMPUTING TIME IS ~2-3 MIN (IN A 2.66GHz-4GBRAM PC UNDER WINDOWS), OR
# LONGER IF THE FOLDER NMMAPS STILL NEED TO BE CREATED AND THE DATA DOWNLOADED
time <- proc.time()
for(i in seq(length(cities))) {

	# LOAD AND PREPARE DATASET
	datatot <- readCity(cities[i], collapseAge = T)
	datatot$tmean <- (datatot$tmpd-32)*5/9
	datatot$time <- 1:nrow(datatot)
	datatot$year <- as.numeric(substr(datatot$date,1,4))
	datatot$month <- as.numeric(substr(datatot$date,6,7))
	datatot$doy <- sequence(tapply(datatot$year,datatot$year,length))
	datatot$dp01 <- filter(datatot$dptp,c(1,1)/2,side=1)
	percentiles <- quantile(datatot$tmean,c(75,97:99)/100,na.rm=T)
	data <- datatot[datatot$month%in%6:9,]

	# SAVE DESCRIPTIVE STATISTICS FOR TEMPERATURE
	descr.tmean[i,1:6] <- summary(data$tmean)[1:6]
	descr.tmean[i,7] <- sum(is.na(data$tmean))
	meantemp[i] <- mean(data$tmean,na.rm=T)

	# CREATE THE CROSSBASIS FOR THE MAIN TEMPERATURE-MORTALITY RELATIONSHIP
	# CENTERED ON 75TH PERCENTILE, REFERENCE VALUE FOR PREDICTED EFFECTS
	range <- round(range(data$tmean,na.rm=T),0)
	ktemp <- range[1] + (range[2]-range[1])/4*1:3
  klag <- logknots(30,3)
	basis <- crossbasis(data$tmean,group=data$year,argvar=list(fun="bs",degree=3,
		knots=ktemp),arglag=list(knots=klag),lag=10)

	#############################################################
	# FIRST ANALYSIS: INDICATOR FOR DIFFERENT HW DEFINITIONS
	#############################################################

	# HW DEFINITIONS
	hw.def <- cbind(rep(percentiles[2:4],2),rep(c(2,4),c(3,3)))
	
	# RUN THE MODEL FOR EACH DEFINITION
	for(k in 1:nrow(hw.def)) {

		# CREATE HEATWAVE INDICATOR FOR THE SPECIFIC HW DEFINITION
		hw <- fun.hw.thr(data$tmean,hw.def[k,1],hw.def[k,2],data$year)
		hw.N[i,k] <- sum(hw)

		# RUN THE MODEL
		model.first <- glm(death ~  hw + basis + dow + ns(year,3) + 
			ns(doy,df=4) + ns(dp01,df=3), family=quasipoisson(), data)
		# SAVE MAIN EFFECT
		if(sum(hw)>0) {
		tmedian <- median(data$tmean[hw==1],na.rm=T)
		pred <- crosspred(basis,model.first,at=c((range[1]+1):(range[2]-1),tmedian),
		  cen=percentiles[1])
		main.eff[i,c(k*2-1,k*2)] <- cbind(pred$allfit,
			pred$allse)[as.character(tmedian),]
		} else main.eff[i,c(k*2-1,k*2)] <- c(NA,NA)
		# SAVE ADDED EFFECT
		added.eff[i,c(k*2-1,k*2)] <- ci.lin(model.first)["hw",1:2]
	}
	
	#############################################################
	# SECOND ANALYSIS: STRATA AND QUAD SPLINE OF CONSECUTIVE HW DAYS
	#############################################################

	# CREATE HEATWAVE INDICATOR AND CONSECUTIVE TERM (97TH PERCENTILE)
	hw <- fun.hw.thr(data$tmean,percentiles[2],2,data$year)
	# CREATE HW CONSECUTIVE DAYS (UP TO 10 DAYS)
	hw.lin <- hw
	for(j in 2:10) {
		hw.lin[apply(Lag(hw,0:(j-1),group=data$year),
			1,sum,na.rm=T)==j] <- j
	}
	# SAVE STATS ON CONSECUTIVE HW DAYS
	hw.cons[i,] <- c(sum(hw),max(hw.lin),sum(hw.lin>3),sum(hw.lin>7))

	# CREATE THE STRATA OF CONSECUTIVE HW DAYS
	strata <- onebasis(c(1:10,hw.lin),fun="strata",
    breaks=c(1,2,4,6,8))[-(1:10),]
  # RUN THE MODEL
	model.strata <- glm(death ~  basis + strata + dow + 
		ns(dp01,df=3) + ns(year,3) + ns(doy,df=4),
		family=quasipoisson(), data)
	# SAVE THE RELATED COEF AND VCOV (INCLUDING MISSING)
	index1 <- grep("strata",names(coef(model.strata)))
	index2 <- (1:length(coef(model.strata)))[is.na(coef(model.strata))]
	index <- index1[!index1%in%index2]
	strata.eff[i,!index1%in%index2] <- ci.lin(model.strata)[index,1]
	strata.vcov[[i]] <- matrix(NA,length(index1),length(index1))
	strata.vcov[[i]][!index1%in%index2,!index1%in%index2] <- 
		vcov(model.strata)[index,index]

	# CREATE THE SPLINE OF CONSECUTIVE HW DAYS
	quad <- onebasis(hw.lin,fun="bs",degree=2,knots=c(2,5,8),Bound=c(0,10))
	# RUN THE MODEL
	model.quad <- glm(death ~  basis + quad + dow + ns(dp01,df=3) +
		ns(year,3) + ns(doy,df=4),family=quasipoisson(), data)
	# SAVE THE RELATED COEF AND VCOV (INCLUDING MISSING)
	index1 <- grep("quad",names(coef(model.quad)))
	index2 <- (1:length(coef(model.quad)))[is.na(coef(model.quad))]
	index <- index1[!index1%in%index2]
	quad.eff[i,!index1%in%index2] <- ci.lin(model.quad)[index,1]
	quad.vcov[[i]] <- matrix(NA,length(index1),length(index1))
	quad.vcov[[i]][!index1%in%index2,!index1%in%index2] <- 
		vcov(model.quad)[index,index]
}
proc.time()-time

###########################################################################

###############################
# MULTIVARIATE META-ANALYSIS
###############################

mv.quad <- mvmeta(quad.eff~1,S=quad.vcov,method="mm")
summary(mv.quad)

mv.strata <- mvmeta(strata.eff~1,S=strata.vcov,method="mm")
summary(mv.strata)

###########################################################################

###############################
# DESCRIPTIVE STATISTICS
###############################

# SUMMARY FOR TMEAN
summary(descr.tmean[,c("Mean","NA's")])

# TOTAL NUMBER OF HW DAYS UNDER DIFFERENT HW DEFINITIONS
summary(hw.N)

# CONSECUTIVE HW DAYS (WITH 97TH PERCENTILE)
# % OF CITIES WITH MAX LENGTH >7 AND >9
sum(hw.cons[,"Max"]>6)/nrow(hw.cons)*100
sum(hw.cons[,"Max"]>9)/nrow(hw.cons)*100
# % OF CONSECUTIVE HW DAYS ABOVE 3 AND 7
colSums(hw.cons[,c(">3",">7")])/sum(hw.cons[,"N"])*100

###############################
# FIRST ANALYSIS
###############################

label <- paste("hw",rep(c(2,4),each=3),rep(c(97,98,99),2),sep=".")
table1 <- matrix(NA,6,7,dimnames=list(label,
	c("N comm","Est.main","95%CI.main","P-het.added","Est.added",
		"95%CI.added","P-het.added")))

for(i in 1:6) {	 

	# SET TO MISSING IF NO ESTIMATE FOR ADDED EFFECT
	added.eff[added.eff[,2*i]==0,c(2*i-1,2*i)] <- NA
	main.eff[is.na(added.eff[,2*i]),c(2*i-1,2*i)] <- NA

	# RUN THE META-ANALYSIS
	pool.main <- rma.uni(yi=main.eff[,2*i-1],sei=main.eff[,2*i])
	pool.added <- rma.uni(yi=added.eff[,2*i-1],sei=added.eff[,2*i])
  
	# FILL TABLE1
	table1[i,] <- c(sum(!is.na(added.eff[,2*i-1])),
		round(exp(pool.main$b)*100-100,1),
		paste(round(exp(pool.main$b-1.96*pool.main$se)*100-100,1),"to",
		round(exp(pool.main$b+1.96*pool.main$se)*100-100,1)),
		round(pool.main$QEp,3),
		round(exp(pool.added$b)*100-100,1),
		paste(round(exp(pool.added$b-1.96*pool.added$se)*100-100,1),"to",
		round(exp(pool.added$b+1.96*pool.added$se)*100-100,1)),
		round(pool.added$QEp,3))
}

# TABLE 1 IN THE MANUSCRIPT
table1

###############################
# SECOND ANALYSIS
###############################

# NB: THE FIGURE IS SLIGHTY DIFFERENT THAN THAT ILLUSTRATED IN THE PAPER:
#     - THE CURVE AND CI ARE AS A DIFFERENT ESTIMATOR IS USED IN THE R PACKAGE
#       mvmeta IF COMPARED TO THAT USED IN STATA AT THE TIME OF THE PUBLICATION
#     - THE STRATA STEPS ARE DIFFERENT FOR AN ERROR IN THE ORIGINAL CODE

# NB: USE OF THE FUNCTIONS IN THE PACKAGE dlnm FOR SIMPLIFYING THE PLOTTING

# CREATE THE SAME TRANSFORMATION USED IN ESTIMATION
x.quad <- onebasis(0:100/10,fun="bs",knots=c(2,5,8),degree=2,
  Bound=c(0,10))
x.strata <-onebasis(0:20/2,fun="strata",breaks=c(1,2,4,6,8))

# PREDICT THE CURVES USING THE ESTIMATED AVERAGE COEF AND VCOV FOR THE SPLINE
pred.quad <- crosspred(x.quad,coef=coef(mv.quad),vcov=vcov(mv.quad),
  model.link="log",by=0.1,from=0,to=10,cen=0)
pred.strata <- crosspred(x.strata,coef=coef(mv.strata),vcov=vcov(mv.strata),
  model.link="log",by=0.1,from=0,to=10)

# PLOT
plot(pred.quad,ylim=c(0.95,1.10),ylab="Percent change %",yaxt="n",col=1,
	xlab="Number of consecutive heat-waves days",frame.plot=F)
axis(2,labels=-1:2*5,at=0.95+0:3*0.05)
lines(pred.strata,type="s",col=1,lty=2)

#
