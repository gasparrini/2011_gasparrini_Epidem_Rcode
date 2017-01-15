####################################################################
# Updated version of the code for the analysis in:
#
#   "The Impact of Heat Waves on Mortality"
#   Gasparrini A, Armstrong B
#   Epidemiology 2011
#   http://www.ag-myresearch.com/epidem2011.html
#
# Update: 14 March 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
####################################################################

#####################################################################
# ADDITIONAL ANALYSIS ON CHICAGO
#	RESIDUALS, PREDICTION, CORRELATIONS
#####################################################################

require(dlnm);require(mvmeta)
require(Epi);require(tsModel)
require(NMMAPSlite);require(metafor);require(foreign)
# ADD THE PACKAGE splines, NOT LOADED ANYMORE BY dlnm
require(splines)

initDB()

# LOAD AND PREPARE DATA
datatot <- readCity("chic", collapseAge = T)
datatot$tmean <- (datatot$tmpd-32)*5/9
datatot$time <- 1:nrow(datatot)
datatot$year <- as.numeric(substr(datatot$date,1,4))
datatot$month <- as.numeric(substr(datatot$date,6,7))
datatot$doy <- sequence(tapply(datatot$year,datatot$year,length))
datatot$dp01 <- filter(datatot$dptp,c(1,1)/2,side=1)
percentiles <- quantile(datatot$tmean,c(75,97:99)/100,na.rm=T)
data <- datatot[datatot$month%in%6:9,]

range <- round(range(data$tmean,na.rm=T),0)
ktemp <- range[1] + (range[2]-range[1])/4*1:3
klag <- logknots(30,3)
basis <- crossbasis(data$tmean,group=data$year,argvar=list(fun="bs",degree=3,
	knots=ktemp),arglag=list(knots=klag),lag=10)

fun.hw.thr <- function(x,thr,dur,group=NULL) {
	as.numeric(apply(Lag(x>=thr,0:(dur-1),group=group),
		1,sum,na.rm=T)>(dur-1))
}

hw <- fun.hw.thr(data$tmean,percentiles[2],2,data$year)
hw.lin <- hw
for(j in 2:10) {
	hw.lin[apply(Lag(hw,0:(j-1),group=data$year),
		1,sum,na.rm=T)==j] <- j
}

# RUN THE MODEL WITH CONTINUOUS HW DAYS
quad <- bs(hw.lin,knots=c(2,5,8),Bound=c(0,10),degree=2)
model.quad <- glm(death ~  basis + quad + dow + ns(dp01,df=3) +
	ns(year,3) + ns(doy,df=4), na.action="na.exclude",
	family=quasipoisson(), data)

#####################################################################
# RESIDUALS
############

# SIMPLE
dres <- residuals(model.quad,type="deviance")

layout(matrix(c(1,3,2,3),2,2))
hist(dres,main="Distribution of deviance residuals",xlab="Pearson residuals")
qqnorm(dres)
plot(data$date,dres,main="Deviance residuals in different years")

# STANDARDIZED

layout(matrix(1:4,2,2))
plot(model.quad)

# ALMOST THE SAME AS...
stdres <- rstandard(model.quad)
layout(matrix(c(1,3,2,3),2,2))
hist(stdres,main="Distribution of std residuals",
	xlab="Pearson residuals")
qqnorm(stdres)
plot(data$date,stdres,main="Series of standardized residuals")

hval <- hatvalues(model.quad)
layout(1)
qqnorm(stdres)
qqnorm(dres/sqrt(1-hval))

#####################################################################
# PREDICTION
##############

pred <- predict(model.quad,type="response")

# PLOT OF OBSERVED AND PREDICTED IN MAJOR HW MONTHS
restr1 <- with(data,(year==1988&month==8))
restr2 <- with(data,(year==1995&month==7))
layout(matrix(1:2,ncol=2))
plot(data$date[restr1],data$death[restr1],pch=19,
	ylab="Deaths",main="Aug 1988")
points(data$date[restr1],pred[restr1],pch=19,col=2)
legend("topright",c("Observed","Predicted"),cex=0.6,pch=19,col=1:2,inset=0.1)
plot(data$date[restr2],data$death[restr2],pch=19,
	ylab="Deaths",main="July 1995")
points(data$date[restr2],pred[restr2],pch=19,col=2)
legend("topright",c("Observed","Predicted"),cex=0.6,pch=19,col=1:2,inset=0.1)

# MEAN OBSERVED AND PREDICTED DURING SPECIFIC HW'S 
cbind(1:31,hw[restr1])
restr3 <- match(paste("1988-08-",12:18,sep=""),as.character(data$date))
cbind(1:31,hw[restr2])
restr4 <- match(paste("1995-07-",13:16,sep=""),as.character(data$date))
# DURING AUG 1988
mean(data$death[restr3]);mean(pred[restr3])
# DURING JULY 1995
mean(data$death[restr4]);mean(pred[restr4])

#####################################################################
# CORRELATION
##############

# CORRELATION WITH HW DAYS
# TOTAL
cor(data$tmean,hw)
layout(1:2)
hist(data$tmean[hw==0],col=2,main="non-HW days",
	xlim=c(5,34),xlab="Mean temperature")
hist(data$tmean[hw==1],col=2,main="HW days",
	xlim=c(5,34),xlab="Mean temperature")

# CORRELATION WITH CONSECUTIVE HW DAYS
# TOTAL
cor(data$tmean,hw.lin)
# IN HW DAYS
cor(data$tmean[hw==1],hw.lin[hw==1])

#
