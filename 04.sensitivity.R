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

#############################################################
# - SENSITIVITY ANALYSIS FOR INDICATOR VARIABLE
#############################################################

require(dlnm);require(Epi);require(tsModel)
require(NMMAPSlite);require(metafor);require(foreign)

# FUNCTION TO CREATE AN HEAT WAVE INDICATOR FOR A TEMPERATURE SERIES
# 	BASED ON THE THRESHOLD AND THE DURATION, BY GROUPS
fun.hw.thr <- function(x,thr,dur,group=NULL) {
	as.numeric(apply(Lag(x>=thr,0:(dur-1),group=group),
		1,sum,na.rm=T)>(dur-1))
}

# INITIALIZE THE DATASET
initDB()
cities <- listCities()

sensitivity <- function(thr,dur,seas,temp,lag) {

	# EMPTY MATRIX
	added.eff <- matrix(NA,length(cities),2,
		dimnames=list(cities,c("est","sd")))

	#############################################################

	# START THE LOOP FOR CITIES
	for(i in seq(length(cities))) {

		# LOAD AND PREPARE DATASET
		datatot <- readCity(cities[i], collapseAge = T)
		datatot$tmean <- (datatot$tmpd-32)*5/9
		datatot$time <- 1:nrow(datatot)
		datatot$year <- as.numeric(substr(datatot$date,1,4))
		datatot$month <- as.numeric(substr(datatot$date,6,7))
		datatot$doy <- sequence(tapply(datatot$year,datatot$year,length))
		datatot$dp01 <- filter(datatot$dptp,c(1,1)/2,side=1)
		percentiles <- quantile(datatot$tmean,c(75,thr)/100,na.rm=T)
		data <- datatot[datatot$month%in%6:9,]
	
		# CREATE THE CROSSBASIS
		range <- round(range(data$tmean,na.rm=T),0)
		ktemp <- range[1] + (range[2]-range[1])/(temp-2)*1:(temp-3)
    klag <- logknots(30,3)
    basis <- crossbasis(data$tmean,group=data$year,argvar=list(fun="bs",
      degree=3,knots=ktemp),arglag=list(knots=klag),lag=10)

		# CREATE HEATWAVE INDICATOR FOR THE SPECIFIC HW DEFINITION
		hw <- fun.hw.thr(data$tmean,percentiles[2],dur,data$year)

		# RUN THE MODEL
		model.first <- glm(death ~  hw + basis + dow + ns(year,3) + 
			ns(doy,df=seas) + ns(dp01,df=3), family=quasipoisson(), data)
		added.eff[i,] <- ci.lin(model.first)["hw",1:2]
	}
	added.eff[added.eff[,1]==0,] <- NA
	pool.added <- rma.uni(yi=added.eff[,1],sei=added.eff[,2])
	eff <- as.numeric(c(round(exp(pool.added$b)*100-100,1),
		paste(round(exp(pool.added$b-1.96*pool.added$se)*100-100,1)),
		paste(round(exp(pool.added$b+1.96*pool.added$se)*100-100,1))))
	names(eff) <- c("Eff","lowCI","HighCI")
	eff
}

#############################################################

table2 <- matrix(NA,8,7)
colnames(table2) <- c("seas","temp","lag","2-97Eff","2-97CI",
	"4-99Eff","4-99CI")
value <- cbind(c(2,6,rep(4,6)),c(6,6,4,4,7,6,6,6),c(5,5,1,5,5,1,3,6))
table2[,1:3] <- value

time <- proc.time()
for(i in 1:nrow(table2)) {
	eff <- sensitivity(97,2,value[i,1],value[i,2],value[i,3])
	table2[i,4:5] <- c(eff[1],paste("(",eff[2]," to ",eff[3],")",sep=""))
	eff <- sensitivity(99,4,value[i,1],value[i,2],value[i,3])
	table2[i,6:7] <- c(eff[1],paste("(",eff[2]," to ",eff[3],")",sep=""))
}
proc.time()-time

table2

#
