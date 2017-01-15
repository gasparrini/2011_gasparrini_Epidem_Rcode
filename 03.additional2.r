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

#####################################################################
# ADDITIONAL ANALYSIS ON THE WHOLE DATASET
# 	CORRELATIONS AND MULTI-COLLINEARITY
#####################################################################

require(dlnm);require(mvmeta)
require(Epi);require(tsModel)
require(NMMAPSlite);require(metafor);require(foreign)
# ADD THE PACKAGE splines, NOT LOADED ANYMORE BY dlnm
require(splines)

# FUNCTION TO CREATE AN HEAT WAVE INDICATOR FOR A TEMPERATURE SERIES
# 	BASED ON THE THRESHOLD AND THE DURATION, BY GROUPS
fun.hw.thr <- function(x,thr,dur,group=NULL) {
	as.numeric(apply(Lag(x>=thr,0:(dur-1),group=group),
		1,sum,na.rm=T)>(dur-1))
}

# INITIALIZE THE DATASET
initDB()
cities <- listCities()

# CREATE VECTORS TO STORE THE RESULTS
R2ind <- R2lin <- corrind <- corrlin <- rep(NA,length(cities))
names(R2ind) <- names(R2lin) <- names(corrind) <- names(corrlin) <- cities

# START THE LOOP FOR CITIES
# COMPUTING TIME IS ~1-2 MIN (IN A 2.66GHz-4GBRAM PC UNDER WINDOWS), OR
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

	# CREATE THE CROSSBASIS FOR THE MAIN TEMPERATURE-MORTALITY RELATIONSHIP
	# CENTERED ON 75TH PERCENTILE, REFERENCE VALUE FOR PREDICTED EFFECTS
	range <- round(range(data$tmean,na.rm=T),0)
	ktemp <- range[1] + (range[2]-range[1])/4*1:3
  klag <- logknots(30,3)
  basis <- crossbasis(data$tmean,group=data$year,argvar=list(fun="bs",degree=3,
    knots=ktemp),arglag=list(knots=klag),lag=10)

	# HW DEFINITIONS (INDICATOR AND CONSECUTIVE) BASED ON 97TH-2DAYS
	hw <- fun.hw.thr(data$tmean,percentiles[2],2,data$year)
	hw.lin <- hw
	for(j in 2:10) {
		hw.lin[apply(Lag(hw,0:(j-1),group=data$year),
			1,sum,na.rm=T)==j] <- j
	}
	
	corrind[i] <- cor(hw,data$tmean,use="pairwise.complete.obs")
	corrlin[i] <- cor(hw.lin,data$tmean,use="pairwise.complete.obs")
	R2ind[i] <- summary(lm(hw ~ basis))$r.squared
	R2lin[i] <- summary(lm(hw.lin ~ basis))$r.squared
}
proc.time()-time

# CORRELATIONS
mean(corrind) ; range(corrind)
mean(corrlin) ; range(corrlin)

# MULTI-COLLINEARITY
mean(R2ind) ; range(R2ind)
mean(R2lin) ; range(R2lin)

#
