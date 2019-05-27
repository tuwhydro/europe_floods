##############################################################################
#This script is meant to reproduce the analysis in the manuscript: 
#Changing climate both increases and decreases European river floods
#Bloeschl et al. (2019)
##############################################################################
#reading in the data
data_europe = read.csv2("europe_data.csv",
                        header=T,stringsAsFactors = F)
##############################################################################
#fitting the trend and evaluating the significance
##############################################################################
#sen-slope and mann-kendall

#two-sided test for trend
mann_kendall_trend <- function(dframe,discharge_name="Qmax",year_name="year",continuity=TRUE) {
  
  #Remove NAs
  if(length(which(is.na(dframe[,discharge_name])))>0) {
    dframe <- dframe[-which(is.na(dframe[,discharge_name])),]
  }
  
  #mann-kendall
  mk <- NA
  for(i in c(1:(nrow(dframe)-1) )) {
    for(j in c((i+1):nrow(dframe))) {
      mk <- c(mk,sign(dframe[j,discharge_name]-dframe[i,discharge_name]) )
    }
  }
  mk <- mk[-1]
  mk <- sum(mk,na.rm=T)
  
  #for variance, need number of ties
  series <- dframe[,discharge_name]
  n <- length(series)
  if(any(duplicated(series))) {
    #indices of elements, that occur multiple times
    duplicate_indices <- which(duplicated(series))
    #number of elements with multiple occurences
    p <- length(unique(series[duplicate_indices]))
    #first occurrence of element with multiple occurrences
    t_index <- sapply(unique(series[duplicate_indices]), function(k) min(which(series==k)) )
    #number of occurrences of a duplicated element
    t <- sapply(unique(series[duplicate_indices]), function(k) length(which(series==k)))
    
    var_mk <- (n*(n-1)*(2*n+5) - sum(sapply(t, function(j) j*(j-1)*(2*j+5))) )/18
  }
  #no ties
  else if(!any(duplicated(series))) {
    var_mk <- (n*(n-1)*(2*n+5))/18
  }
  #actual test-statistic
  if(continuity==FALSE) {
    Z <- mk/sqrt(var_mk)
  }
  else if(continuity==TRUE) {
    Z <- sign(mk)*(abs(mk)-1)/sqrt(var_mk)
  }
  #p-value with normal distribution, assuming two-sided test here
  p_val <- (1-pnorm(abs(Z)))*2
  
  return(c(sign_sum=mk,test_stat=Z,p_val=p_val))
}

#calculating sen-slope on original scale and percentage of mean per decade
sen_slope <- function(dframe,discharge_name="Qmax",year_name="year") {
  
  #check for duplicated values 
  if(length(which(duplicated(dframe)))>0) {
    dframe <- dframe[-which(duplicated(dframe)),]
  }
  #check for NAs
  if(length(which(is.na(dframe[,discharge_name])))>0) {
    dframe <- dframe[-which(is.na(dframe[,discharge_name])),]
  }
  #slopes
  mk <- 0
  for(i in c(1:(nrow(dframe)-1) )) {
    for(j in c((i+1):nrow(dframe))) {
      mk <- c(mk,ifelse(!is.na((dframe[j,discharge_name]-dframe[i,discharge_name])/(dframe[j,year_name]-dframe[i,year_name]) ),
                        (dframe[j,discharge_name]-dframe[i,discharge_name])/(dframe[j,year_name]-dframe[i,year_name]) ,
                        NA))
      
    }
  }
  mk <- mk[-1]
  sen <- median(mk,na.rm=TRUE)
  
  mean_flow <- mean(dframe[,discharge_name],na.rm=TRUE)
  #slope as percentage if mean per decade
  sen_pct_of_mean_per_dec = (sen/mean_flow)*100*10
  return(c(sen=sen,sen_pct_of_mean_per_dec=sen_pct_of_mean_per_dec))
}

#applying to data
data_europe$mk_pval = rep(NA,nrow(data_europe))
data_europe$sen_original = rep(NA,nrow(data_europe))
data_europe$pct_of_mean_per_decade = rep(NA,nrow(data_europe))

for(j in which(data_europe$F1==1) ) {
  #at least 3 observations
  if(sum(is.na(data_europe[j,8:58]))<49) {
    #p-value of mann-kendall test
    data_europe$mk_pval[j] = mann_kendall_trend(dframe=data.frame(year=c(1960:2010),Qmax=as.numeric(data_europe[j,c(8:58)]) ) )[3] 
    #sen-slope on scale of data and as percentage of mean per decade
    sen_tmp = sen_slope(dframe=data.frame(year=c(1960:2010),Qmax=as.numeric(data_europe[j,c(8:58)]) ) )
    data_europe$sen_original[j] = sen_tmp[1]
    data_europe$pct_of_mean_per_decade[j] = sen_tmp[2]
  }
}

##############################################################################
#interpolating the sen-slopes via kriging
##############################################################################

library("sp")
library("rworldmap")
newmap <- getMap(resolution = "high")
europe_frame <- data.frame(country=NA,long=NA,lat=NA,plot_group=NA)
for(j in c(1:nrow(newmap@data))) {
  if(newmap@data$REGION[j]!="Europe"|is.na(newmap@data$REGION[j])) {next}
  else {
    n_polys = length(newmap@polygons[[j]]@Polygons)
    for(k in c(1:n_polys)) {
      tmp_frame <- data.frame(country=rep(newmap@polygons[[j]]@ID,nrow(newmap@polygons[[j]]@Polygons[[k]]@coords)),
                              long=newmap@polygons[[j]]@Polygons[[k]]@coords[,1],
                              lat=newmap@polygons[[j]]@Polygons[[k]]@coords[,2],
                              plot_group = rep(paste(newmap@polygons[[j]]@ID,k,sep="-"),nrow(newmap@polygons[[j]]@Polygons[[k]]@coords)))
      europe_frame <- rbind(europe_frame,tmp_frame)
    }
  }
}
europe_frame <- europe_frame[-1,]

#in order to perform kriging we need a different coordinate system than long-lat and a grid
spatial_krig = SpatialPointsDataFrame(coords=data_europe[,c("LON", "LAT")], data=data_europe[,c("pct_of_mean_per_decade","mk_pval")],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
spatial_krig_trafo = spTransform(spatial_krig,CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))

data_europe$lambert_x = spatial_krig_trafo@coords[,1]
data_europe$lambert_y = spatial_krig_trafo@coords[,2]

europe_krig = SpatialPointsDataFrame(coords=europe_frame[,c("long", "lat")], data=europe_frame[,c(1:4)],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
europe_krig_trafo = spTransform(europe_krig,CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))

europe_frame$lambert_x = europe_krig_trafo@coords[,1]
europe_frame$lambert_y = europe_krig_trafo@coords[,2]

#cover of europe for spatial grid
europe_cover = data.frame(x=c(-2700,-2700,2000,2000),y=c(-2500,2000,2000,-2500),data=rep("data",4))
coordinates(europe_cover) =~x+y
my_big_grid = makegrid(europe_cover,cellsize=20)
keep_point = rep(0,nrow(my_big_grid))
relevant_countries = unique(europe_frame$plot_group)

#for the grid, we check which points are contained in any country
for(country in relevant_countries)  {
  tmp_country = europe_frame[which(europe_frame$plot_group==country),]
  pt_in_country_tmp = point.in.polygon(point.x=my_big_grid[,1], point.y=my_big_grid[,2], pol.x=c(tmp_country[,"lambert_x"],tmp_country[1,"lambert_x"]), pol.y=c(tmp_country[,"lambert_y"],tmp_country[1,"lambert_y"]), mode.checked=FALSE)
  if(length(which(pt_in_country_tmp==1))>0) {
    keep_point[which(pt_in_country_tmp==1)] = 1
  }
  print(which(country==relevant_countries)/length(relevant_countries))
}

my_grid = my_big_grid[which(keep_point==1),]
names(my_grid) = c("lambert_x","lambert_y")

grid_krig = SpatialPointsDataFrame(coords=my_grid[,c(1,2)], data=my_grid[,c(1,2)],proj4string=CRS("+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))

#kriging
library("automap")
kriging_result <- autoKrige(pct_of_mean_per_decade~1, input_data=spatial_krig_trafo[which(data_europe$F1==1),], new_data = grid_krig)


##############################################################################
#estimating q100
##############################################################################
# Fig. 3 | Specific 100year floods ((m3/s)/km^2) in Europe (>= 30 years of data, n = 3738), 
#   where larger points indicate 90% confidence intervals smaller than 60% of the estimate.

# quantile function for the GEV distribution
qGEV <- function (p, xi=1, sigma=1, kappa=-.5) {
  xi + sigma * (1 - (-log(p))^kappa)/kappa
}

library(rstan)
# inference function for the stationary GEV model (3 parameters)
cod3 <- 'functions {
real gev_lpdf (real y, real xi, real sigma, real kappa) {
real z;
z = log(1 - kappa * (y - xi)/sigma)/(-kappa);
return -log(sigma) - (1 - kappa) * z - exp(-z);
}
}
data {
int<lower=0> N;
real y[N];
}
parameters {
real xi;
real log_sigma;
real kappa;
}
transformed parameters {
real sigma;
sigma = exp(log_sigma);
}
model {
// priors on component parameters
// xi ~ uniform
// log_sigma ~ uniform
kappa ~ normal(-0.1, 0.122);  // Martins and Stedinger (2000)
for(i in 1:N) {
y[i] ~ gev(xi, sigma, kappa);
}
}
generated quantities {
real log_lik[N];
for (i in 1:N) {
log_lik[i] = gev_lpdf(y[i] | xi, sigma, kappa);
}
}'
modelstGEV3 <- stan_model(model_code=cod3, model_name='modelstGEV3')


q_q100s <- matrix(NA, nrow=sum(data_europe$F3 > 0), ncol=7)
rownames(q_q100s) <- which(data_europe$F3>0)
colnames(q_q100s) <- c('q025','q050','q250','q500','q750','q950','q975')
for (ii in c(1:nrow(q_q100s)) ) {
  Q <- data_europe[as.numeric(rownames(q_q100s)[ii]),8:58]
  Q <- Q[!is.na(Q)]
  fit <- sampling(modelstGEV3,
                  data=list(N=length(Q), y=Q),
                  iter=10000, chains=4, thin=5)
  parfit <- extract(fit, pars=c('xi','sigma','kappa'))
  q100 <- qGEV(0.99, xi=parfit$xi, sigma=parfit$sigma, kappa=parfit$kappa)
  q_q100s[ii,] <- quantile(q100, prob=c(.025,.05,.25,.5,.75,.95,.975))
  print(ii/nrow(q_q100s))
}
med_q100 <- q_q100s[,'q500']  # posterior median estimate = color in Fig. 3
range_q100 <- q_q100s[,'q950'] - q_q100s[,'q050']
unc_q100 <- (range_q100/med_q100 > 0.6)  
# if TRUE the points in Fig. 3 are small
# if FALSE the points in Fig. 3 are big

data_europe$post_median_q100 <- NA
data_europe$ci_small <- NA
data_europe$post_median_q100[which(data_europe$F3==1)] = med_q100
data_europe$ci_small[which(data_europe$F3==1)] = ifelse(unc_q100,"small","big")

