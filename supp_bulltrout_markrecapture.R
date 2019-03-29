# Supporting information for Whoriskey et al. Statistical methods for detection data.

# Analyzing the bull trout dataset with mark-recapture analysis. 
# Author: Kim Whoriskey

######## 
# packages
require(RMark)
require(R2ucare)

########
# load the detection data
dets <- read.csv("Jan_2011_BT.csv")
# get the unique fish ids
ids <- unique(dets$FISHID)
# add a posix time column
dets$date <- paste(dets$Date, dets$Time.Hr)
dets$date <- as.POSIXct(as.character(dets$date), "%d/%m/%Y %H", tz="GMT")


########
# need to get the proper response variable for the mark-recapture analysis
# summarize the data as presence every week

# sequence of all the days
firstday <- min(dets$date)
lastday <- max(dets$date) + 86400
days <- seq(firstday, lastday, by=86400*7) #summarize every 7 days


# accumulators
ch <- numeric()     # the response variable
ch_matx <- matrix(0, nrow=88, ncol=(length(days)-1))
# only take the first four weeks for consistency
id <- numeric()     # id of the fish
sex <- numeric()    # sex of the fish, males are 0, females are 1
len <- numeric()    # length of the fish

for(i in 1:length(ids)){
  # ith individual
  sub <- dets[(dets$FISHID==ids[i]),]
  # create an encounter matrix
  for(j in 1:dim(ch_matx)[2]){
    # jth date interval
    sub2 <- sub[(sub$date >= days[j]) & (sub$date < days[j+1]),]
    if(dim(sub2)[1]>0) ch_matx[i,j] <- 1
  }
  
  # append everything to the vectors
  ch <- append(ch, paste(ch_matx[i,], collapse=""))
  id <- append(id, sub$FISHID[1])
  sex <- append(sex, ifelse(sub$Sex[1]=="m", 0, 1))
  len <- append(len, sub$Total.Length[1])
  
}

# make the data frame
dat <- data.frame(ch, id, sex, len)
dat$ch <- as.character(dat$ch)

# the 88th fish was not detected in the first four weeks
# (but was detected on the 31). Omit. 
dat <- dat[-88,]

# take a look
head(dat)


########
# fit the mark-recapture models in RMark

# set up all of the possible model formulas
int <- list(formula=~1) #intercept only model 
time <- list(formula=~time) #temporally varying model
sexlen <- list(formula=~sex+len) #sex and length model, additive effects
sexxlen <- list(formula=~sex*len) #sex and length model, interactive effects
sex <- list(formula=~sex) #just sex
len <- list(formula=~len) #just length


# fit each of the models
trout.simple <- mark(dat, model="CJS", model.parameters = list(Phi=int, p=int))
trout.ptime <- mark(dat, model="CJS", model.parameters = list(Phi=int, p=time))
trout.phisex <- mark(dat, model="CJS", model.parameters = list(Phi=sex, p=int), groups="sex")
trout.philen <- mark(dat, model="CJS", model.parameters = list(Phi=len, p=int))
trout.phisexlen <- mark(dat, model="CJS", model.parameters = list(Phi=sexlen, p=int), groups="sex")
trout.phisexxlen <- mark(dat, model="CJS", model.parameters = list(Phi=sexxlen, p=int), groups="sex")
trout.phisex.ptime <- mark(dat, model="CJS", model.parameters = list(Phi=sex, p=time), groups="sex")
trout.philen.ptime <- mark(dat, model="CJS", model.parameters = list(Phi=len, p=time))
trout.phisexlen.ptime <- mark(dat, model="CJS", model.parameters = list(Phi=sexlen, p=time), groups="sex")
trout.phisexxlen.ptime <- mark(dat, model="CJS", model.parameters = list(Phi=sexxlen, p=time),groups="sex")

# all model results
bulltrout.cjs.results <- collect.models()
bulltrout.cjs.results


# the results from the best model, according to AIC
trout.simple$results


# check assumptions using R2ucare
# Gimenez et al. 2018. https://doi.org/10.1111/2041-210X.13014
overall_CJS(ch_matx, freq=rep(1, dim(ch_matx)[1])) #not significant




