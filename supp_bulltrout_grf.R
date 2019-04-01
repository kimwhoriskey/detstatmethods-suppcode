# Supporting information for Whoriskey et al. Statistical methods for detection data.

# Analyzing the bull trout dataset with a Gaussian Random Field. 
# Author: Marie Auger-Méthé & Kim Whoriskey


library(plyr)
library(sp)
library(rgdal)
library(rgeos)
library(INLA)
library(TMB)
library(lattice)


########
# Read data
# Receivers data
recDD <- read.csv("Kinbasket_Receiver_groups.csv", stringsAsFactors = FALSE)

# Fish data
fishD <- read.csv("Jan_2011_BT.csv", stringsAsFactors = FALSE)

# Simple analysis summary for January - and only counting each individual once
# Only need the first columns
fishDs <- unique(fishD[,1:6])
fishDsCount <- count(fishDs,c("Receiver","x_plane", "y_plane"))

# Add the receivers with no detections
noDet <- recDD[is.na(match(recDD$receiver_sn, fishDsCount$Receiver)), c("receiver_sn", "x_plane", "y_plane")]
noDet$freq <- 0
colnames(noDet) <- colnames(fishDsCount)
fishDsCount <- rbind(fishDsCount,noDet)
rownames(fishDsCount) <- 1:nrow(fishDsCount)


#########
# Read in and dissolve the lake shapefile
uunit <- readOGR(dsn="shpfile/.", layer='kinbasket')
uunitM <- gUnaryUnion(uunit)
plot(uunitM)
with(fishDsCount,points(x_plane,y_plane,cex=freq/10, col=rgb(1,0,0,0.4), pch=19))
# Check overlapping - one receiver not in the lake
fishDsCountsp <- fishDsCount
coordinates(fishDsCountsp) <- ~x_plane+y_plane
proj4string(fishDsCountsp) <- proj4string(uunitM)
inLake <- gIntersects(fishDsCountsp,uunitM, byid=TRUE)
# Remove the receiver not in the Lake
fishDsCount[!inLake,]
fishDsCount <- fishDsCount[inLake,]


#########
# Calculate the distance from the dam to each receiver
coordScale <- 1
coords <- fishDsCount[,c("x_plane","y_plane")]/coordScale
# according to wikipedia the location of the dam is: 52°04′40″N 118°33′59″W
ptDam <- as.data.frame(cbind(as.numeric(char2dms("118d33'59\"W")),as.numeric(char2dms("52d04'40\"N"))))
coordinates(ptDam) <- ptDam
proj4string(ptDam) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
ptDamC <- spTransform(ptDam, proj4string(uunit))
ptInt <- coordinates(ptDamC)
distFromPT <- sqrt((coords[,1] -ptInt[1])^2 + (coords[,2] -ptInt[2])^2)/1000 # distance from point of interest, in km
data <- data.frame(y=fishDsCount$freq, coords=coords, distPT=distFromPT)


#########
# INLA mesh
# Calculate inla mesh
bound1 <- inla.sp2segment(uunitM)
mesh <- inla.mesh.2d(coords, bound=bound1, max.edge=1000, cutoff=100)

# Plot the mesh
sizingCex <- 10
plot(mesh, main='',asp=1, col="aliceblue", edge.color="steelblue", draw.segments=FALSE)
with(fishDsCount,points(x_plane/coordScale,y_plane/coordScale,cex=freq/sizingCex, col=rgb(1,0,0,0.4), pch=19))
lines(c(392000, 392656.3), c(5771600.5, 5770846.5), lwd=10, col="goldenrod", pch=20, cex=5) #dam
with(recDD,points(x_plane,y_plane,pch=19, cex=0.1))

# Caclulate estimated expected surface - INLA
mesh_proj <- inla.mesh.projector(mesh, xlim=range(coords[,1]), ylim=range(coords[,2]), 
                                 dims=c(100,100))




#########
# Fit the GRF with TMB 

# Fixed effects 
X <- model.matrix(~1 + distPT, data=data)
data_tmb <- list(counts=data$y, meshidxloc=mesh$idx$loc-1, X=as.matrix(X))

# SPDE part
spde_param <- inla.spde2.matern(mesh, alpha=2)$param.inla # alpha =2 -> nu = 1
# need to get all SparseMatrix elements of the SPDE object - to calculate the
# inverse-covariance matrix of GRF c("M0","M1","M2") is G0, G1, G2 in eqn (10)
# in Lindgren 2011
data_tmb$spde <- spde_param[c("M0","M1","M2")]
n_s <- nrow(data_tmb$spde$M0)


# Compile the grf likelihood template
compile("supp_grf.cpp")
dyn.load(dynlib("supp_grf"))

# Starting values for the parameters 
init_par <- list(beta=c(0,0),log_tau=0,log_kappa=0,x=rep(0.0,n_s))

# Construct likelihood object
obj <- MakeADFun(data_tmb,init_par,random="x",DLL="supp_grf")

# Optimize the likelihood function
opt <- nlminb(obj$par,obj$fn,obj$gr)

# Look at the convergence
opt$convergence
opt$message

# Look at the parameters, and calculate their standard errors and confidence intervals
opt$par
parGRF <- summary(sdreport(obj))
pparGRF <- data.frame(parGRF[rownames(parGRF) != "x",])
pparGRF$par <- rownames(parGRF[rownames(parGRF) != "x",])
pparGRF$LCI <- pparGRF[,1] + qnorm(0.025)*pparGRF[,2]
pparGRF$UCI <- pparGRF[,1] + qnorm(0.975)*pparGRF[,2]
pparGRF


pl <- obj$env$parList()
xEst <- pl$x/parGRF["tau",1]

#distance to the dam
distG <- sqrt((mesh$loc[,1] -ptInt[1])^2 + (mesh$loc[,2] -ptInt[2])^2)/1000 # distance from point of interest



# The predicted values (eta)
etaGRF <- inla.mesh.project(mesh_proj, xEst+opt$par[1]+opt$par[2]*distG)
etaP <- levelplot(etaGRF, xlab='', ylab='', main=substitute(eta==log(lambda)), 
                  col.regions=rev(topo.colors(99)), 
                  at=seq(min(etaGRF, na.rm=TRUE),max(etaGRF, na.rm=TRUE), length=100),
                  scales=list(draw=FALSE))
etaP

# The exp() of the predicted values, i.e., the predicted fish counts (lambda)
expEtaGRF <- inla.mesh.project(mesh_proj, exp(xEst+opt$par[1]+opt$par[2]*distG))
lambdaP <- levelplot(expEtaGRF, xlab='', ylab='', main=substitute(lambda), 
                     col.regions=rev(topo.colors(99)), 
                     at=seq(min(expEtaGRF, na.rm=TRUE),max(expEtaGRF, na.rm=TRUE), length=100),
                     scales=list(draw=FALSE))
lambdaP

# The effect of the dam (beta_1)
beta1d <- inla.mesh.project(mesh_proj, opt$par[2]*distG)
beta1dP <- levelplot(beta1d, xlab='', ylab='', main=substitute(beta[1]*dam), 
                     col.regions=rev(topo.colors(99)), 
                     at=seq(min(beta1d, na.rm=TRUE),max(beta1d, na.rm=TRUE), length=100),
                     scales=list(draw=FALSE))
beta1dP

# The random field values (xi)
GRFest <- inla.mesh.project(mesh_proj, xEst)
xiP <- levelplot(GRFest, xlab='', ylab='', main=substitute(xi), 
                  col.regions=rev(topo.colors(99)), at= seq(min(GRFest, na.rm=TRUE),max(GRFest, na.rm=TRUE), length=100), 
                  scales=list(draw=FALSE))
xiP

