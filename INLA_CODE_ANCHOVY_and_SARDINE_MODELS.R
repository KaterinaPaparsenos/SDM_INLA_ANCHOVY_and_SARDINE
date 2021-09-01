
#Libraries
library(sp)
library(rgdal) 
library(carData)
library(car)
library(nlme)
library(gstat)
library(sf)
library(spData)
library(spdep)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(Hmisc)
library(raster) 
library(leaflet)
library(GGally)
library(maptools)
library(rgeos)
library(maptools) 
library(rJava)
library(Matrix)
library(parallel)
library(foreach)
library(INLA) 
library(dotCall64)
library(grid)
library(spam)
library(fields)

#dataset
E.encrasicolus<-read.csv("Engraulius_encrasicolus.csv")
str(E.encrasicolus)

####### CREATE MESH #############

##1. Load the countries shapefiles and join them ----------------------------------

#download use "download =T"
spain <- getData('GADM',country="ESP",level=0, download =T)
france<-getData("GADM", country="FRA", level=0, download =T) 
portugal<-getData("GADM", country="PRT", level=0, download =T)

#load 
spain <- readRDS("gadm36_ESP_0_sp.rds")
france<- readRDS("gadm36_FRA_0_sp.rds")
portugal<- readRDS("gadm36_PRT_0_sp.rds")

#join shapefiles
ab<-union(spain, france)
paises<-union(ab, portugal)

#cut out study area
ext<-extent(-13,0,41,46)
cat <- crop(paises, ext) 
plot(cat)

#save the area as a file .rds
cat<- readRDS("area_mesh.rds");
plot(cat)

##2. Define polygon containing data ----------------------------------

poligono.x<-c(min(E.encrasicolus$longitud) -1,
              min(E.encrasicolus$longitud) +1.7, 
              min(E.encrasicolus$longitud) +1.7,
              max(E.encrasicolus$longitud) +0.7, 
              max(E.encrasicolus$longitud) +0.7,
              min(E.encrasicolus$longitud) -1)

poligono.y<-c(min(E.encrasicolus$latitud)-1,
              min(E.encrasicolus$latitud)-1, 
              min(E.encrasicolus$latitud)+1,
              min(E.encrasicolus$latitud)+1,
              max(E.encrasicolus$latitud)+0.7,
              max(E.encrasicolus$latitud)+0.7)

xym<- as.matrix(data.frame(x = poligono.x,  
                           y = poligono.y))


p = Polygon(xym)
ps = Polygons(list(p),1) 
sps = SpatialPolygons(list(ps))

#cut the map
map_rec<-crop(cat, sps)

proj4string(sps)<-proj4string(map_rec)

# select polygon with the data
coast <- gDifference(sps, map_rec )
plot(coast)
points(E.encrasicolus[,6:7], pch=20)


##3. Define MESH y plot it ----------------------------------

boundary=inla.nonconvex.hull(as.matrix(E.encrasicolus[,7:8]))
mesh<-inla.mesh.2d(boundary=boundary, max.edge=c(0.4, 0.9),
                   cutoff=0.05,  offset=c(-0.1, -0.3))

#pdf("mesh.pdf", width = 10, height = 10)
plot(mesh)
plot(map_rec, add=TRUE, col="gray")
plot(mesh, add=TRUE)
points(E.encrasicolus[,7:8][which(E.encrasicolus$anchoa==0),], col="blue", pch=20)
points(E.encrasicolus[,7:8][which(E.encrasicolus$anchoa==1),], col="red", pch=20)
#dev.off()


####### MODEL PREPARATION #############

### -- Define SPDE -- ###
spde <- inla.spde2.matern(mesh)

### -- Matrix link data with MESH -- ###
A.est <- inla.spde.make.A(mesh, loc=cbind(E.encrasicolus$longitud, E.encrasicolus$latitud))

### -- inla.stack for estimations -- ###
stk.est<-inla.stack(data=list(y=E.encrasicolus$anchoa),
                    A=list(A.est, 1),
                    effects=list(spatial=1:spde$n.spde,
                                 data.frame(beta0=1, E.encrasicolus)),
                    tag='est')


#######  BEST MODEL SELECCION ####### 

source("Bdiclcpomodel_stack.R")

variables <- c("batimetria", "clorofila","temperatura","salinidad", "f(spatial, model=spde)")

resp=E.encrasicolus$anchoa

# call the function 
models_bin<-Bdiclcpomodel_stack(resp=resp, variables=variables, datos=inla.stack.data(stk.est),
                                n=20,family="binomial",control.predictor=list(compute=TRUE,A=inla.stack.A(stk.est)), 
                                control.compute = list(config=TRUE, dic=TRUE, cpo=TRUE,waic=TRUE),num.threads=3,
                                control.inla=list(strategy="gaussian"),verbose=FALSE)

saveRDS(models_bin, "best_model_anchoa.rds") 
models_bin<- readRDS("best_model_anchoa.rds")


models_bin$`Modelos dic`[1:10,]
models_bin$`Modelos waic`[1:10,]
models_bin$`Modelos lcpo`[1:10,]


###################################################################
#################### INLA MODELIZATION #############################
###################################################################


####### MODEL ESTIMATION #############

#fitted model

formula.1 <- y~-1 + beta0 + batimetria + clorofila + temperatura  +
  f(spatial, model=spde) + f(year, model="iid")

formula.2 <- y~-1 + beta0 + batimetria + clorofila + temperatura  +
  f(spatial, model=spde) + f(year, model="iid")

formula.3 <- y~-1 + beta0 + batimetria + clorofila + temperatura  +
  f(spatial, model=spde) + f(year, model="rw1")

formula.4 <- y~-1 + beta0 + batimetria + clorofila + temperatura  +
  f(spatial, model=spde) + f(year, model="rw2")

formula.5 <- y~-1 + beta0 + batimetria + clorofila + temperatura  +
  f(spatial, model=spde) + f(year, model="ar", order=2)

model.est <- inla(formula.1, 
                  data=inla.stack.data(stk.est), family="binomial" ,
                  control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE, return.marginals=TRUE), 
                  control.predictor=list(A=inla.stack.A(stk.est), compute=TRUE, 
                                         quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)),                  
                  #control.family=list(quantile=c(0.025)),
                  #control.inla=list(strategy ="laplace"),
                  num.threads = 3,
                  verbose=T)

saveRDS(model.est, "model_1_anchovy.rds")
model.est<- readRDS("model_1_anchoy.rds")
summary(model.est)



####### MODEL EVALUATION #############


#### plots of CPO and PIT values
plot(model.est) 


#### posterior distributions and probabilities ----------------------------------

## environmental variables ----------------------------------


1-inla.pmarginal(0, model.est$marginals.fixed$batimetria) 
1-inla.pmarginal(0, model.est$marginals.fixed$clorofila) 
1-inla.pmarginal(0, model.est$marginals.fixed$temperatura) 

par(mfrow=c(1,3))
plot(model.est$marginals.fixed$batimetria, type="l",main="Dist.posterior de la pendiente de batimetría")
abline(v=0, col="red", lwd=2)
plot(model.est$marginals.fixed$clorofila, type="l",main="Dist. posterior de la pendiente de Clorofila A")
abline(v=0, col="red", lwd=2)
plot(model.est$marginals.fixed$temperatura, type="l",main="Dist. posterior de la pendiente de la temperatura")
abline(v=0, col="red", lwd=2)
par(mfrow=c(1,1))


## spatial effect ----------------------------------


### --- Check the range is smaller than the offset --- ###
spde.result = inla.spde2.result(model.est, "spatial", spde, do.transform=TRUE)

range<-inla.emarginal(function(x) x, spde.result$marginals.range.nominal[[1]]) 

### --- Check the range is smaller than the offset --- ###
range < max(diff(range(E.encrasicolus[7])), diff(range(E.encrasicolus[,8])))*0.40 #Yes!!!
#[1] TRUE

### --- Plot --- ###
par(mfrow=c(1,1))
#par(mfrow=c(2,2), mar=c(3,3,1,0.5)) 
plot(spde.result$marginals.range.nominal[[1]], type='l')  #¿Qué era esto?
#marginals.range.nominal =	Marginal densities for range

### --- Distribucion a posteriori hyperparametros theta 1 y 2 --- ##

prec.marg1 <- model.est$marginals.hyperpar$`Theta1 for spatial`
prec.marg2 <- model.est$marginals.hyperpar$`Theta2 for spatial`
plot(1/prec.marg1[,1],prec.marg1[,2],xlim=c(0,2),type='l',xlab="Dist. a posteriori de la Varianza de Theta1", ylab="")
plot(1/prec.marg2[,1],prec.marg2[,2],xlim=c(0,2),type='l',xlab="Dist. a posteriori de la Varianza de Theta2", ylab="")

### --- Interpolate the posterior mean and sd --- ##

# plot in a grid m X m 

### --- Customize the grid to predict --- ###
par(mfrow=c(1,1))
bbox(coast)
(dxy <- apply(bbox(coast),1, diff))
(r <- dxy[1]/dxy[2])
m<-150
proj.grid.mat <- 
  inla.mesh.projector(mesh, 
                      xlim=bbox(coast)[1,],
                      ylim=bbox(coast)[2,] ,
                      dims=c(r, 1)*m)

plot(coast)
points(proj.grid.mat$lattice$loc, pch=20, cex=0.5)

### --- clean (set NA to the values outside boundary) --- ###
ov <- over(SpatialPoints(proj.grid.mat$lattice$loc, coast@proj4string),
           coast)

### --- check grid points inside the map --- ###
i.map <- is.na(ov)

### Plot the points where we will predict ###
par(mar=c(0,0,0,0))
plot(sps)
points(proj.grid.mat$lattice$loc[!i.map,], col="red", cex=0.2)
points(proj.grid.mat$lattice$loc[i.map,], col="blue", cex=0.2)

### --- consider only those inside map --- ###
proj.grid.mat$lattice$loc[i.map, ]

### --- Project the values of the mean and sd of the spatial effect --- ###
mean.g <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$mean)
sd.g <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$sd)
quantile_0.025 <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$`0.025quant`)
quantile_0.975 <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$`0.975quant`)

sd.g[i.map] <- mean.g[i.map] <- quantile_0.025[i.map] <- quantile_0.975[i.map] <- NA


#pdf("spatial_effect.pdf", width=10, height = 10) #dev.off()
### --- Spatial effect --- ###
par(mfrow=c(2,2))
par(mar=c(2,3,3,6))


### --- Posterior mean --- ###
plot(sps, col="gray", main="Spatial mean")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           mean.g, add=TRUE)
plot(map_rec, add=TRUE)

### --- Posterior sd --- ###
plot(sps, col="gray", main="Spatial sd")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           sd.g, add=TRUE)
plot(map_rec, add=TRUE)


### --- Posterior q0.025 --- ###
plot(sps, col="gray", main="Spatial q0.025")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           quantile_0.025, add=TRUE)
plot(map_rec, add=TRUE)

### --- Posterior q0.975 --- ###
plot(sps, col="gray", main="Spatial q0.975")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           quantile_0.975, add=TRUE)
plot(map_rec, add=TRUE)





#### Predictions with present data ----------------------------------

sst<- raster("./envar_variables/sst_mean.tif")
bati<- raster("./envar_variables/batimetria_mean.tif")
cla<- raster("./envar_variables/clor_A_mean.tif")


ext<-extent(-16.95833,4.0,25.0,55.0)
sst<-crop(sst,ext)
bati<-crop(bati,ext)
clor<-crop(clor,ext)

bati<-scale(bati)
clor<-scale(clor)
sst<-scale(sst)

predictors<-stack(bati,clor,sst)
names(predictors) <- c("bathymetry","clorophyll A","temperature")
plot(predictors)


### --- Matrix which link the mesh with coordinates to predict--- ###
A.pred <- inla.spde.make.A(mesh, loc=proj.grid.mat$lattice$loc[!i.map, ])

### --- Stack to predict --- ###
stk.pred <- inla.stack(data=list(y=NA),
                       A=list(A.pred, 1), 
                       effects=list(spatial=1:spde$n.spde,
                                    data.frame(beta0 = 1, 
                                               extract(predictors, 
                                                       proj.grid.mat$lattice$loc[!i.map, ]))),
                       tag='pred')

stk <- inla.stack(stk.est, stk.pred)


#### --- model --- ###
model.pred <- inla(formula.1, 
                   data=inla.stack.data(stk), family="binomial",
                   control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1), #link:link is a vector of
                   #length given by the size of the response variable with values 1 if the corresponding
                   #data is missing and NA otherwise
                   control.inla=list(strategy = "simplified.laplace"), # Strategy
                   control.mode=list(theta=model.est$mode$theta, restart=TRUE), #Mode 
                   control.results=list(return.marginals.random=FALSE,
                                        return.marginals.predictor=FALSE), # Avoid some marginals
                   num.threads = 3,
                   verbose=FALSE)

saveRDS(model.pred, "model_pred.rds")
model.pred<-readRDS("model_pred.rds")


### index for the prediction data
idx <- inla.stack.index(stk, 'pred')$data

summary(model.pred$summary.fitted.val$mean[idx])

### --- Organize probabilities into a matrix to visualize --- ###
prob.mean <- prob.sd <- prob.0.025<- prob.0.975 <- matrix(NA, proj.grid.mat$lattice$dims[1],
                                                          proj.grid.mat$lattice$dims[2])
prob.mean[!i.map] <- c(model.pred$summary.fitted.val$mean[idx])
prob.sd[!i.map] <- c(model.pred$summary.fitted.val$sd[idx])
prob.0.025[!i.map] <- c(model.pred$summary.fitted.val$`0.025quant`[idx])
prob.0.975[!i.map] <- c(model.pred$summary.fitted.val$`0.975quant`[idx])

#### --- plot --- ###

#pdf("predictive.pdf", width=10, height = 10)
### --- Spatial effect --- ###
par(mfrow=c(2,2))
par(mar=c(2,3,3,6))

### --- posterior predictive mean --- ###
plot(sps, col="gray", main="Predictive mean")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.mean, add=TRUE, zlim=c(0,1))

plot(map_rec, add=TRUE)

# points(data[,2:3][which(data$presence==1),], col="red", pch=20)
# points(data[,2:3][which(data$presence==0),], col="blue", pch=20)

### --- posterior predictive sd --- ###
plot(sps, col="gray", main="Predictive sd")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.sd, add=TRUE)
# 
# points(data[,2:3][which(data$presence==1),], col="red", pch=20)
# points(data[,2:3][which(data$presence==0),], col="blue", pch=20)

plot(map_rec, add=TRUE)


### --- posterior predictive q0.025 --- ###
plot(sps, col="gray", main="Predictive q0.025")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.0.025, add=TRUE, zlim=c(0,1))
plot(map_rec, add=TRUE)

### --- posterior predictive q0.975 --- ###
plot(sps, col="gray", main="Predictive q0.975")
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.0.975, add=TRUE, zlim=c(0,1))
plot(map_rec, add=TRUE)

#dev.off()












