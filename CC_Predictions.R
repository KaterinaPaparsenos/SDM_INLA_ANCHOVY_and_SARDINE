######################################
###### Climate Change forcasts ######
#####################################

library(sp)
library(raster)
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
library(Matrix)
library(parallel)
library(foreach)
library(INLA)
library(dotCall64)
library(grid)
library(spam)
library(fields)
library("rworldmap")
library(dismo)
data(wrld_simpl)


#datasets

E.encrasicolus<-read.csv("Engraulius_encrasicolus.csv")
S.pichardus<-read.csv("Sardina_pichardus.csv")


paises<- readRDS("area_mesh.rds")
ext<-extent(-13,0,41,46)
cat <- crop(paises, ext) 


#Define polygon with data 

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

map_rec<-crop(cat, sps) #cut map
proj4string(sps)<-proj4string(map_rec)

coast <- gDifference(sps, map_rec )

#Define MESH 
boundary=inla.nonconvex.hull(as.matrix(E.encrasicolus[,7:8]))
mesh<-inla.mesh.2d(boundary=boundary, max.edge=c(0.4, 0.9),
                   cutoff=0.05,  offset=c(-0.1, -0.3))

#Define SPDE 
spde <- inla.spde2.matern(mesh)

#Matriz link 
A.est <- inla.spde.make.A(mesh, loc=cbind(E.encrasicolus$longitud, E.encrasicolus$latitud))

#inla.stack  
stk.est<-inla.stack(data=list(y=E.encrasicolus$anchoa),
                    A=list(A.est, 1),
                    effects=list(spatial=1:spde$n.spde,
                                 data.frame(beta0=1, E.encrasicolus)),
                    tag='est')


formula.1 <- y~-1 + beta0 + bathymetry  + temperatura  + f(spatial, model=spde)
model.est1<- readRDS("model_1_anchoy.rds")
summary(model.est1)


# plot in a grid m X m 

### --- Customize the grid to predict --- ###
par(mfrow=c(1,1))
bbox(coast)
dxy <- apply(bbox(coast),1, diff)
r <- dxy[1]/dxy[2]
m<-150
proj.grid.mat <- 
  inla.mesh.projector(mesh, 
                      xlim=bbox(coast)[1,],
                      ylim=bbox(coast)[2,] ,
                      dims=c(r, 1)*m)

### --- clean (set NA to the values outside boundary) --- ###
ov <- over(SpatialPoints(proj.grid.mat$lattice$loc, coast@proj4string),
           coast)

### --- check grid points inside the map --- ###
i.map <- is.na(ov)

### --- consider only those inside map --- ###
proj.grid.mat$lattice$loc[i.map, ]


##### Variables: RPCs #####

### RPCs for sst
t50_45<- raster("./predic_variables/2050/2050.RCP45.Surface.Temperature.Mean..tif")
t50_85<- raster("./predic_variables/2050/2050.RCP85.Surface.Temperature.Mean.tif")
t100_45<- raster("./predic_variables/2100/2100.RCP45.Surface.Temperature.Mean.tif")
t100_85<- raster("./predic_variables/2100/2100.RCP85.Surface.Temperature.Mean.tif")

bati<- raster("./envar_variables/batimetria_mean.tif")

ext<-extent(-13,0,41,46)

t5.45<-crop(t50_45,ext)
t5.85<-crop(t50_85,ext)
t10.45<-crop(t100_45,ext)
t10.85<-crop(t100_85,ext)
bati<-crop(bati,ext)


#Standarize
bati<-scale(bati)
t5.45<-scale(t5.45)
t5.85<-scale(t5.85)
t10.45<-scale(t10.45)
t10.85<-scale(t10.85)


### Form dataframes for each RCP:

# RPC4.5 2050
predict_4.5_2050<-stack(bati,t5.45)
names(predict_4.5_2050) <- c("batthymetry","temp 2050 rcp 4.5")
plot(predict_4.5_2050)

# RPC8.5 2050
predict_8.5_2050<-stack(bati,t5.85)
names(predict_8.5_2050) <- c("bathymetry","temp 2050 rcp 8.5")
plot(predict_8.5_2050)

# RPC4.5 2100
predict_4.5_2100<-stack(bati,t10.45)
names(predict_4.5_2100) <- c("bathymetry","temp 2100 rcp 4.5")
plot(predict_4.5_2100)

# RPC8.5 2100
predict_8.5_2100<-stack(bati,t10.85)
names(predict_8.5_2100) <- c("bathymetry","temp 2100 rcp 8.5")
plot(predict_8.5_2100)



### -------- Forcast example -------- ####


##### PREDICCION 2050  RCP 4.5 #####


### --- Matrix which link the mesh with coordinates to predict--- ###
A.pred <- inla.spde.make.A(mesh, loc=proj.grid.mat$lattice$loc[!i.map, ])

### --- Stack to predict --- ###
stk.pred <- inla.stack(data=list(y=NA),
                       A=list(A.pred, 1), 
                       effects=list(spatial=1:spde$n.spde,
                                    data.frame(beta0 = 1, 
                                               extract(predict_4.5_2050, 
                                                       proj.grid.mat$lattice$loc[!i.map, ]))),
                       tag='pred')

stk <- inla.stack(stk.est, stk.pred)


#### --- model --- ###

model.pred.2050.4 <- inla(formula.1, 
   data=inla.stack.data(stk), family="binomial",
   control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1), #link:link is a vector of
#length given by the size of the response variable with values 1 if the corresponding
#data is missing and NA otherwise
 control.inla=list(strategy = "simplified.laplace"), # Strategy
 control.mode=list(theta=model.est1$mode$theta, restart=TRUE), #Mode 
   control.results=list(return.marginals.random=FALSE,
   return.marginals.predictor=FALSE), # Avoid some marginals
   num.threads = 4,
   verbose=FALSE)

saveRDS(model.pred.2050.4, "anchovy_pred2050_45.rds")

model.pred.2050.4<-readRDS("./anchovy_pred2050_45.rds")

### index for the prediction data
idx <- inla.stack.index(stk, 'pred')$data

summary(model.pred.2050.4$summary.fitted.val$mean[idx])



##### Plot predictions #####

### --- Organize probabilities into a matrix to visualize --- ###
prob.mean <- prob.sd <- prob.0.025<- prob.0.975 <- matrix(NA, proj.grid.mat$lattice$dims[1],
                                                          proj.grid.mat$lattice$dims[2])
prob.mean[!i.map] <- c(model.pred.2050.4$summary.fitted.val$mean[idx])
prob.sd[!i.map] <- c(model.pred.2050.4$summary.fitted.val$sd[idx])
prob.0.025[!i.map] <- c(model.pred.2050.4$summary.fitted.val$`0.025quant`[idx])
prob.0.975[!i.map] <- c(model.pred.2050.4$summary.fitted.val$`0.975quant`[idx])

#### --- plot --- ###

#pdf("predictive_2050_4.pdf", width=10, height = 10)
### --- Spatial effect --- ###
par(mfrow=c(2,2))
par(mar=c(2,3,3,6))

depth=raster("./envar_variables/batimetria_mean.tif")#read bathymtery map
depth=abs(depth) #change negative values to positive
depth=crop(depth, ext)

#trasfrom in a dataframe the bathymetry raster
matrix<- cbind(coordinates(depth), depth=getValues(depth))
I <- is.na(matrix[,3])
matrix<- matrix[!I,] #delete NA
matrix<-as.data.frame(matrix)
new<- subset(matrix, matrix$depth <=500)#crop at 350 meters
xy <- cbind(new$x, new$y)
rast<- raster(xmn=-10.74185, xmx=-0.590278, ymn=40.85342, ymx=45.58385 , nrows=200, ncols=431)#same of raster
p<- rasterize(xy, rast, new$depth, fun=max,na.rm=F)
p<-resample(p,depth)#predictors must have same resolution
e<-extent(p)



### --- posterior predictive mean --- ###
prob.mean.raster<-raster(list(x = proj.grid.mat$x, 
                              y = proj.grid.mat$y,
                              z = prob.mean ))
ext<-extent(-10, 0, 40, 46)
sp<-crop(prob.mean.raster,e) 
sp=resample(sp,p)
sp<-raster::mask(sp,p)

plot(sp,col=tim.colors(100)[1:100],main="Predictive mean", axes=T)
br.sp <- getMap(resolution = "high") 
plot(br.sp, col='dark grey',add=T)
contour(depth,col='dark grey',add=T)


#writeRaster(sp, filename="pres_mean_anchovy_2050_4.asc", format="ascii", overwrite=TRUE)

### --- posterior predictive sd --- ###
prob.sd.raster<-raster(list(x = proj.grid.mat$x, 
                            y = proj.grid.mat$y,
                            z = prob.sd ))

sp<-crop(prob.sd.raster,e) 
sp=resample(sp,p)
sp<-raster::mask(sp,p)

plot(sp,col=tim.colors(100)[1:100],main="Predictive sd", axes=T)
br.sp <- getMap(resolution = "high") 
plot(br.sp, col='dark grey',add=T)
contour(depth,col='dark grey',add=T)

### --- posterior predictive q0.025 --- ###
prob.q0.025.raster<-raster(list(x = proj.grid.mat$x, 
                                y = proj.grid.mat$y,
                                z = prob.0.025 ))

sp<-crop(prob.q0.025.raster,e) 
sp=resample(sp,p)
sp<-raster::mask(sp,p)

plot(sp,col=tim.colors(100)[1:100],main="Predictive q0.025", axes=T)
br.sp <- getMap(resolution = "high") 
plot(br.sp, col='dark grey',add=T)
contour(depth,col='dark grey',add=T)


### --- posterior predictive q0.975 --- ###
prob.q0.975.raster<-raster(list(x = proj.grid.mat$x, 
                                y = proj.grid.mat$y,
                                z = prob.0.975 ))

sp<-crop(prob.q0.975.raster,e) 
sp=resample(sp,p)
sp<-raster::mask(sp,p)

plot(sp,col=tim.colors(100)[1:100],main="Predictive q0.975", axes=T)
br.sp <- getMap(resolution = "high") 
plot(br.sp, col='dark grey',add=T)
contour(depth,col='dark grey',add=T)

#dev.off()



### -------- Difference in forecasts -------- ####

anchoa_pred_actual <- resample(anchoa_pred_actual,anchoa_2050_4)

DFanchoa_2050_4<- overlay(anchoa_2050_4, anchoa_pred_actual,
                          fun = function(r1, r2) { return( r1 - r2) })

DFanchoa_2050_8<- overlay(anchoa_2050_8, anchoa_pred_actual,
                          fun = function(r1, r2) { return( r1 - r2) })

DFanchoa_2100_4<- overlay(anchoa_2100_4, anchoa_pred_actual,
                          fun = function(r1, r2) { return( r1 - r2) })

DFanchoa_2100_8<- overlay(anchoa_2100_8, anchoa_pred_actual,
                          fun = function(r1, r2) { return( r1 - r2) })

DFanchoa_4.5<- overlay(anchoa_2100_4, anchoa_2050_4,
                       fun = function(r1, r2) { return( r1 - r2) })

DFanchoa_8.5<- overlay(anchoa_2100_8, anchoa_2050_8,
                       fun = function(r1, r2) { return( r1 - r2) })


par(mfrow=c(2,3))
par(adj=0.5)
par(mar=c(3, 4, 1.5, 4), mgp=c(3, 1, 0), las=1)


plot(DFanchoa_2050_4,col=tim.colors(100)[1:100],main="RCP4.5 2050", axes=T)
plot(crop(br.sp, ext), col='dark grey',add=T)
#contour(depth,col='dark grey', add=T)

plot(DFanchoa_2100_4,col=tim.colors(100)[1:100],main="RCP4.5 2100", axes=T)
plot(crop(br.sp, ext), col='dark grey',add=T)
#contour(depth,col='dark grey',add=T)

plot(DFanchoa_4.5,col=tim.colors(100)[1:100],main="Difference RCP4.5 2050-2100", axes=T)
plot(crop(br.sp, ext), col='dark grey',add=T)
#contour(depth,col='dark grey',add=T)


plot(DFanchoa_2050_8,col=tim.colors(100)[1:100],main="RCP8.5 2050", axes=T)
plot(crop(br.sp, ext), col='dark grey',add=T)
#contour(depth,col='dark grey',add=T)
#plot(br.sp, col='dark grey') #,add=T

plot(DFanchoa_2100_8,col=tim.colors(100)[1:100],main="RCP8.5 2100", axes=T)
plot(crop(br.sp, ext), col='dark grey',add=T)
#contour(depth,col='dark grey',add=T)

plot(DFanchoa_8.5,col=tim.colors(100)[1:100],main="Difference RCP8.5 2050-2100", axes=T)
plot(crop(br.sp, ext), col='dark grey',add=T)
#contour(depth,col='dark grey',add=T)
















