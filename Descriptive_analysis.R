
#The code is presented only for the case of anchovy. For sardine, the code is the same.


library(sp) 
library(raster)
library(carData)
library(car)
library(lattice)
library(latticeExtra)
library(rasterVis)
library(maptools)
library(rgeos)
library(nlme)
library(dismo)
library(gstat)
library(spData)
library(sf)
library(spdep)
library(survival)
library(Formula)
library(ggplot2)
library(Hmisc)
library(leaflet)
library(GGally)
library(maptools)
library(corrplot)
library(sdmpredictors)
library(tidyverse)

library(rworldmap)
data(wrld_simpl)


#### PREPARE DATABASES ####
######################################################################

pelacus<-read.csv("PELACUSrawNASCgrid2013-19.csv")
str(pelacus)

##ANCHOVY DATASET --------------------

#prepare dataframe with presence/absence and coordinates for the anchovy
anchoa.dat<-as.data.frame(cbind(pelacus$anchoa, pelacus$x, pelacus$y))
colnames(anchoa.dat)<-c("anchoa","x","y")
View(anchoa.dat)
str(anchoa.dat)
dupl <- duplicated(anchoa.dat);dupl
sum(dupl)
anch_newd <-  anchoa.dat[!dupl, ]
str(anch_newd)



#Environmental Variables
library(rgdal)

bati<- raster("./envar_variables/batimetria_mean.tif")
clor<- raster("./envar_variables/clor_A_mean.tif")
sst<- raster("./envar_variables/sst_mean.tif")
sss<- raster("./envar_variables/sss_mean.tif")


#standardize
bati<-scale(bati)
sss<-scale(sss)
sst<-scale(sst)
clor<-scale(clor)


#verify that all of them have mean 0
round(summary(bati),2)
round(summary(sss),2)
round(summary(sst),2)
round(summary(clor),2)

#extract coordinates of variables in the studied area
coordenadas<-cbind(anch_newd$x, anch_newd$y)
colnames(coordenadas)<-c("x","y")

pres_bati <- extract(bati, coordenadas)
pres_clor <- extract(clor, coordenadas)
pres_log_clor <- extract(log_clor, coordenadas)
pres_sst <- extract(sst, coordenadas)
pres_sss <- extract(sss, coordenadas)

andata <- data.frame(cbind(anch_newd$anchoa,pres_bati,pres_clor,pres_log_clor,pres_sst,pres_sss, coordenadas))
colnames(andata)<- c("anchoa","batimetria","clorofila","logClorofila","temperatura","salinidad","longitud","latitud")
head(andata)
View(andata)


#Remove NAs
summary(andata)
to.remove <- which(!complete.cases(andata))
andata <- andata[-to.remove,]
summary(andata)


#Save the new dataframe as .csv
write.csv(andata, file="Engraulius_encrasicolus.csv", row.names = F )


##SARDINE DATASET --------------------

#same as the previous one. In this case, the dataframe is called  
write.csv(sardata, file="Sardina_pichardus.csv", row.names = F)



######################################################################
#### Descriptive statistics ####
##################################


#Step 1# CORRELATION between explanatory variables
#Highly correlated variables should be eliminated (=>0.7) 

matrix<-rcorr(as.matrix(E.encrasicolus[,c(1:6)]), type = "pearson")

# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

corrplot(matrix$r, type="lower", tl.col = "black",method="number",
         p.mat = matrix$P, sig.level = 0.05)



#Step 2# Check if the response variables and the explicative variables have linear relationship

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

g = ggpairs(E.encrasicolus[,c(1:6)], lower = list(continuous = my_fn))
g


#Step 3# COLINEALITY study 

#values higher than 3 should be consider collinear.

#to extract variance inflation factors (VIF) run their HighStatLib.r code and use the function corvif. 
source("HighstatLib.r")
corvif(E.encrasicolus[,c(1:6)])
corvif(E.encrasicolus[,c(2,4,5)])















