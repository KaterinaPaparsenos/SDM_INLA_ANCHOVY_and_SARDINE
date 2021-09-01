
  #### ENVIRONMENTAL VARIABLES PREPARATION ####
#################################################

#All the environmental variables have been downloaded from official online sites.  

#SEA SUFACE TEMPERATURE (SST) and SEA SURFACE SALINITY (SSS) have been obtained by
# Copernicus Marine resources: https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=IBI_MULTIYEAR_PHY_005_002

#Clorofila A (cla) have been obtained by 
#Ocean Color: https://www.oceancolour.org/browser/index.php?product=chlor_a&page=0&period=monthly&version=5&limit=24&from=2013-03-01&to=2019-06-01&mode=browse#%5B


##Libraries  
library(ncdf4)
library(sp)
library(raster)
library(rgdal)
library(scales)
library(ggplot2)
library(reshape)
library(dplyr)
library(lattice) 
library(RColorBrewer)


### Open netCDF file -------------------------------

  nc_file="variablefile.nc" #variable file netCDF
nc_data <- nc_open(nc_file)
print(nc_data) #query file metadata

lat_variable = 'latitude'
lon_variable = 'longitude'
nc_variable = '' #variable name sss or sst

variable = ncvar_get(nc_data,nc_variable) #extraccion valores variable
lats = ncvar_get(nc_data,lat_variable) 
lons = ncvar_get(nc_data,lon_variable)

dims_variable = dim(variable) #variable dimensions
mat = matrix(variable, nrow=dims_variable[1],ncol=dims_variable[2])

### levelplot
grid <- expand.grid(lon=lons, lat=lats)
cutpts <- c(-10,0,10,20,30,40)
levelplot(sst_slice ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
          col.regions=(rev(brewer.pal(10,"RdBu"))))




### SST ------------------------------------------

print(nc_data)

var_names <- c('thetao')

for (var_name in var_names) {
  
  # Create raster stack
  sst <- stack(
    raster('sst_3_2013.nc', varname = var_name),
    raster('sst_4_2013.nc', varname = var_name),
    raster('sst_5_2013.nc', varname = var_name),
    raster('sst_3_2014.nc', varname = var_name),
    raster('sst_4_2014.nc', varname = var_name),
    raster('sst_5_2014.nc', varname = var_name),
    raster('sst_3_2015.nc', varname = var_name),
    raster('sst_4_2015.nc', varname = var_name),
    raster('sst_5_2015.nc', varname = var_name),
    raster('sst_3_2016.nc', varname = var_name),
    raster('sst_4_2016.nc', varname = var_name),
    raster('sst_5_2016.nc', varname = var_name),
    raster('sst_3_2017.nc', varname = var_name),
    raster('sst_4_2017.nc', varname = var_name),
    raster('sst_5_2017.nc', varname = var_name),
    raster('sst_3_2018.nc', varname = var_name),
    raster('sst_4_2018.nc', varname = var_name),
    raster('sst_5_2018.nc', varname = var_name),
    raster('sst_3_2019.nc', varname = var_name),
    raster('sst_4_2019.nc', varname = var_name),
    raster('sst_5_2019.nc', varname = var_name))
  
  # Name each layer
  names(sst) <- c("03_2013", "04_2013", "05_2013",
                  "03_2014", "04_2014", "05_2014",
                  "03_2015", "04_2015", "05_2015",
                  "03_2016", "04_2016", "05_2016",
                  "03_2017", "04_2017", "05_2017",
                  "03_2018", "04_2018", "05_2018",
                  "03_2019", "04_2019", "05_2019" ) 
  
  writeRaster(x = sst, 
              filename = paste0(var_name, '_variable.nc'),
              overwrite = TRUE, 
              format = 'CDF')
}

str(sst)
dim(sst)
nlayers(sst)
crs(sst) 
res(sst) 
yres(sst) ; xres(sst) 


#Once we have created our RasterStack, we can visualize them. 
#We can use the ggplot() command to create a multi-panelled plot showing each band in our RasterStack.
#First we need to create a data frame object.
#Because there are multiple bands in our data, 
#we will reshape (or “melt”) the data so that we have a single column with the NDVI observations. 
#We will use the function melt() from the reshape package to do this:

# Grafics
sst_stack_df <- as.data.frame(sst, xy = TRUE) %>%
  melt(id.vars = c('x','y'))

sst_mes<-ggplot() +
  geom_raster(data = sst_stack_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable)+
  xlab("longitud") + ylab("latitude")+
  ggtitle("Temperatura superficial del mar 2013-2019") +
  theme(plot.title = element_text(lineheight=.8, face="bold",size = 20)) +
  theme(text = element_text(size=10))

sst_mes


hist_sst<-ggplot(sst_stack_df) +
  geom_histogram(aes(value)) +
  facet_wrap(~variable)

hist_sst


### SSS ------------------------------------------

print(nc_data)

var_names <- c('so')

for (var_name in var_names) {
  
  # Create raster stack
  sss <- stack(
    raster('sss_3_2013.nc', varname = var_name),
    raster('sss_4_2013.nc', varname = var_name),
    raster('sss_5_2013.nc', varname = var_name),
    raster('sss_3_2014.nc', varname = var_name),
    raster('sss_4_2014.nc', varname = var_name),
    raster('sss_5_2014.nc', varname = var_name),
    raster('sss_3_2015.nc', varname = var_name),
    raster('sss_4_2015.nc', varname = var_name),
    raster('sss_5_2015.nc', varname = var_name),
    raster('sss_3_2016.nc', varname = var_name),
    raster('sss_4_2016.nc', varname = var_name),
    raster('sss_5_2016.nc', varname = var_name),
    raster('sss_3_2017.nc', varname = var_name),
    raster('sss_4_2017.nc', varname = var_name),
    raster('sss_5_2017.nc', varname = var_name),
    raster('sss_3_2018.nc', varname = var_name),
    raster('sss_4_2018.nc', varname = var_name),
    raster('sss_5_2018.nc', varname = var_name),
    raster('sss_3_2019.nc', varname = var_name),
    raster('sss_4_2019.nc', varname = var_name),
    raster('sss_5_2019.nc', varname = var_name))
  
  # Name each layer
  names(sss) <- c("03_2013", "04_2013", "05_2013",
                  "03_2014", "04_2014", "05_2014",
                  "03_2015", "04_2015", "05_2015",
                  "03_2016", "04_2016", "05_2016",
                  "03_2017", "04_2017", "05_2017",
                  "03_2018", "04_2018", "05_2018",
                  "03_2019", "04_2019", "05_2019" ) 
  
  writeRaster(x = sss, 
              filename = paste0(var_name, '_variable.nc'),
              overwrite = TRUE, 
              format = 'CDF')
}

dim(sss)
plot(sss)
str(sss)
nlayers(sss)
crs(sss) 
res(sss) 


#Grafics

sss_stack_df <- as.data.frame(sss, xy = TRUE) %>%
  melt(id.vars = c('x','y'))

sss_mes<-ggplot() +
  geom_raster(data = sss_stack_df , aes(x = x, y = y, fill = value)) +
  facet_wrap(~ variable)+
  xlab("longitud") + ylab("latitude")+
  ggtitle("Salinidad del mar 2013-2019") +
  theme(plot.title = element_text(lineheight=.8, face="bold",size = 20)) +
  theme(text = element_text(size=10))

sss_mes

hist_sss<-ggplot(sss_stack_df) +
  geom_histogram(aes(value)) +
  facet_wrap(~variable)

hist_sss



### ClA ------------------------------------------


setwd("./clorofila a")

nc_file="./clorofila a/CHLOR_A_2014_3.nc"
nc_data <- nc_open(nc_file)

print(nc_data)
variable<-(list.files("./clorofila a", full.names=T, pattern=".nc"))#change directory


cla <- stack(variable, varname= "chlor_a")

var_name <- c("chlor_a")

names(cla) <- c("03_2013", "04_2013", "05_2013",
                "03_2014", "04_2014", "05_2014",
                "03_2015", "04_2015", "05_2015",
                "03_2016", "04_2016", "05_2016",
                "03_2017", "04_2017", "05_2017",
                "03_2018", "04_2018", "05_2018",
                "03_2019", "04_2019", "05_2019"  )

plot(cla)

#cut the area of the study.
ext<-extent(-30,10,20,60)
cla<-crop(cla,ext)
plot(cla)
extend(cla, ext) #information of the cut made

#cuidado con la resolución que no coincide con las demás variables
#hay que cambiarlo

#Careful with the resolution, must match with the resolution of the other variables
#need to be changed
res(cla)
#aggregate from 0.04166667x0.04166667 resolution to 120x120 (factor = 2)
cla <- aggregate(cla, fact=2)
res(cla)
#[1] 0.08333333 0.08333333



writeRaster(x = cla, 
            filename = paste0(var_name, '_variable.nc'),
            overwrite = TRUE, 
            format = 'CDF')


str(cla)
dim(cla)
nlayers(cla)
dim(cla)
crs(cla) 
res(cla) 


## Grafics

cla_stack_df <- as.data.frame(cla, xy = TRUE) %>% melt(id.vars = c('x','y'))


cla_mes<-ggplot() +
  geom_raster(data = cla_stack_df  , aes(x = x, y = y, fill = values)) +
  facet_wrap(~ variable)+
  xlab("longitud") + ylab("latitude")+
  ggtitle("Concentración clorofila a") +
  theme(plot.title = element_text(lineheight=.8, face="bold",size = 20)) +
  theme(text = element_text(size=10))


cla_mes ##esta grafica no me gusta, hay que cambiarla


hist_cla<-ggplot(sst_stack_df) +
  geom_histogram(aes(value)) +
  facet_wrap(~variable)

hist_cla




### Bathymetry ------------------------------------------

library("sdmpredictors")
library("leaflet")
bathy<- load_layers(c("MS_bathy_5m"))
summary(bathy)
summary(bathy, maxsamp = ncell(bathy))
plot(bathy)
ext<-extent(-30,10,20,60)
bathy<-crop(bathy,ext)
plot(bathy)
extend(bathy, ext)


writeRaster(x = bathy, 
            filename = "batimetria_mean",
            overwrite = TRUE, 
            format = 'GTiff',
            NAflag=-9999) #CDF files requier library(ncdf4)); GTiff files requiere library(rgdal)


GDALinfo("batimetria_mean.tif")#displaying file attributes
bati<- raster("batimetria_mean.tif")
summary(bati)
summary(bati, maxsamp = ncell(bati))
plot(bati)

#same image with ggplot2
info_df <- raster::as.data.frame(bati, xy = TRUE)
str(info_df)
ggplot() +
  geom_raster(data = info_df , aes(x = x, y = y, fill = batimetria_mean)) +
  coord_equal() +
  scale_fill_gradientn(colours = terrain.colors(10))+
  ggtitle("Batimetría") +
  theme(plot.title = element_text(lineheight=.8, face="bold",size = 20)) +
  theme(text = element_text(size=10))




###########################################################
##--- OBTAIN THE MEAN VALUE OF EACH VARIABLE ---##
###########################################################

###1. SST ----------------------------------------------

nc_file="./envar_variables/thetao_variable.nc"
nc_data <- nc_open(nc_file)
print(nc_data) 

#in this netCFD file my variable (sst) is called "variable" and the time is shown as "z".
#z=21 (each month we saved)

##Name of variables I am interested in: 
lat_variable = 'latitude'
lon_variable = 'longitude'
time_variable ='z'
nc_variable = 'variable'

variable = ncvar_get(nc_data,nc_variable) #extraction variable values
lats = ncvar_get(nc_data,lat_variable) #get latitude
lons = ncvar_get(nc_data,lon_variable) #get longitude 
nlon<-dim(lons)
nlat<- dim(lats)
print(c(nlon,nlat))
times = ncvar_get(nc_data,time_variable) #get time (months)
nt<- dim(times)
dims_variable = dim(variable) #variable dimensions

sst_mat = matrix(variable, nrow=dims_variable[1],ncol=dims_variable[2]) 

#create a matrix with the variable values according to longitude and latitude in order to represent them.


#Data stored in netCDF files can have different types:
#2-D raster slabs (example: longitude and latitude layers).
#3_D bricks (example long-lat-time)
#4_D arrays (longitude-latitude_height-time)

#However, in R we usually work with 2-D "tidy" data frame variables.
#Therefore, before performing the analysis, the 3-D and 4_D arrays in netCDF files 
#have to be "flattened" and converted into 2-D arrays.

#Convert the whole array to data frame and calculate mean

# reshape the array into vector
variable_long <- as.vector(variable)
length(variable_long)

# reshape the vector into a matrix
tmp_mat <- matrix(variable_long, nrow=nlon*nlat, ncol=nt)
dim(tmp_mat)
head(na.omit(tmp_mat))

# create dataframe
lonlat <- as.matrix(expand.grid(lons,lats))
tmp_df02 <- data.frame(cbind(lonlat,tmp_mat))
names(tmp_df02) <- c("lon","lat","tmpMar2013","tmpAbr2013","tmpMay2013",
                     "tmpMar2014","tmpAbr2014","tmpMay2014",
                     "tmpMar2015","tmpAbr2015","tmpMay2015",
                     "tmpMar2016","tmpAbr2016","tmpMay2016",
                     "tmpMar2017","tmpAbr2017","tmpMay2017",
                     "tmpMar2018","tmpAbr2018","tmpMay2018",
                     "tmpMar2019","tmpAbr2019","tmpMay2019")
# options(width=96)
head(na.omit(tmp_df02, 20))
tail(na.omit(tmp_df02, 10))
dim(tmp_df02)

# get the annual mean
tmp_df02$mean <- apply(tmp_df02[3:23],1,mean) # calculates the average per row (per pixel in the maps)
head(na.omit(tmp_df02$mean))
dim(na.omit(tmp_df02))

#Plot of the mean value (levelplot)
grid <- expand.grid(lon=lons, lat=lats)
cutpts <- c(0,3,6,9,12,15,18,21)
levelplot( tmp_df02$mean~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
           col.regions=(rev(brewer.pal(10,"RdBu"))))

# write out the dataframe as a .csv file
#csvpath <- "./var_ambientales"
#csvname <- "sst_mean.csv"
#csvfile <- paste(csvpath, csvname, sep="")
#write.table(na.omit(tmp_df02),csvfile, row.names=FALSE, sep=",")

#CREATE RASTER WITH LATITUDE, LONGITUD AND MEAN INFORMATION AND PLOT IT:
# from x,y,z-matrix
sst<- data.frame(cbind(tmp_df02$lon,tmp_df02$lat,tmp_df02$mean))
names(sst)<-c("lon","lat","sst_mean")
View(sst)

sst<- rasterFromXYZ(sst) #Convert first two columns as lon-lat and third as value  
plot(sst)
sst

#save the raster with the mean as a separate file
writeRaster(x = sst, 
            filename = "sst_mean",
            overwrite = TRUE, 
            format = 'GTiff',
            NAflag=-9999) 

#check if the raster with the mean values we have created is correct
library(rgdal)
GDALinfo("sst_mean.tif")
SST<- raster("sst_mean.tif")
summary(SST)
summary(SST, maxsamp = ncell(SST))
plot(SST)

#The same plot with ggplot2
info_df <- raster::as.data.frame(SST, xy = TRUE)
str(info_df)
ggplot() +
  geom_raster(data = info_df , aes(x = x, y = y, fill = sst_mean)) +
  coord_equal() +
  scale_fill_gradientn(colours = terrain.colors(10))+
  ggtitle("Temperatura media") +
  theme(plot.title = element_text(lineheight=.8, face="bold",size = 20)) +
  theme(text = element_text(size=10))


##2. SSS--------------------------------------------------------------------------------------

nc_file="./envar_variables/so_variable.nc"
nc_data <- nc_open(nc_file)
print(nc_data) 

lat_variable = 'latitude'
lon_variable = 'longitude'
time_variable ='z'
nc_variable = 'variable'

variable = ncvar_get(nc_data,nc_variable) 
lats = ncvar_get(nc_data,lat_variable) 
lons = ncvar_get(nc_data,lon_variable) 
nlon<-dim(lons)
nlat<- dim(lats)
print(c(nlon,nlat))
times = ncvar_get(nc_data,time_variable) 
nt<- dim(times)
#dims_variable = dim(variable) 

#sss_mat = matrix(variable, nrow=dims_variable[1],ncol=dims_variable[2]) 

# reshape the array into vector
variable_long <- as.vector(variable)
length(variable_long)

# reshape the vector into a matrix
salinity_mat <- matrix(variable_long, nrow=nlon*nlat, ncol=nt)
dim(salinity_mat)
head(na.omit(salinity_mat))

# create dataframe
lonlat <- as.matrix(expand.grid(lons,lats))
sal_df <- data.frame(cbind(lonlat,salinity_mat))
names(sal_df) <- c("lon","lat","salMar2013","salAbr2013","salMay2013",
                   "salMar2014","salAbr2014","salMay2014",
                   "salMar2015","salAbr2015","salMay2015",
                   "salMar2016","salAbr2016","salMay2016",
                   "salMar2017","salAbr2017","salMay2017",
                   "salMar2018","salAbr2018","salMay2018",
                   "salMar2019","salAbr2019","salMay2019")
# options(width=96)
head(na.omit(sal_df, 10))
dim(sal_df)

# get the annual mean
sal_df$mean <- apply(sal_df[3:23],1,mean) 
head(na.omit(sal_df$mean))

#levelplot
grid <- expand.grid(lon=lons, lat=lats)
cutpts <- c(0,5,10,15,20,25,30,35,40,45,50)
levelplot( sal_df$mean~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
           col.regions=(rev(brewer.pal(10,"RdBu"))))

# write out the dataframe as a .csv file
#csvpath <- "./var_ambientales"
#csvname <- "sst_mean.csv"
#csvfile <- paste(csvpath, csvname, sep="")
#write.table(na.omit(tmp_df02),csvfile, row.names=FALSE, sep=",")


#Create the raster with mean values
# from x,y,z-matrix
sss<- data.frame(cbind(sal_df$lon,sal_df$lat,sal_df$mean))
names(sss)<-c("lon","lat","sst_mean")
View(sss)

sss<- rasterFromXYZ(sss) #Convert first two columns as lon-lat and third as value  
plot(sss)


writeRaster(x = sss, 
            filename = "sss_mean",
            overwrite = TRUE, 
            format = 'GTiff',
            NAflag=-9999) 


GDALinfo("sss_mean.tif")
SSS<- raster("sss_mean.tif")
summary(SSS)
summary(SST, maxsamp = ncell(SST))
plot(SSS)

#ggplot2
info_df <- raster::as.data.frame(SSS, xy = TRUE)
str(info_df)
ggplot() +
  geom_raster(data = info_df , aes(x = x, y = y, fill = sss_mean)) +
  coord_equal() +
  scale_fill_gradientn(colours = terrain.colors(10))+
  ggtitle("Salinidad media") +
  theme(plot.title = element_text(lineheight=.8, face="bold",size = 20)) +
  theme(text = element_text(size=10))



##3. CHLOR_A--------------------------------------------------------------------------------------------

nc_file="./envar_variables/chlor_a_variable.nc"
nc_data <- nc_open(nc_file)
print(nc_data) 

lat_variable = 'latitude'
lon_variable = 'longitude'
time_variable ='z'
nc_variable = 'variable'

variable = ncvar_get(nc_data,nc_variable)
lats = ncvar_get(nc_data,lat_variable) 
lons = ncvar_get(nc_data,lon_variable) 
nlon<-dim(lons)
nlat<- dim(lats)
print(c(nlon,nlat))
times = ncvar_get(nc_data,time_variable) 
nt<- dim(times)

# reshape the array into vector
variable_long <- as.vector(variable)
length(variable_long)

# reshape the vector into a matrix
cl_A_mat <- matrix(variable_long, nrow=nlon*nlat, ncol=nt)
dim(cl_A_mat)
head(na.omit(cl_A_mat))

# create dataframe
lonlat <- as.matrix(expand.grid(lons,lats))
cla_df <- data.frame(cbind(lonlat,cl_A_mat))
names(cla_df ) <- c("lon","lat",
                    "clorAMar2013","clorAAbr2013","clorAMay2013",
                    "clorAMar2014","clorAAbr2014","clorAMay2014",
                    "clorAMar2015","clorAAbr2015","clorAMay2015",
                    "clorAMar2016","clorAAbr2016","clorAMay2016",
                    "clorAMar2017","clorAAbr2017","clorAMay2017",
                    "clorAMar2018","clorAAbr2018","clorAMay2018",
                    "clorAMar2019","clorAAbr2019","clorAMay2019")
# options(width=96)
head(na.omit(cla_df , 5))
dim(cla_df )

# get the annual mean
cla_df $mean <- apply(cla_df [3:23],1,mean) 
head(na.omit(cla_df $mean))


#levelplot
grid <- expand.grid(lon=lons, lat=lats)
cutpts <- c(0,0.2,0.4,0.6,0.8,1)
levelplot( cla_df$mean~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
           col.regions=(rev(brewer.pal(10,"RdBu"))))

# write out the dataframe as a .csv file
#csvpath <- "./var_ambientales"
#csvname <- "sst_mean.csv"
#csvfile <- paste(csvpath, csvname, sep="")
#write.table(na.omit(tmp_df02),csvfile, row.names=FALSE, sep=",")


# from x,y,z-matrix
cl_A<- data.frame(cbind(cla_df$lon,cla_df$lat,cla_df$mean))
names(cl_A)<-c("lon","lat","clor_A_mean")
View(cl_A)

cl_A<- rasterFromXYZ(cl_A) #Convert first two columns as lon-lat and third as value  
plot(cl_A)
res(cl_A)


writeRaster(x = cl_A, 
            filename = "clor_A_mean",
            overwrite = TRUE, 
            format = 'GTiff',
            NAflag=-9999)


GDALinfo("clor_A_mean.tif")
clorof<- raster("clor_A_mean.tif")
summary(clorof)
summary(clorof, maxsamp = ncell(clorof))
plot(clorof)

#ggplot2
info_df <- raster::as.data.frame(clorof, xy = TRUE)
str(info_df)
ggplot() +
  geom_raster(data = info_df , aes(x = x, y = y, fill = clor_A_mean)) +
  coord_equal() +
  scale_fill_gradientn(colours = terrain.colors(10))+
  ggtitle("media conc clorofila a") +
  theme(plot.title = element_text(lineheight=.8, face="bold",size = 20)) +
  theme(text = element_text(size=10))


