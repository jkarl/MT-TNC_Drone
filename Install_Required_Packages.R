# Install R library packages. 
# RTools should be downloaded and installed before running
# this script. Download from CRAN.

# You may need to update the path to the RTools bin directory (line 16)
# If you did not install it to the default c:/RTools location. 

install.packages("devtools")
library(devtools)

#Add Rtools to system path. 
Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), 
    "C:/RBuildTools/bin", #update as needed
    sep = ";"
  )
)

install.packages("tictoc", version= "1.2")
install.packages("parallel") #sometimes fails when installed via devtools::install_version
install_version("sp", version = "1.4-7")
install_version("rgdal", version = "1.4-3")
install_version("raster", version = "3.6-14")
install_version("terra", version = "1.5-21")
install_version("rgeos", version = "0.4-2")
install_version("abind", version = "1.4-5")
install_version("doParallel", version = "1.0.14")
install_version("foreach", version = "1.4.4")
install_version("gdalUtils", version = "2.0.1.14")
install_version("maptools", version = "0.9-5")
#install_version("spatial.tools", version = )
install_version("shotGroups", version = "0.7.4")
install_version("gam", version = "1.16")
install_version("TileManager", version = "0.4.1")
install_version("whitebox", version="2.2.0")
whitebox::install_whitebox()


install_version("tensorflow", version = 1.9)
install_version("keras", version = "2.2.0")

library(keras)
install_keras(tensorflow = "1.9.0-cpu") # if this fails, refer to Rstudio Tensorflow/Keras Documentation. 

#load packages to ensure libraries are functioning properly. 

library(abind)
library(doParallel)
library(foreach)
library(gam)
library(gdalUtils)
library(maptools)
library(parallel)
library(raster)
library(rgdal)
library(rgeos)
library(shotGroups)
library(sp)
library(TileManager)

library(keras)
library(tensorflow)
use_condaenv("r-tensorflow")
hello <- tf$constant('Hello, TensorFlow!')
sess <- tf$Session()
sess$run(tf$global_variables_initializer())
sess$run(hello)

