"""
Script: R_Workflow_Modified.R
Author: Jason Karl, modified from original (non-functional) work of Scott Morford
Date: 14 Dec 2023
########################################################################################################
Description:
This second part of the UAV Data Processing Workflow uses the products created in the first part to create
estimates of rangeland monitoring indicators. First, the Orthomosaic and CHM are merged together and then
split into 1024x1024 pixel tiles. Tiling of the input layers is necessary for fast/efficient model
predictions with TensorFlow. Once the tiles have been created, parallel processing is set up to allow
multiple tiles to be processed at the same time. Each tile is passed to two TensorFlow models, one that
predicts sagebrush cover, and the other that predicts bare ground. The resulting model prediction tiles
are then reassembled into full site rasters and clipped to the AOI boundary. Model probability values are
converted to binary classes (present or absent) to create final sagebrush and bare ground models. Areas
predicted to be bare ground but that have a canopy height greater than zero are reclassified as herbaceous
vegetation. Then the land cover class rasters (sagebrush, bare ground, herbaceous) are merged to create the
final land cover raster. Areas designated as herbaceous vegetation are then run through the herbaceous
biomass generalized additive model (GAM) to estimate herbaceous biomass. Finally, summary statistics and
tables are created from the final products.
########################################################################################################
Pre-requisites:
  1. This script requires that the Metashape processing be completed and the canopy height model
     be constructed before running.
  2. Run the Install_Required_Packages.R script to install the requisite libraries and the
     correct versions.
  3. This script will check for the correct directory structure that was created when running the
     Metashape processing script. If the directory structure is not present, the script will abort.
########################################################################################################
Notes
  1. This script is somewhat fragile and can crash for various reasons. I have found the best luck
     running it in sections from within R Studio and not trying to batch/source the entire thing at once.
  2. Memory management is a big deal with this script. Multiple efforts have been made to reduce the
     memory usage of the script as it runs, but it will still accumulate a large amount of temporary files
     in R's temporary storage location (on my computer, C:/Users/jkarl/AppData/Local/R). Make sure to
     run the clean_temp(tempdir()) function when you are done running a plot to remove those temp files.
########################################################################################################
"""

## -----------------------------------------------------------------------------------
# Set the default project directory and name if not using command line arguments
## Uncomment these and adjust if running within RStudio
##
override.dir <- "C:/Users/jkarl/Downloads/2022/LH1"
override.name <- "LH1"

# Directory where scripts are tensor flow models are stored.
programDir <- "C:/Users/jkarl/OneDrive - University of Idaho/Documents/GitHub/MT-TNC_Drone"


## -----------------------------------------------------------------------------------
# Set up the project environment and parse input arguments ---------------------------
if (length(commandArgs(trailingOnly=TRUE))>0) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) < 2) {
    cat("You ran the program with ", args,"\n")
    quit("Usage: R_Workflow.R <project.directory> <project.name>", save="no")
  }
} else{
  ## Set default values for arguments
  args <- c(override.dir,override.name)
}

project.dir <- args[1]
project.name <- args[2]

## Check to see if project directory exists
if (!file.exists(project.dir)) {
  stop(paste("Project directory",project.dir,"does not exist. Exiting.",sep=" "))
}

## Check for scratch directory. Create if does not exist
scratch.dir <- file.path(project.dir,"/scratch")
if (!file.exists(scratch.dir)) {
  print(paste0("Project directory ",project.dir," does not exist. Creating it."))
  dir.create(scratch.dir)
}
setwd(project.dir)


## -----------------------------------------------------------------------------------
# Set modeling parameters
## Control how many threads are used for processing. This will vary by process due to
## memory requirements. The preset values have been tested on machine with
## 12 threads, 64 GB RAM, and a GTX980 TI GPU.

tile.threads <- 6 # Default 11
unet.threads <- 8 #Default

unet.batchsize <- 4 #default 4 with GPU (perferred). Can be 6 with CPU. CPU will be slower.


binDir <- paste(programDir,"/models", sep="")  #Binaries for processing. Do not change.
srcDir <- paste(programDir,"/models", sep="")  #Scripts for processing.  Do not change.

## Script and Model Locations
path.gam.forb.model       <- paste(binDir,"/forb_gam.RData", sep="") # GAM model for estimating forb biomass.
path.keras.bg.model       <- paste(binDir,"/KerasBarePred_Unet.hd5",sep="") # HD5 Keras/Tensorflow model for bare ground
path.keras.sage.model     <- paste(binDir,"/KerasSagePred_Unet.hd5",sep="") # HD5 Keras/Tensorflow model for sagebrush

## Memory limit for raster processing
mem.limit <- 10e+09 #10e+09


# Load Functions ----------------------------------------------------------
library(keras)
use_condaenv("r-tensorflow")

source(paste(programDir,"/functions.R", sep=""))
if(find("tile_data") != ".GlobalEnv") {
  print("Functions not loaded. Check source.")
  quit(save = no)
} else {
  print("Functions loaded...")
}

# Start timing
library(tictoc)
tic("Total elapsed time")




# ---------------------------------------------------------------------------------

# Splits ortho and chm into 1024 x 1024 tiles for DL prediction. -----------
print("Tiling imagery and CHM to pass to deep learning frameworks...")
tic()
tile_data(project.dir,project.name,tile.threads)
toc()

# Run Shrub and Bareground Model predictions ---------------------
print("Running TensorFlow models to predict shrub and bare ground cover...")
tic()
reproject_aoi(project.dir)
pred_unet(scratch.dir, unet.threads, unet.batchsize)
toc()

# Clean the shrub and bare ground rasters ---------------------
print("Cleaning up the TensorFlow model outputs...")
tic()
clean_models(project.dir)
toc()

# Construct initial cover type class raster ------------------
print("Assembling the cover type class raster")
tic()
gen_cover(project.dir)
class.stats <- get_class_stats(project.dir) #Gets % cover for bare ground and sagebrush
class.list <- list(class.stats = class.stats)
toc()

# Update CHM based on class predictions ------------------
print("Updating the canopy height model for TensorFlow bare ground...")
tic()
updateCHM(project.dir, project.name)
toc()

# Grass Height Estimation ----------------------------
print("Estimating grass height...")
tic()
grass_ht(project.dir)
grass.stats <- get_grass_stats(project.dir)
toc()

# DID NOT IMPLEMENT:: Sagebrush count --------------------------
#print("Identifying individual sagebrush (as best we can)...")

# DID NOT IMPLEMENT:: Sagebrush height ---------------------------
#print("Estimating sagebrush height...")
#tic()
#get_sagebrush_plot(project.dir)
#toc()

# Herbaceous biomass ---------------------------
print("Estimating herbaceous biomass...")
tic()
biomass.stats <- plot_biomass(project.dir)
toc()

# Summary statistics -------------------------
print("Calculating summary statistics and final output report...")
tic()
gen_table_data(project.dir)
toc()
# Use raster::freq on the raster classification to get % cover of classes

# Final Clean-up of Scratch files
print("Cleaning up temporary files...")
tic()
clean_scratch(scratch.dir, layers=TRUE, tiles=TRUE)
clean_temp(tempdir()) # remove any left over temp .tif files that eat system memory. Run this if the script crashes.
rm(list=ls())
gc()
toc()

# Finished!
print(paste("Finished running plot",override.name,sep=" "))
toc()
