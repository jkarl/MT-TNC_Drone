## UAV data processing pipeline v0.5. April 17 2017
## Developed for TNC Montana by Terra Analytics
## Contact scott@terra-analytics.net

# Program Parameters -----------------------------------------------------

## Control how many threads are used for processing. This will vary by process due to
## memory requirements. The preset values have been tested on machine with
## 12 threads, 64 GB RAM, and a GTX980 TI GPU.

tile.threads <- 11 # Default 11
lastools.threads <- 10 #Default 10
unet.threads <- 8 #Default

unet.batchsize <- 4 #default 4 with GPU (perferred). Can be 6 with CPU. CPU will be slower.

remove.photogrammetry = TRUE #Set to FALSE if you want photogrammetry included in output.

# Working Directories -----------------------------------------------------

programDir <- "z:/Projects/TNC/2018Matador/Products/UAV_Pipeline" # Master Project Directory
binDir <- paste(programDir,"/bin", sep="")  #Binaries for processing. Do not change.
srcDir <- paste(programDir,"/src", sep="")  #Scripts for processing.  Do not change.

inputDirMaster <- "c:/work/uavInput" #input folder for raw UAV imagery. Must exists
outputDirMaster <- "c:/work/uavOutput" #output folder for processed data. Must Exist
scratchDir <- "c:/work/uavScratch" #temporary directory for processing. Should be local. Must exist.

# Define Executable Dependency Directories --------------------------------------------------------

path.photoscan  <- "c:/Program Files/Agisoft/PhotoScan Pro/" #directory for photoscan executable.
path.lastools <- "C:/Program Files/LasTools/bin/" # directory for lastool  executables.
path.arcgis.python <- "C:/Python27/ArcGIS10.6/" # directory with arcgis python executable. 

#add locations to %PATH%
Sys.setenv(
  PATH = paste(
    Sys.getenv("PATH"), 
    path.photoscan,
    path.lastools,
    path.arcgis.python,
    sep = ";"
  )
)

# Load Functions ----------------------------------------------------------
source(paste(srcDir,"/0_functions.R", sep=""))
if(find("check_lastools") != ".GlobalEnv") {
  print("Functions not loaded. Check source.")
  quit(save = no)
} else {
  print("Functions loaded...")
}

# Script and Model Locations --------------------------------------------------------
path.photoscan.script     <- paste(srcDir,"/1_AutoprocessingDrone.py", sep="") #Photoscan Processing
path.lastools.script      <- paste(srcDir,"/2_lasProcessing.bat",sep="") #Point Cloud Processing
path.arcpy.dem.script     <- paste(srcDir,"/3_dem_chm.py",sep="") #Create DEM
path.arcpy.sage.script    <- paste(srcDir,"/4_SageClass.py",sep="") #Sage AI classification post processing
path.arcpy.bg.script      <- paste(srcDir,"/5_BgClass.py",sep="") # Bare AI ground classification post processing
path.arcpy.updateCHM.script   <- paste(srcDir,"/6_UpdateCHM.py",sep="") #Update CHM from bareground AI model.
path.arcpy.rat.script      <- paste(srcDir,"/7_AddRAT_Classification.py",sep="") #Add RAT to final Classification 
path.arcpy.countsage.script   <- paste(srcDir,"/8_SageCount.py",sep="") #Count sagebrush and estimate height.
path.arcpy.grassht.script     <- paste(srcDir,"/9_GrassHt_Est.py",sep="") #Estimate Grass Height.
path.arcpy.buildstats.script  <- paste(srcDir,"/10_BuildStats.py",sep="") #Build final GIS Stats
path.gam.forb.model       <- paste(binDir,"/forb_gam.RData", sep="") # GAM model for estimating forb biomass.
path.keras.bg.model       <- paste(binDir,"/KerasBarePred_Unet.hd5",sep="") # HD5 Keras/Tensorflow model for bare ground
path.keras.sage.model     <- paste(binDir,"/KerasSagePred_Unet.hd5",sep="") # HD5 Keras/Tensorflow model for sagebrush

# Files/Depedendicies Check ------------------------------------------------------

#Fail out if bin directory is missing
if(!dir.exists(srcDir)) {
  print("ERROR: bin directory does not exist.")
  print(paste("Check",inputDirMaster))
  print("Exiting")
  quit(save = "no")
} else {
  print("bin directory exists...")
}

#Fail out if master input directory is missing
if(!dir.exists(inputDirMaster)) {
  print("ERROR: Master input directory does not exist.")
  print(paste("Check",inputDirMaster))
  print("Exiting")
  quit(save = "no")
} else {
  print("Master input directory exists...")
}

#Create master output directory if missing.
if(!dir.exists(outputDirMaster)) {
  print("ERROR: Master OutPut Directory Does not exist.")
  print(paste("Creating",outputDirMaster))
  dir.create(outputDirMaster)
} else {
  print("Master output directory exists...")
}

#Create scratch directory if missing.
if(!dir.exists(scratchDir)) {
  print("ERROR: Scratch directory does not exist.")
  print(paste("Creating",scratchDir))
  dir.create(scratchDir)
} else {
  print("Scratch directory exists....")
}

#Fail out if photoscan is missing
if(!file.exists(paste(path.photoscan,"photoscan.exe",sep=""))) {
  print("ERROR: Photoscan executable not found.")
  print("expecting it at:")
  print(path.photoscan)
  print("Check location or update master variables. exiting")
  quit(save = "no")
} else {
  print("Found Photoscan executable...")
}

#Fail out if photoscan script is missing
if(!file.exists(path.photoscan.script)) {
  print("ERROR: Photoscan script not found.")
  print("expecting it at:")
  print(path.photoscan.script)
  print("Check location or update master variables. exiting")
  quit(save = "no")
} else {
  print("Found Photoscan processing script...")
}

#Fail out if lastools is missing
lastools_exists <- check_lastools(path.lastools)
if(length(lastools_exists) > 0) {
  print("ERROR: The following LasTools executables were not found:")
  print(lastools_exists)
  print("Check location or update dependency list")
  quit(save = "no")
} else { 
  print("Found Lastools....")
}
rm(lastools_exists) # cleanup

#Fail out if lastools script is missing
if(!file.exists(path.lastools.script)) {
  print("ERROR: Lastools script not found.")
  print("expecting it at:")
  print(path.lastools.script)
  print("Check location or update master variables. exiting")
  quit(save = "no")
} else {
  print("Found LasTools processing script...")
}

#Fail out if Arc Python is missing
if(!file.exists(paste(path.arcgis.python,"python.exe",sep=""))) {
  print("ERROR: Arcpy Python executable not found.")
  print("expecting it at:")
  print(path.arcgis.python)
  print("Check location or update master variables. exiting")
  quit(save = "no")
} else {
  print("Found Arcpy Python executable...")
}

#Fail out if Arcpy dem/chm script is missing.
if(!file.exists(path.arcpy.dem.script)) {
  print("ERROR: Arcpy dem/chm script not found.")
  print("expecting it at:")
  print(path.arcpy.dem.script)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found Arcpy dem/chm processing script...")
}

#Fail out if hd5 model files for keras are missing.
if(!file.exists(path.keras.bg.model)) {
  print("ERROR: KerasBarePred_Unet.hd5 not found.")
  print("expecting it at:")
  print(binDir)
  print("Check location or update binDir. Exiting")
  stop()
} else {
  print("Found KerasBarePred_Unet.hd5 model file...")
}

if(!file.exists(path.keras.sage.model)) {
  print("ERROR: KerasSagePred_Unet.hd5 not found.")
  print("expecting it at:")
  print(binDir)
  print("Check location or update binDir. Exiting")
  stop()
} else {
  print("Found KerasSagePred_Unet.hd5 model file...")
}

#Fail out if Arcpy SageClass script is missing.
if(!file.exists(path.arcpy.sage.script)) {
  print("ERROR: Arcpy sage classification script not found.")
  print("expecting it at:")
  print(path.arcpy.sage.script)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found Arcpy sage classification processing script...")
}

#Fail out if Arcpy bareground class script is missing.
if(!file.exists(path.arcpy.bg.script)) {
  print("ERROR: Arcpy bareground classification script not found.")
  print("expecting it at:")
  print(path.arcpy.bg.script)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found Arcpy bareground classification processing script...")
}

#Fail out if Arcpy update CHM script is missing.
if(!file.exists(path.arcpy.updateCHM.script)) {
  print("ERROR: Arcpy update CHM script not found.")
  print("expecting it at:")
  print(path.arcpy.updateCHM.script)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found Arcpy Update CHM processing script...")
}

#Fail out if Arcpy update RAT script is missing.
if(!file.exists(path.arcpy.rat.script)) {
  print("ERROR: Arcpy update RAT script not found.")
  print("expecting it at:")
  print(path.arcpy.rat.script)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found Arcpy Update RAT processing script...")
}

#Fail out if Arcpy Count Sage script is missing.
if(!file.exists(path.arcpy.countsage.script)) {
  print("ERROR: Arcpy Count Sage script not found.")
  print("expecting it at:")
  print(path.arcpy.countsage.script)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found Arcpy Count Sage processing script...")
}

#Fail out if Arcpy Grass Height script is missing.
if(!file.exists(path.arcpy.grassht.script)) {
  print("ERROR: Arcpy Grass Height script not found.")
  print("expecting it at:")
  print(path.arcpy.grassht.script)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found Arcpy grass height processing script...")
}

#Fail out if Arcpy Build Stats script is missing.
if(!file.exists(path.arcpy.buildstats.script)) {
  print("ERROR: Arcpy Grass Height script not found.")
  print("expecting it at:")
  print(path.arcpy.buildstats.script)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found Arcpy Build Stats processing script...")
}

#Fail out if forb GAM model is missing.
if(!file.exists(path.gam.forb.model)) {
  print("ERROR: Forb GAM model not found.")
  print("expecting it at:")
  print(path.gam.forb.model)
  print("Check location or update master variables. Exiting")
  quit(save = "no")
} else {
  print("Found forb GAM model...")
}


# 0. Pipeline. Start. Check inputs, copy files ----------------------------------------------------------------

#Check input directory for folders.
#The Subfolder name will be output project folder.
#raw imagery in the master input directory will not be accepted.
setwd(inputDirMaster)
dlist.in <- list.dirs(path = ".")[-1] #Get list of directories, minus ".".
if(length(dlist.in) == 0 ) {
  print("No input directories found. Aborting.")
  print("raw imagery must be put into a subdirectory for processing")
  quit(save = "no")
}

print(paste("Found",length(dlist.in),"inputs. Processing sequentially"))

#Iterate through inputs sequentially using n....
for(n in seq_along(dlist.in)) {
  setwd(inputDirMaster)
  project.name <- gsub("./", "", dlist.in[n])
  print(paste("Opening Project Folder:", project.name))
  setwd(dlist.in[n])
  
  project.subdirs <- list.dirs(path = ".")[-1] #Check for subdirectories. Break if found
  if(length(project.subdirs) > 0) {
    print(paste("Subdirectories found. All imagery must be in folder", project.subdirs))
    print(paste("Place all imagery in", project.subdirs, ", delete subdirectories, and restart pipeline."))
    quit(save = "no")
  }
  rm(project.subdirs) #cleanup
  
  raw.files <- grep(".raw$|.RAW$",x = dir()) #check for raw image files. Return error if found.
  if(length(raw.files) > 0) {
    print("ERROR: RAW image files found. These images may not be compatible with Photoscan")
    print("Proceeding with processing, but expect photogrammetry to break")
  }
  
  tif.files <- grep(".TIF$|.tif$|.tiff$|.TIFF$",x = dir()) #check for tif image files. Return warning if found.
  if(length(tif.files) > 0) {
    print(paste(length(tif.files),"Tif files found. Processing pipeline built for jpeg inputs."))
    print("The pipeline should handle tifs, but if it breaks, may need to revisit photoscan script.")
  }
  
  jpeg.files <- grep(".JPG$|.jpg$|.jpeg$|.JPEG$",x = dir())
  if(length(jpeg.files) > 0) {
    print(paste(length(jpeg.files),"jpeg files found. Proceeding."))
  }
  
  all.image.files <-  dir()[as.numeric(c(raw.files,tif.files,jpeg.files))]
  rm(raw.files,tif.files,jpeg.files) #cleanup
  
  if(length(all.image.files) == 0) {
    print("No Standard Image Files Found. Aborting")
    quit(save = "no")
  }
  
  ## Next, check to see if output folder exists.
  ## If exists, break and tell user to remove.
  project.outdir <- paste(outputDirMaster,"/",project.name,sep="")
  if(dir.exists(project.outdir)) {
    print("Output Folder exists:")
    print(project.outdir)
    print ("Remove output directory or rename input.")
    quit(save = "no")
  }
  
  ## Next, check to see if scratch/processing directory exists. 
  ## If it exists, remove, give warning.
  
  project.scratch <-  paste(scratchDir,"/",project.name,sep="")
  if(dir.exists(project.scratch)) {
    print("Found intermediate data. Deleting.")
    unlink(project.scratch, recursive = TRUE, force = TRUE)
  }
  
  #Create Scratch Dir and required subdirs
  dir.create(project.scratch)
  dir.create(paste(project.scratch,"/rawImagery",sep=""))
  dir.create(paste(project.scratch,"/scratch",sep=""))
  dir.create(paste(project.scratch,"/photogrammetry",sep=""))
  dir.create(paste(project.scratch,"/reports",sep=""))
  dir.create(paste(project.scratch,"/logs",sep=""))
  dir.create(paste(project.scratch,"/products",sep=""))
  
  #copy input images to scratch directory for processing.
  files.dest <- paste(project.scratch,"/rawImagery/",all.image.files,sep="")
  print("Copying input imagery to scratch location...")
  file.copy(all.image.files,files.dest)
  rm(files.dest,all.image.files)
  
  # 1. Pipeline. Run Photoscan on raw imagery inputs ------------------------
  
  setwd(project.scratch)
  photoscan.args  <- paste("-r ",
                           path.photoscan.script,
                           " ",
                           project.scratch,
                           "/",
                           sep="")
  
  print(paste("(1 of 16) Starting Photogrammetry proceessing at",date(),"... this could take a while"))
  file.create("logs/photoscanProcessing.txt")
  file.create("logs/photoscanProcessingError.txt")
  
  system2("photoscan.exe",
          args = photoscan.args, 
          stdout = "logs/photoscanProcessing.txt",
          stderr = "logs/photoscanProcessingError.txt")  
  print(paste("Photogrammetry complete at",date())) 
  
  rm(photoscan.args) #cleanup
  
  # 2. Pipeline. Run LasTools (process pointcloud). -------------------------
  print("(2 of 16) Starting Point Cloud Filtering")
  setwd(project.scratch)
  
  las.scratch <- paste(project.scratch,"/scratch/las",sep="")
  dir.create(las.scratch)
  
  print("Copying pointcloud")
  file.copy("photogrammetry/pointcloud.laz","scratch/las/pointcloud.laz")
  
  setwd(las.scratch)
  lastools.args <- paste("pointcloud.laz",lastools.threads)
  
  file.create("../../logs/lasToolsProcessing.txt")
  file.create("../../logs/lasToolsProcessingError.txt")
  
  print(paste("Starting LasTools at",date(),"... this could take a while"))
  system2(path.lastools.script,
          args = lastools.args,
          stdout = "../../logs/lasToolsProcessing.txt",
          stderr = "../../logs/lasToolsProcessingError.txt")
  
  #Move processed files to products directory.
  file.rename("pointcloud_clean.laz","../../products/pointcloud.laz")
  file.rename("pointcloud_ground.laz","../../products/pointcloud_ground.laz")
  
  #remove raw pointcloud file
  file.remove("pointcloud.laz")
  
  #copy remaining files to scratch ("up one directory")
  files.remain <- dir()
  files.remain2 <- paste("../",files.remain,sep="")
  file.rename(files.remain,files.remain2)
  rm(files.remain,files.remain2)
  
  #Remove las directory
  setwd(project.scratch)
  unlink(las.scratch, recursive = TRUE, force = TRUE)
  rm(las.scratch,lastools.args) #clean up
  print(paste("Point clould filtering complete at",date(),"."))
  
  
  
  # 3. Pipeline. Create AOI Polygon. ----------------------------------------
  setwd(project.scratch)
  print("(3 of 16) Create AOI Polygon")
  
  generate_aoi(project.scratch)
  
  setwd(project.scratch)
  if(!file.exists("scratch/aoiPoly.shp")) {
    print("ERROR: aoipoly.shp not found in scratch directory. Aborting.")
    quit(save = "no")
  } else {
    print("AOI Polygon created. Moving on...")
  }
  
  # 4. Pipeline. Create Bare Earth and Canopy Height Models. ----------------
  #This script also clips the photogrammetry products to the aoi.
  setwd(project.scratch)
  
  #Call ArcGIS python executable directly, as conda will likely override.
  arcpy.exec.string <- paste(path.arcgis.python,"python.exe",sep="") 
  arcpy.args <- paste(path.arcpy.dem.script,project.scratch)
  print("(4 of 16) Generating CHM/DEM...")
  
  file.create("logs/Arcpy_dem_chm_Processing.txt")
  file.create("logs/Arcpy_dem_chm_processing_errors.txt")
  
  system2(arcpy.exec.string,
          args = arcpy.args,
          stdout = "logs/Arcpy_dem_chm_Processing.txt",
          stderr = "logs/Arcpy_dem_chm_processing_errors.txt")
  
  prod_check <- check_prelim_products(project.scratch)
  if(length(prod_check) > 0) {
    print("ERROR: Not all products were found. Missing:")
    print(prod_check)
    print("check dem/chm processing log errors. Exiting.")
    quit(save = "no")
  } else {
    print("DEM/CHM creation complete. Moving on...")
  }
  
  rm(arcpy.exec.string,arcpy.args,prod_check) #clean up
  
  # 5. Pipeline. Tile imagery and chm for Deep Learning ---------------------
  setwd(project.scratch)
  print("(5 of 16) Tiling imagery and CHM to pass to deep learning frameworks...")
  
  # Splits ortho and chm into 1024 x 1024 tiles for DL prediction.
  tile_data(project.scratch,tile.threads)
  
  # 6. Pipeline. Deep Learning Prediciton: Sagebrush and Bareground---------------------
  setwd(project.scratch)
  print("(6 of 16) Generating sagebrush and bareground cover models...")
  pred_unet(project.scratch,unet.threads,unet.batchsize)
  
  # 7. Pipeline. Postprocessing Sagebrush and Bareground---------------------
  setwd(project.scratch)
  print("(7 of 16) Post processing sagebrush and bareground models ...")
  
  arcpy.exec.string <- paste(path.arcgis.python,"python.exe",sep="") 
  arcpy.args <- paste(path.arcpy.sage.script,project.scratch)
  
  file.create("logs/Arcpy_SageClass_Processing.txt")
  file.create("logs/Arcpy_SageClass_processing_errors.txt")
  
  system2(arcpy.exec.string,
          args = arcpy.args,
          stdout = "logs/Arcpy_SageClass_Processing.txt",
          stderr = "logs/Arcpy_SageClass_processing_errors.txt")
  
  arcpy.exec.string <- paste(path.arcgis.python,"python.exe",sep="") 
  arcpy.args <- paste(path.arcpy.bg.script,project.scratch)
  
  file.create("logs/Arcpy_BareClass_Processing.txt")
  file.create("logs/Arcpy_BareClass_processing_errors.txt")
  
  system2(arcpy.exec.string,
          args = arcpy.args,
          stdout = "logs/Arcpy_BareClass_Processing.txt",
          stderr = "logs/Arcpy_BareClass_processing_errors.txt")
  
  rm(arcpy.exec.string,arcpy.args) #clean up
  
  # 8. Pipeline. Update CHM ---------------------
  setwd(project.scratch)
  print("(8 of 16) Update canopy height model from Deep Learning Models ...")
  
  arcpy.exec.string <- paste(path.arcgis.python,"python.exe",sep="") 
  arcpy.args <- paste(path.arcpy.updateCHM.script ,project.scratch)
  
  file.create("logs/Arcpy_UpdateCHM_Processing.txt")
  file.create("logs/Arcpy_UpdateCHM_processing_errors.txt")
  
  system2(arcpy.exec.string,
          args = arcpy.args,
          stdout = "logs/Arcpy_UpdateCHM_Processing.txt",
          stderr = "logs/Arcpy_UpdateCHM_processing_errors.txt")
  
  file.copy("scratch/chm2.tif","products/chm2.tif")
  rm(arcpy.exec.string, arcpy.args) #cleanup
  
  # 9. Pipeline. Classify Cover---------------------
  setwd(project.scratch)
  print("(9 of 16) Generate final cover model ...")
  gen_cover(project.scratch) # combine cover class models into one raster.
  
  #add Raster Attribute Table
  #add table files to scratch dir before executing join
  src.files <- dir(binDir,pattern="^rat_class", full.names = TRUE)
  dest.files <- paste("scratch/",dir(binDir,pattern="rat_class"),sep="")
  file.copy(src.files,dest.files,overwrite = TRUE)
  rm(src.files,dest.files) #cleanup
  
  arcpy.exec.string <- paste(path.arcgis.python,"python.exe",sep="") 
  arcpy.args <- paste(path.arcpy.rat.script, project.scratch)
  
  system2(arcpy.exec.string,
          args = arcpy.args,
          stdout = "logs/Arcpy_addRAT_Processing.txt",
          stderr = "logs/Arcpy_addRAT_processing_errors.txt")
  
  rm(arcpy.exec.string, arcpy.args) #cleanup
  
  copy_cover(project.scratch) #copy cover class to products
  
  # 10. Pipeline. Get class and sagebrush Statistics  ---------------------   
  print("(10 of 16) Extract cover statistics ...")
  class.stats <- get_class_stats(project.scratch) #Gets % cover for bare ground and sagebrush
  
  arcpy.exec.string <- paste(path.arcgis.python,"python.exe",sep="") 
  arcpy.args <- paste( path.arcpy.countsage.script, paste(project.scratch,"/scratch",sep=""))
  
  #Generates a point shapefile that estimates number of sagebrush plants and their height.
  system2(arcpy.exec.string,
          args = arcpy.args,
          stdout = "logs/Arcpy_countSage_Processing.txt",
          stderr = "logs/Arcpy_countSage_processing_errors.txt")
  sage.stats <- count_sage(project.scratch) # returns vector of sagebrush heights.
  sage.density <- round(length(sage.stats)/ sum(class.stats$area_m2) * 4046.86,0) # plants per acre
  class.list <- list(class.stats = class.stats, 
                     sage.ht = sage.stats, 
                     sage.density = sage.density) #Roll up stats into list)
  rm(class.stats, sage.stats,sage.density) #clean up
  
  # 11. Pipeline. Estimate Herb. Biomass.  ---------------------     
  setwd(project.scratch)
  print("(11 of 16) Generate biomass model and statistics ...")
  biomass.stats <- plot_biomass(project.scratch)
  
  # 12. Pipeline. Estimate Grass Height.  ---------------------     
  setwd(project.scratch)
  print("(12 of 16) Generate grass height model and statistics ...")
  
  arcpy.exec.string <- paste(path.arcgis.python,"python.exe",sep="") 
  arcpy.args <- paste( path.arcpy.grassht.script, paste(project.scratch,"/scratch",sep=""))
  
  #Generates a point shapefile that estimates number of sagebrush plants and their height.
  system2(arcpy.exec.string,
          args = arcpy.args,
          stdout = "logs/Arcpy_grassht_Processing.txt",
          stderr = "logs/Arcpy_grassht_processing_errors.txt")
  
  grass.stats <- get_grass_stats(project.scratch)
  
  # 13. Pipeline. Build Final Statistics for GIS Files  ---------------------
  setwd(project.scratch)
  print("(13 of 16) Build final GIS statistics ...")
  
  arcpy.exec.string <- paste(path.arcgis.python,"python.exe",sep="") 
  arcpy.args <- paste(path.arcpy.buildstats.script, paste(project.scratch,"/products",sep=""))
  
  #Build pyramids and statistics for all raster layers.
  system2(arcpy.exec.string,
          args = arcpy.args,
          stdout = "logs/Arcpy_buildstats_Processing.txt",
          stderr = "logs/Arcpy_buildstats_processing_errors.txt")
  
  #copy cameras and aoi poly
  copy_cameras(project.scratch)
  
  #Copy Project map file for ArcGIS Users
  file.copy(paste(binDir,"/Project.mxd", sep=""),paste(project.scratch,"/Products/Project.mxd",sep=""),overwrite = TRUE)
  
  rm(arcpy.exec.string,arcpy.args)
  
  # 14. Pipeline. Finalize Reports  ---------------------
  print("(14 of 16) Finalize Reports ...")
  
  get_sagebrush_plot(project.scratch)
  gen_table_data(project.scratch)
  
  setwd(project.scratch)
  file.rename("logs/photoscan.pdf","reports/photogrammetry.pdf")
  
  
  # 15. Pipeline. Remove Scratch and  photogrammetry---------------------
  print("(15 of 16) Remove Temporary Files ...")
  setwd(project.scratch)
  unlink("scratch/", recursive = TRUE, force=TRUE)
  
  if(remove.photogrammetry == TRUE) {
    unlink("photogrammetry", recursive = TRUE, force = TRUE)
  }
  
  # 16. Pipeline. Copy project to final output folder ---------------------
  print("(16 of 16) Copy Files to final location ...")
  setwd(project.scratch)
  setwd("..")
  fc.result <- file.copy(dlist.in[n], recursive = TRUE, outputDirMaster)
  
  #remove scratch and input directory
  if(fc.result == TRUE) {
    setwd(inputDirMaster)
    unlink(dlist.in[n], recursive = TRUE, force=TRUE)
  } else { 
    print("Something went wrong we copying to Output Directory. Aborting...")
    quit(save = "no")    
  }
  
  rm(fc.result)
  print(paste("Plot Processing Completed at",date(),"Moving on..."))
}

print(paste("All processing complete at",date(),"Exiting."))










