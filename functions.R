###################################################################
## functions.R
## This script contains the main functions that the data analysis
## workflow uses to process the photogrammetry products into
## indicator estimates and produce the output graphs and reports
## This script should be located in the same directory as the
## R_Workflow.R script
###################################################################


# FUNCTION: reproject the AOI to the same as the project rasters, save in scratch directory -----------------
reproject_aoi <- function (x, epsg="32100") { # x is project directory, epsg is output coordinate system
  
  # Load the AOI and reproject it
  aoi <- terra::vect(paste(x,"photogrammetry","aoi.shp",sep="/"),layer="aoi")
  terra::crs(aoi) <- "epsg:4326"
  aoi.proj <- terra::project(aoi, paste("EPSG:",epsg,sep=""))
  #aoi.proj <- sp::spTransform(aoi,sp::CRS(paste("+init=epsg:",epsg,sep="")))
  terra::writeVector(aoi.proj,paste(x,"scratch",sep="/"),"aoi_proj",filetype="ESRI Shapefile",overwrite=TRUE)
  
  # get the raster extent and resolution from the bgpred.tif model
  #template.raster <- raster(paste(x,"scratch","bgpred.tif",sep="/"))
  
  # create the aoi raster and write it out
  #aoi.raster <- raster(ext=extent(template.raster),resolution=res(template.raster), crs=CRS(paste("+init=epsg:",epsg,sep="")), vals=1)
  #aoi.raster <- mask(aoi.raster,aoi.proj)
  #writeRaster(aoi.raster,paste(x,"scratch","aoi.tif",sep="/"))
  
  return(TRUE)
}

# FUNCTION: Raster if/else function (e.g., ArcGIS con function) --------------
Con=function(condition, trueValue, falseValue){
  return(condition * trueValue + (!condition)*falseValue)
}

# FUNCTION: Tiles RGB-CHM into 1024x1024x4 (RGBZ) tiles to pass to keras. -----------------------------------------------------------------
tile_data <- function(x,name,y) { #x is project.dir, name is project.name, y is number of cores for parallel processing
  require(terra)
  #require(raster)
  #require(gdalUtils)
  #require(spatial.tools)
  require(TileManager)
  #require(doParallel)  
  #require(foreach) 
  #require(rgdal)
  require(sp)
  
  #setwd(x)
  nCores <- y
  #rasterOptions(tmpdir = scratch.dir, maxmemory = mem.limit)
  
  print("tile_data: Creating raster brick")
  ortho <- c(terra::rast(paste(x,"/products/",name,"_ortho.tif",sep="")))
  ortho <- terra::subset(ortho, 1:3)
  
  print("tile_data: Converting nodata values to zero")
  #TensorFlow does not like NA values.
  ortho[is.na(ortho)] <- 0 
  gc()
  
  # realign chm. Sometimes extent is a little wonky coming out of PDAL.
  print("tile_data: Aligning CHM to Ortho")
  chm <- terra::resample(terra::rast(paste(x,"/products/",name,"_chm.tif",sep="")), ortho, method="bilinear")
  terra::writeRaster(chm,paste(x,"/scratch/chm_align.tif",sep=""),overwrite=T)
  #gc()
  
  #rescale chm 0 to 1 with max height (1) equal 2m.
  print("tile_data: rescaling CHM")
  chm <- chm  / 2
  chm[chm > 1] <- 1
  
  ## test add to address mislabelling at edge
  chm[is.na(chm)] <- 0
  gc()
  
  
  #combine rasters into a rgbz brick.
  print("tile_data: Combining all rasters into the brick")
  combined_raster <- c(ortho[[1]],ortho[[2]],ortho[[3]],chm)
  rm(chm,ortho)
  #unlink(paste(x,"/scratch/chm_align.tif",sep=""), force = TRUE)
  gc()
  #print("test1")
  #combined_raster <- readAll(combined_raster)
  #gc()
  print("Setting extent")
  ext <- terra::ext(combined_raster)
  r <- terra::res(combined_raster)
  ext <- terra::extend(ext,c(0,2000*r[1],2000*r[2],0))
  combined_raster <- terra::extend(combined_raster,ext)
  
  #combined_raster <- modify_raster_margins(combined_raster, extent_delta = c(0,2000,2000,0), value = NA)
  
  
  
  gc()
  
  
  print("tile_data: setting up tile scheme")
  ts.ortho <- TileManager::tileScheme(raster::raster(combined_raster),
                                      tiledim=c(1000,1000),
                                      cells=TRUE,
                                      buffer = 12,
                                      bufferspill = FALSE,
                                      removeEmpty = TRUE)
  
  buff.poly <- ts.ortho@buffs #      $buffPolygons
  pred.poly <- ts.ortho@tiles #      $tilePolygons
  n.tiles <- length(buff.poly) #dim(buff.poly)[1]
  poly.id <- as.character(names(buff.poly))  #as.character(buff.poly[[1]])
  
  if(dir.exists(paste(x,"/scratch/tiles",sep=""))) {
    unlink(paste(x,"/scratch/tiles",sep=""), recursive=TRUE, force=TRUE)
  }
  
  dir.create(paste(x,"/scratch/tiles",sep=""))
  save(buff.poly,pred.poly, file = paste(x,"/scratch/TilePolys.RData",sep=""))
  
  print("tile_data: tiling raster brick. This may take a while...")
  
  for (j in 1:n.tiles) {
    print(paste("Clipping tile ",j))
    t.extent <- bbox(buff.poly[[j]])
    
    t.rast <- terra::crop(combined_raster,c(t.extent[1,1],t.extent[1,2],t.extent[2,1],t.extent[2,2]))
    terra::writeRaster(t.rast, filename = paste("scratch/tiles/",poly.id[j],".tif",sep=""),
                       filetype = "GTiff", overwrite=TRUE, datatype="FLT4S")
    rm(t.extent,t.rast)
    gc()
  }
  
  #setwd("scratch/tiles")
  unlink(dir(path=paste(x,"/scratch/tiles",sep=""), pattern=".aux.xml$"))
  
  print("tile_data: Done tiling. Whew!")
}

# FUNCTION U-net 1024 - prep layers for Tensorflow modeling -------------------------------------
get_unet_1024 <- function(input_shape = c(1024, 1024, 4),  #input shape is 1024,1024,4 - RGB(CHM)
                          num_classes = 1) { 
  
  inputs <- layer_input(shape = input_shape)
  # 1024
  
  print("get_unet_1024: down1")
  down1 <- inputs %>%
    layer_conv_2d(filters = 8, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 8, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") 
  down1_pool <- down1 %>%
    layer_max_pooling_2d(pool_size = c(2, 2), strides = c(2, 2))
  # 512
  
  print("get_unet_1024: down2")
  down2 <- down1_pool %>%
    layer_conv_2d(filters = 16, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 16, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") 
  down2_pool <- down2 %>%
    layer_max_pooling_2d(pool_size = c(2, 2), strides = c(2, 2))
  # 256
  
  print("get_unet_1024: down3")
  down3 <- down2_pool %>%
    layer_conv_2d(filters = 32, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 32, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") 
  down3_pool <- down3 %>%
    layer_max_pooling_2d(pool_size = c(2, 2), strides = c(2, 2))
  # 128
  
  print("get_unet_1024: down4")
  down4 <- down3_pool %>%
    layer_conv_2d(filters = 64, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 64, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") 
  down4_pool <- down4 %>%
    layer_max_pooling_2d(pool_size = c(2, 2), strides = c(2, 2))
  # 64
  
  print("get_unet_1024: down5")
  down5 <- down4_pool %>%
    layer_conv_2d(filters = 128, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 128, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") 
  down5_pool <- down5 %>%
    layer_max_pooling_2d(pool_size = c(2, 2), strides = c(2, 2))
  #32
  
  print("get_unet_1024: down6")
  down6 <- down5_pool %>%
    layer_conv_2d(filters = 256, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 256, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") 
  down6_pool <- down6 %>%
    layer_max_pooling_2d(pool_size = c(2, 2), strides = c(2, 2))
  #16
  
  print("get_unet_1024: down7")
  down7 <- down6_pool %>%
    layer_conv_2d(filters = 512, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 512, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") 
  down7_pool <- down7 %>%
    layer_max_pooling_2d(pool_size = c(2, 2), strides = c(2, 2))
  # 8
  
  print("get_unet_1024: center")
  center <- down7_pool %>%
    layer_conv_2d(filters = 1024, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 1024, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") 
  # center
  
  print("get_unet_1024: up7")
  up7 <- center %>%
    layer_upsampling_2d(size = c(2, 2)) %>%
    {layer_concatenate(inputs = list(down7, .), axis = 3)} %>%
    layer_conv_2d(filters = 512, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 512, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 512, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu")
  # 16
  
  print("get_unet_1024: up6")
  up6 <- up7 %>%
    layer_upsampling_2d(size = c(2, 2)) %>%
    {layer_concatenate(inputs = list(down6, .), axis = 3)} %>%
    layer_conv_2d(filters = 256, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 256, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 256, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu")
  # 32
  
  print("get_unet_1024: up5")
  up5 <- up6 %>%
    layer_upsampling_2d(size = c(2, 2)) %>%
    {layer_concatenate(inputs = list(down5, .), axis = 3)} %>%
    layer_conv_2d(filters = 128, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 128, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 128, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu")
  # 64
  
  print("get_unet_1024: up4")
  up4 <- up5 %>%
    layer_upsampling_2d(size = c(2, 2)) %>%
    {layer_concatenate(inputs = list(down4, .), axis = 3)} %>%
    layer_conv_2d(filters = 64, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 64, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 64, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu")
  # 128
  
  print("get_unet_1024: up3")
  up3 <- up4 %>%
    layer_upsampling_2d(size = c(2, 2)) %>%
    {layer_concatenate(inputs = list(down3, .), axis = 3)} %>%
    layer_conv_2d(filters = 32, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 32, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 32, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu")
  #256
  
  print("get_unet_1024: up2")
  up2 <- up3 %>%
    layer_upsampling_2d(size = c(2, 2)) %>%
    {layer_concatenate(inputs = list(down2, .), axis = 3)} %>%
    layer_conv_2d(filters = 16, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 16, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 16, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu")
  #512
  
  print("get_unet_1024: up1")
  up1 <- up2 %>%
    layer_upsampling_2d(size = c(2, 2)) %>%
    {layer_concatenate(inputs = list(down1, .), axis = 3)} %>%
    layer_conv_2d(filters = 8, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 8, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu") %>%
    layer_conv_2d(filters = 8, kernel_size = c(5, 5), padding = "same") %>%
    layer_batch_normalization() %>%
    layer_activation("relu")
  #1024
  
  print("get_unet_1024: classify")
  classify <- layer_conv_2d(up1,
                            filters = num_classes, 
                            kernel_size = c(1, 1),
                            activation = "sigmoid")
  
  model <- keras_model(
    inputs = inputs,
    outputs = classify
  )
  
  # model %>% compile(
  #   optimizer = optimizer_rmsprop(lr = 0.00001),
  #   loss = bce_dice_loss,
  #   metrics = custom_metric("dice_coef", dice_coef)
  # )
  
  return(model)
}

# FUNCTION: pred_unet - Run Tensorflow model ------------------------------------------------------------
pred_unet <- function(x,y,z) { #X is scratch diretory and Y is number of threads, z is keras batch size
  require(keras)
  require(terra)
  require(abind)
  #require(parallel)
  #require(doParallel)
  #require(foreach)
  #require(spatial.tools)
  require(TileManager)
  require(gdalUtils)
  require(rgdal)
  require(sp)
  
  use_condaenv("r-tensorflow")
  
  ### Dont need this for prediction. Make sure loss function is not 
  ### defined in model()
  # # Loss Functions for U-net Architecture
  # # Dont need for predict().
  K <- backend()
  
  # dice_coef <- function(y_true, y_pred, smooth = 1.0) {
  #   y_true_f <- k_flatten(y_true)
  #   y_pred_f <- k_flatten(y_pred)
  #   intersection <- k_sum(y_true_f * y_pred_f)
  #   result <- (2 * intersection + smooth) /
  #     (k_sum(y_true_f) + k_sum(y_pred_f) + smooth)
  #   return(result)
  # }
  # 
  # bce_dice_loss <- function(y_true, y_pred) {
  #   result <- loss_binary_crossentropy(y_true, y_pred) +
  #     (1 - dice_coef(y_true, y_pred))
  #   return(result)
  # }
  
  setwd(x)
  setwd("tiles")
  
  pred_samples <- length(dir(pattern="[0-9]{1,2}.tif$"))
  pred_names <- dir(pattern="[0-9]{1,2}.tif$")
  pred_id <- gsub(".tif$","",pred_names)
  
  print("Loading TensorFlow Models...")
  sage_model <- get_unet_1024()
  load_model_weights_hdf5(sage_model, paste(binDir,"/KerasSagePred_Unet.hd5", sep=""))  
  bare_model <- get_unet_1024()
  load_model_weights_hdf5(bare_model, paste(binDir,"/KerasBarePred_Unet.hd5", sep=""))  
  
  # Load the tile polygons for clipping the prediction rasters 
  print("Loading Tile Polygons...")
  load("../TilePolys.RData")
  
  for (j in seq_along(pred_names)) {
    
    print(paste("Processing tile ",pred_names[j]))
    t.img <- rast(pred_names[j])
    rgb <- subset(t.img,1:3)
    rgb <- rgb/255
    rgb1 <- as.matrix(rgb[[1]], wide=T)
    rgb2 <- as.matrix(rgb[[2]], wide=T)
    rgb3 <- as.matrix(rgb[[3]], wide=T)
    
    chm <-t.img[[4]]
    chm <- as.matrix(chm, wide=T)
    
    t.array <- abind(rgb1,rgb2,rgb3,chm,along=3)
    pred.array <- array(as.numeric(),dim=c(1,1024,1024,4))
    pred.array[1,,,] <- t.array
    
    # Run the TensorFlow model
    sage.out <- predict(sage_model,pred.array,batch_size = z)
    bare.out <- predict(bare_model,pred.array,batch_size = z)
    
    # Trim the result to the tile boundary
    t.extent <- bbox(buff.poly[[pred_id[j]]])
    t.crs <- crs(buff.poly[[pred_id[j]]])
    sage.rast <- terra::rast(sage.out[1,,,])
    ext(sage.rast) <- c(t.extent[1,1],t.extent[1,2],t.extent[2,1],t.extent[2,2])
    crs(sage.rast) <- t.crs
    bare.rast <- terra::rast(bare.out[1,,,])
    ext(bare.rast) <- c(t.extent[1,1],t.extent[1,2],t.extent[2,1],t.extent[2,2])
    crs(bare.rast) <- t.crs
    t.extent <- bbox(pred.poly[[pred_id[j]]])
    sage.crop <- terra::crop(sage.rast,ext(c(t.extent[1,1],t.extent[1,2],t.extent[2,1],t.extent[2,2])))
    bare.crop <- terra::crop(bare.rast,ext(c(t.extent[1,1],t.extent[1,2],t.extent[2,1],t.extent[2,2])))
    
    # Write out the model predictions
    terra::writeRaster(sage.rast, filename = paste(x,"/tiles/",pred_id[j],"_sage.tif",sep=""),
                       filetype = "GTiff", overwrite=TRUE, datatype="FLT4S")
    terra::writeRaster(bare.rast, filename = paste(x,"/tiles/",pred_id[j],"_bare.tif",sep=""),
                       filetype = "GTiff", overwrite=TRUE, datatype="FLT4S")
    
  }
  gc()
  
  # Reasemble the tiles into the full raster models. Trim to AOI
  
  print("Reassembling the tiles and trimming to the AOI...")
  
  # Build the lists of sage and bare ground tiles and assemble into a virtual raster
  sage_tiles <- length(dir(pattern="sage.tif$"))
  sage_list <- dir(pattern="sage.tif$")
  sage_vrt <- terra::vrt(sage_list)
  
  bare_tiles <- length(dir(pattern="bare.tif$"))
  bare_list <- dir(pattern="bare.tif$")
  bare_vrt <- terra::vrt(bare_list)
  
  # Crop the rasters to the AOI
  aoi <- terra::vect(paste("../aoi_proj.shp",sep="/"),layer="aoi_proj")
  sage_vrt <- terra::mask(sage_vrt, aoi)
  bare_vrt <- terra::mask(bare_vrt, aoi)
  
  # Export the prediction rasters
  terra::writeRaster(sage_vrt, filename = "../sage_pred.tif",
                     filetype = "GTiff", overwrite=TRUE, datatype="FLT4S")
  terra::writeRaster(bare_vrt, filename = "../bg_pred.tif",
                     filetype = "GTiff", overwrite=TRUE, datatype="FLT4S")
  
  gc()
  print("Finally... Finished running the @!#%@!#$@ TensorFlow models!")
}

# FUNCTION: Clean/simplify the tensorflow models ---------------------------------------------------
clean_models <- function(x) {  #x is the project directory
  require(terra)
  require(fasterize)
  require(dplyr)
  require(sf)
  
  rast.list <- c("sage","bg")
  for (r in rast.list) {
    # load the raster
    model.rast <- terra::rast(paste(x,"/scratch/",r,"_pred.tif",sep=""))
    # convert prediction probabilities to binary
    mod1 <- Con(model.rast>0.5,1,0)
    # majority filtering to reduce singletons
    mod2 <- terra::focal(mod1,w=5,fun="modal")
    mod2[mod2==0] <- NA # Change zeros to nodata
    # region group and remove any region with area < 2m^2
    mod3 = terra::as.polygons(mod2[[1]]) |>
      sf::st_as_sf() |>
      sf::st_cast("POLYGON", warn = FALSE) |>
      dplyr::mutate(id = dplyr::row_number())
    mod3$area <- sf::st_area(mod3)
    mod4 <- mod3[mod3$area>1,]
    if(nrow(mod4)>0) {
      mod5 = fasterize::fasterize(sf = mod4, raster = raster::raster(mod2[[1]])) |>
        terra::rast(mod5)
      # Write out the cleaned raster as the class
    } else {
      mod5 <- terra::rast(extent=ext(model.rast), resolution=res(model.rast),vals=NA)
    }
    terra::writeRaster(mod5,paste(x,"/scratch/",r,"_class.tif",sep=""),overwrite=T)
    
    rm(list=c('mod1','mod2','mod3','mod4','mod5'))
  }
  gc()
}

# FUNCTION: Update the canopy height model base don the sage and bare ground model classes -----------------
updateCHM <- function(x, name) { # x is project.dir, name is project.name
  require(terra)
  require(gdalUtils)
  
  bare <- terra::rast("scratch/bg_class.tif")
  chm <- terra::rast("scratch/chm_align.tif")
  
  # crop the bare class raster to the aligned chm to ensure same extents
  bare <- terra::crop(bare, chm)
  
  # 1. extract canopy height values for cells classified as bare ground
  bare.chm <- terra::mask(chm, bare)
  
  # 2. Set canopy height max value to 0.015 for any bare ground cell with height > 0.015m
  bare.mod <- Con(bare.chm>0.15,0.15,bare.chm)
  
  # 3. Merge the updated values back into the CHM
  bare.mod2 <- terra::merge(bare.mod, chm)
  
  terra::writeRaster(bare.mod2,paste(x,"/scratch/","chm_update.tif",sep=""),overwrite=T)
}

# FUNCTION: Generate the cover class predictions ---------------------------------------------------------
gen_cover <- function(x) { # x is dir
  #Merge sage and bareground classification rasters.
  # 0 is no class, 
  # 1 is bareground
  # 2 is sagebrush
  
  require(terra)
  require(gdalUtils)
  
  setwd(x)
  
  bare <- rast("scratch/bg_class.tif")
  sage <- rast("scratch/sage_class.tif")
  aoi <- terra::vect(paste("scratch/aoi_proj.shp",sep="/"),layer="aoi_proj")
  aoi.rast <- rast(extent=ext(sage),resolution=res(sage))
  aoi.rast <- rasterize(aoi, aoi.rast, field=0)
  writeRaster(aoi.rast, "scratch/aoi.tif", datatype = "INT1U", filetype="GTiff", overwrite=TRUE)
  
  sage <- terra::classify(sage, 
                          rcl = matrix(c(0,0,1,2), nrow=2,byrow = TRUE))
  rast.class <- terra::merge(bare,sage,aoi.rast)
  #  writeRaster(rast.class,"scratch/prelim_class.tif", datatype = "INT1U", filetype="GTiff", overwrite=TRUE)
  #  align_rasters("scratch/prelim_class.tif","scratch/aoi.tif","scratch/prelim_class2.tif")
  #  rast.class <- raster("scratch/prelim_class2.tif")
  rast.class <- mask(rast.class,aoi)
  writeRaster(rast.class,"products/cover_class.tif", datatype = "INT1U", filetype="GTiff", overwrite=TRUE)
}

# FUNCTION: Generate cover statistics. -----------------------------------------------------------------
get_class_stats <- function(x) { # x is project.dir
  require(terra)
  setwd(x)
  rast.class <- rast("products/cover_class.tif")
  class.stats <- as.data.frame(freq(rast.class))
  class.stats <- class.stats[!is.na(class.stats$value),]
  class.stats$area_m2 <- class.stats$count * res(rast.class)[1] * res(rast.class)[2]
  total_area <- sum(class.stats$area_m2)
  class.stats$area_pct <- round(class.stats$area_m2 / total_area * 100,1)
  return(class.stats)
}


# FUNCTION: Estimate grass and sage heights -----------------------------------------------------------------
grass_ht <- function(x) { # x = project.dir
  require(terra)
  setwd(x)
  
  # Load CHM and sage_class rasters
  chm <- terra::rast("scratch/chm_align.tif")
  sage_class <- terra::rast("products/cover_class.tif")
  crs(sage_class) <- crs(chm)
  sage_class2 <- terra::crop(sage_class, chm) # Crop to enforce same extents
  
  # 1. find heights less than 0.5m
  con1 <- Con(chm>0.5,0,chm)
  
  # 2. expand sage classes and find areas not sagebrush
  sage_class2 <- Con(sage_class2==2,1,0)
  expand <- terra::focal(sage_class2,w=3,fun="max")
  con2 <- Con(expand==1,0,1)
  
  # 3. Get areas less than 0.5m high and not sagebrush
  times <- con1 * con2
  
  # 4. calc max height within a circular area around each cell
  grass_ht <- focal(times, w=11, fun="max")
  
  # 5. multiply result times 1.5 ??why?? and return result as grass height.
  grass_ht <- grass_ht * 1.5001
  
  # 6. Write out result
  writeRaster(grass_ht,"scratch/grass_ht.tif", datatype = "FLT8S", filetype="GTiff", overwrite=TRUE)
  
  # 7. Get height for sagebrush areas
  sage_ht <- expand * chm
  
  # 8. Write out sagebrush height result    
  writeRaster(sage_ht,"scratch/sage_ht.tif", datatype = "FLT8S", filetype="GTiff", overwrite=TRUE)
}

# FUNCTION: get_grass_stats - Creates plot/PDF of grass heights --------------------------------------
get_grass_stats <- function(x) {
  require(terra)
  
  setwd(x)
  grass <- rast("scratch/grass_ht.tif")
  grass.v <- as.vector(grass)
  grass.v <- grass.v[!is.na(grass.v)]
  agg.factor <- 1/res(grass)[1]
  grass.1m <- aggregate(grass,fact = agg.factor,fun = mean)
  grass.1v <- as.vector(grass.1m)
  grass.1v <- grass.1v[!is.na(grass.1v)]
  
  writeRaster(grass.1m, "products/grassht_est.tif", filetype="GTiff", overwrite=TRUE)
  
  pdf("reports/grass_height.pdf",width = 6,height = 4)
  layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(3,8))
  
  par(mar=c(0, 5, 1.1, 5))
  boxplot(grass.v , horizontal=TRUE , outline = FALSE , 
          ylim=c(0,max(grass.v)), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F)
  
  par(mar=c(4, 5, 0, 5))
  h <- hist(grass.v , breaks=40 , col="lightblue", border="black" , main="" , xlab="Grass height (meters)", xlim=c(0,max(grass.v)), axes = FALSE)
  axis(2)
  axis(1,seq(0,1,by=0.1))
  mtext("Grass Height", side = 3, outer= TRUE, padj = 3)
  
  par(new = T)
  
  ec <- ecdf(grass.v)
  plot(x = h$mids, y=ec(h$mids)*max(h$counts), col = rgb(0,0,0,alpha=0), axes=F, xlab=NA, ylab=NA)
  lines(x = h$mids, y=ec(h$mids)*max(h$counts), col ='red', lwd=1.5)
  axis(4, at=seq(from = 0, to = max(h$counts), length.out = 11), labels=seq(0, 1, 0.1), col = 'red', col.axis = 'red')
  mtext(side = 4, line = 3, 'Cumulative Density', col = 'red')
  
  dev.off()
  grass.dist <- summary(grass.v)
  return(grass.dist)
}

# FUNCTION: Count Sagebrush Plants ---------------------------------------------
count_sage <- function(x) { # x is project.dir
  require(terra)
  setwd(x)
  
  sage_chm <- terra::rast("scratch/sage_ht.tif")
  sage_chm <- terra::subst(sage_chm, 0, NA)
  
  sage_neg <- sage_chm * -1
  inv_chm <- sage_neg + 10
  writeRaster(sage_ht,"scratch/inv_chm.tif", filetype="GTiff", overwrite=TRUE)
  whitebox::wbt_d8_flow_accumulation("scratch/inv_chm.tif", "scratch/sage_flow.tif")  
  flow <- rast("scratch/sage_flow.tif")
  flow_focal <- focal(flow,w=51,fun='max')
}



get_sagebrush_plot <- function(x) {
  
  sage_ht <- class.list$sage.ht
  
  setwd(x)
  pdf("reports/sagebrush_height.pdf",width = 6,height = 4)
  layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(3,8))
  
  par(mar=c(0, 5, 1.1, 5))
  boxplot(sage_ht , horizontal=TRUE , outline = FALSE , 
          ylim=c(0,max(sage_ht)), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F)
  
  par(mar=c(4, 5, 0, 5))
  h <- hist(sage_ht , breaks=40 , col="lightblue", border="black" , main="" , xlab="Sagebrush height (meters)", xlim=c(0,max(sage_ht)), axes = FALSE)
  axis(2)
  axis(1,seq(0,1,by=0.1))
  mtext("Sagebrush Height", side = 3, outer= TRUE, padj = 3)
  
  par(new = T)
  
  ec <- ecdf(sage_ht)
  plot(x = h$mids, y=ec(h$mids)*max(h$counts), col = rgb(0,0,0,alpha=0), axes=F, xlab=NA, ylab=NA)
  lines(x = h$mids, y=ec(h$mids)*max(h$counts), col ='red', lwd=1.5)
  axis(4, at=seq(from = 0, to = max(h$counts), length.out = 11), labels=seq(0, 1, 0.1), col = 'red', col.axis = 'red')
  mtext(side = 4, line = 3, 'Cumulative Density', col = 'red')
  
  dev.off()
}

# FUNCTION: Estimate standing biomass --------------------------------------------------------------
plot_biomass <- function(x) { # x is project.dir
  require(terra)
  require(gdalUtils)
  require(gam)
  setwd(x)  
  
  #rescale aoi and sage and bareground
  #align_rasters("products/aoi.tif","scratch/chm2.tif","scratch/aoi.tif", overwrite=TRUE)
  #align_rasters("scratch/sage_class.tif","scratch/chm2.tif","scratch/sage_class2.tif",overwrite=TRUE)
  
  #load rescaled aoi and sage and bareground
  aoi <- terra::rast("scratch/aoi.tif")
  cover <- terra::rast("products/cover_class.tif")
  
  #Load cHM
  chm <- terra::rast("scratch/chm_update.tif")
  
  #reclassify sage for multiplication
  cover <- classify(cover,
                    matrix(c(0,1,1,1,2,0),byrow = TRUE,ncol=2))
  
  #Set max height at 0.5meter , deals with greasewood or other large shrubs/trees.
  cover.trim <- crop(cover,chm)
  crs(cover.trim) <- crs(chm)
  chm.herb <- cover.trim * chm
  chm.herb[chm.herb > 0.5] <- 0.5
  
  #Figure out how many cells to aggregate across
  agg.factor <- 1/mean(res(chm.herb))
  
  #aggregate rasters to 1m
  chm.herb.1m <- aggregate(chm.herb,agg.factor)
  aoi.1m <- aggregate(aoi,agg.factor)
  
  
  #Convert rasters to vectors for analysis
  chm.herb.v <- as.vector(chm.herb.1m)
  aoi.v <- as.vector(aoi.1m)
  aoi.idx <- which(aoi.v == 0)
  
  biomass_pred <- function(x) {
    1496.867 * x + 8.032
  }
  
  #simulate variance for total biomass estimate
  biomass_rnorm <- function(x) {
    rnorm(n=1e4,x,sd=48.63)
  }
  
  chm.herb.v.test <- chm.herb.v[aoi.idx]
  chm.herb.v.test <- unlist(lapply(chm.herb.v.test,biomass_pred))
  chm.herb.v.test[chm.herb.v.test < 0] <- 0
  
  chm.biomass.all.sum <- lapply(chm.herb.v.test,biomass_rnorm)
  chm.biomass.all.sum <- do.call(rbind,chm.biomass.all.sum)
  chm.biomass.all.sum[chm.biomass.all.sum < 0] <- 0
  chm.biomass.all.sim.sum <- colSums(chm.biomass.all.sum, na.rm = TRUE)
  chm.biomass.all.sim.sum.kg.ha <- chm.biomass.all.sim.sum / (length(chm.herb.v.test) * 0.0001) * 1e-3
  biomass.prob.density.kg.ha <- as.numeric(round(quantile(chm.biomass.all.sim.sum.kg.ha, probs=c(0,0.5,1))))
  biomass.prob.density.lb.ac <- round(0.892179 * biomass.prob.density.kg.ha)
  rm(chm.biomass.all.sum) # cleanup
  gc()
  
  #grass biomass function
  gbiomass_pred <- function(x) {
    1628.61 * x - 17.07
  }
  
  #simulate variance for grass biomass estimate
  gbiomass_rnorm <- function(x) {
    rnorm(n=1e4,x,sd=49.79)
  }
  
  chm.grass.v.test <- chm.herb.v[aoi.idx]
  chm.grass.v.test <- unlist(lapply(chm.grass.v.test,gbiomass_pred))
  chm.grass.v.test[chm.grass.v.test < 0] <- 0
  
  chm.grass.all.sum <- lapply(chm.grass.v.test,gbiomass_rnorm)
  chm.grass.all.sum <- do.call(rbind,chm.grass.all.sum)
  chm.grass.all.sum[chm.grass.all.sum < 0] <- 0
  chm.grass.all.sim.sum <- colSums(chm.grass.all.sum, na.rm = TRUE)
  chm.grass.all.sim.sum.kg.ha <- chm.grass.all.sim.sum / (length(chm.grass.v.test) * 0.0001) * 1e-3
  grass.prob.density.kg.ha <- as.numeric(round(quantile(chm.grass.all.sim.sum.kg.ha, probs=c(0,0.5,1))))
  grass.prob.density.lb.ac <- round(0.892179 * grass.prob.density.kg.ha)
  rm(chm.grass.all.sum) # cleanup
  gc()
  
  #Forb biomass model
  load(path.gam.forb.model)
  
  #simulate variance for forb biomass estimate
  fbiomass_rnorm <- function(x) {
    rnorm(n=1e4,x,sd=6)
  }
  
  chm.forb.v.test <- chm.herb.v[aoi.idx]
  chm.forb.v.test <-  predict(forb_gam, newdata=data.frame(x=chm.forb.v.test))
  chm.forb.v.test[chm.forb.v.test < 0] <- 0
  
  chm.forb.all.sum <- lapply(chm.forb.v.test,fbiomass_rnorm)
  chm.forb.all.sum <- do.call(rbind,chm.forb.all.sum)
  chm.forb.all.sum[chm.forb.all.sum < 0] <- 0
  chm.forb.all.sim.sum <- colSums(chm.forb.all.sum, na.rm = TRUE)
  chm.forb.all.sim.sum.kg.ha <- chm.forb.all.sim.sum / (length(chm.forb.v.test) * 0.0001) * 1e-3
  forb.prob.density.kg.ha <- as.numeric(round(quantile(chm.forb.all.sim.sum.kg.ha, probs=c(0,0.5,1))))
  forb.prob.density.lb.ac <- round(0.892179 * forb.prob.density.kg.ha)
  rm(chm.forb.all.sum) # cleanup
  gc()
  
  biomass.list <- list()
  biomass.list[[1]] <- project.name
  biomass.list[[2]] <- biomass.prob.density.kg.ha
  biomass.list[[3]] <- biomass.prob.density.lb.ac
  biomass.list[[4]] <- grass.prob.density.kg.ha 
  biomass.list[[5]] <- grass.prob.density.lb.ac 
  biomass.list[[6]] <- forb.prob.density.kg.ha 
  biomass.list[[7]] <- forb.prob.density.lb.ac 
  
  save(biomass.list, file = "scratch/biomass.list.RData")
  plotName <- project.name
  
  #Distribution of total biomass by 1meter squares in plot
  plot_biomass_average <- function(){
    list_density <- max(density(chm.herb.v.test,na.rm=TRUE)$y) * 1.2
    list_histo <- max(hist(chm.herb.v.test,plot = FALSE)$density) * 1.2
    max.den <- max(c(list_histo,list_density))
    hist(chm.herb.v.test, # histogram
         col = "lightblue", # column color
         border = "black", 
         prob = TRUE, # show densities instead of frequencies
         xlim = c(0,500),
         ylim = c(0,max.den),
         xlab = expression("Biomass per area (grams /" ~m^2 *")"),
         main = paste("Plot ",plotName,": Herbaceous biomass (dry)",sep=""))
    lines(density(chm.herb.v.test,na.rm=TRUE), # density plot
          lwd = 2, # thickness of line
          col = "red")
    abline(v = mean(chm.herb.v.test, na.rm=TRUE),
           col = "black",
           lwd = 2)
    abline(v = median(chm.herb.v.test, na.rm=TRUE),
           col = "darkblue",
           lwd = 2)
    
    legend(x = "topright", # location of legend within plot area
           c("Density plot", 
             paste("Mean:",round(mean(chm.herb.v.test, na.rm=TRUE)),"grams"), 
             paste("Median:",round(median(chm.herb.v.test, na.rm=TRUE)),"grams")),
           col = c("red", "black","darkblue"),
           lwd = c(2, 2, 2))
  }
  
  #Estimate Total Biomass in Plot on Per Area Basis
  plot_biomass_total <- function(){
    list_density <- max(density(chm.biomass.all.sim.sum.kg.ha,na.rm=TRUE)$y) * 1.2
    list_histo <- max(hist(chm.biomass.all.sim.sum.kg.ha,plot = FALSE)$density) * 1.2
    max.den <- max(c(list_histo,list_density))
    hist(chm.biomass.all.sim.sum.kg.ha, # histogram
         col = "lightblue", # column color
         border = "black", 
         prob = TRUE, # show densities instead of frequencies
         #xlim = c(0,500),
         ylim = c(0,max.den),
         xlab = "Biomass per area (kg / ha)",
         main = paste("Plot ",plotName,": Estimated total herbaceous biomass (dry)",sep=""))
    lines(density(chm.biomass.all.sim.sum.kg.ha,na.rm=TRUE), # density plot
          lwd = 2, # thickness of line
          col = "red")
    abline(v = mean(chm.biomass.all.sim.sum.kg.ha, na.rm=TRUE),
           col = "black",
           lwd = 2)
    
    legend(x = "topright", # location of legend within plot area
           c("Density plot", 
             "Mean"
           ),
           col = c("red", "black"),
           lwd = c(2, 2))
    
    legend(x = "topleft", # location of legend within plot area
           c(paste(biomass.prob.density.kg.ha[2],
                   " kg / ha (",
                   biomass.prob.density.kg.ha[1],
                   " to ",
                   biomass.prob.density.kg.ha[3],
                   ")",sep=""),
             paste(biomass.prob.density.lb.ac[2],
                   " lbs / acre (",
                   biomass.prob.density.lb.ac[1],
                   " to ",
                   biomass.prob.density.lb.ac[3],
                   ")",sep=""))
    )
  }
  
  #Distribution of grass biomass by 1meter squares in plot
  plot_grass_average <- function(){
    list_density <- max(density(chm.grass.v.test,na.rm=TRUE)$y) * 1.2
    list_histo <- max(hist(chm.grass.v.test,plot = FALSE)$density) * 1.2
    max.den <- max(c(list_histo,list_density))
    hist(chm.grass.v.test, # histogram
         col = "lightblue", # column color
         border = "black", 
         prob = TRUE, # show densities instead of frequencies
         xlim = c(0,500),
         ylim = c(0,max.den),
         xlab = expression("Biomass per area (grams /" ~m^2 *")"),
         main = paste("Plot ",plotName,": Grass biomass (dry)",sep=""),
         yaxs='i',xaxs='i')
    lines(density(chm.grass.v.test,na.rm=TRUE), # density plot
          lwd = 2, # thickness of line
          col = "red")
    abline(v = mean(chm.grass.v.test, na.rm=TRUE),
           col = "black",
           lwd = 2)
    abline(v = median(chm.grass.v.test, na.rm=TRUE),
           col = "darkblue",
           lwd = 2)
    
    legend(x = "topright", # location of legend within plot area
           c("Density plot", 
             paste("Mean:",round(mean(chm.grass.v.test, na.rm=TRUE)),"grams"), 
             paste("Median:",round(median(chm.grass.v.test, na.rm=TRUE)),"grams")),
           col = c("red", "black","darkblue"),
           lwd = c(2, 2, 2))
  }
  
  #Estimate grass Biomass in Plot on Per Area Basis
  plot_grass_total <- function(){
    list_density <- max(density(chm.grass.all.sim.sum.kg.ha,na.rm=TRUE)$y) * 1.2
    list_histo <- max(hist(chm.grass.all.sim.sum.kg.ha,plot = FALSE)$density) * 1.2
    max.den <- max(c(list_histo,list_density))
    hist(chm.grass.all.sim.sum.kg.ha, # histogram
         col = "lightblue", # column color
         border = "black", 
         prob = TRUE, # show densities instead of frequencies
         #xlim = c(0,500),
         ylim = c(0,max.den),
         xlab = "Biomass per area (kg / ha)",
         main = paste("Plot ",plotName,": Estimated total grass biomass (dry)",sep=""))
    lines(density(chm.grass.all.sim.sum.kg.ha,na.rm=TRUE), # density plot
          lwd = 2, # thickness of line
          col = "red")
    abline(v = mean(chm.grass.all.sim.sum.kg.ha, na.rm=TRUE),
           col = "black",
           lwd = 2)
    
    legend(x = "topright", # location of legend within plot area
           c("Density plot", 
             "Mean"
           ),
           col = c("red", "black"),
           lwd = c(2, 2))
    
    legend(x = "topleft", # location of legend within plot area
           c(paste(grass.prob.density.kg.ha[2],
                   " kg / ha (",
                   grass.prob.density.kg.ha[1],
                   " to ",
                   grass.prob.density.kg.ha[3],
                   ")",sep=""),
             paste(grass.prob.density.lb.ac[2],
                   " lbs / acre (",
                   grass.prob.density.lb.ac[1],
                   " to ",
                   grass.prob.density.lb.ac[3],
                   ")",sep=""))
    )
  }
  
  #Distribution of forb biomass by 1meter squares in plot
  plot_forb_average <- function(){
    list_density <- max(density(chm.forb.v.test,na.rm=TRUE)$y) * 1.2
    list_histo <- max(hist(chm.forb.v.test,plot = FALSE)$density) * 1.2
    max.den <- max(c(list_histo,list_density))
    hist(chm.forb.v.test, # histogram
         col = "lightblue", # column color
         border = "black", 
         prob = TRUE, # show densities instead of frequencies
         xlim = c(0,50),
         ylim = c(0,max.den),
         xlab = expression("Biomass per area (grams /" ~m^2 *")"),
         main = paste("Plot ",plotName,": Forb biomass (dry)", sep=""),
         yaxs='i',xaxs='i')
    lines(density(chm.forb.v.test,na.rm=TRUE), # density plot
          lwd = 2, # thickness of line
          col = "red")
    abline(v = mean(chm.forb.v.test, na.rm=TRUE),
           col = "black",
           lwd = 2)
    abline(v = median(chm.forb.v.test, na.rm=TRUE),
           col = "darkblue",
           lwd = 2)
    
    legend(x = "topright", # location of legend within plot area
           c("Density plot", 
             paste("Mean:",round(mean(chm.forb.v.test, na.rm=TRUE)),"grams"), 
             paste("Median:",round(median(chm.forb.v.test, na.rm=TRUE)),"grams")),
           col = c("red", "black","darkblue"),
           lwd = c(2, 2, 2))
  }
  
  #Estimate forb Biomass in Plot on Per Area Basis
  plot_forb_total <- function(){
    list_density <- max(density(chm.forb.all.sim.sum.kg.ha,na.rm=TRUE)$y) * 1.2
    list_histo <- max(hist(chm.forb.all.sim.sum.kg.ha,plot = FALSE)$density) * 1.2
    max.den <- max(c(list_histo,list_density))
    hist(chm.forb.all.sim.sum.kg.ha, # histogram
         col = "lightblue", # column color
         border = "black", 
         prob = TRUE, # show densities instead of frequencies
         #xlim = c(0,500),
         ylim = c(0,max.den),
         xlab = "Biomass per area (kg / ha)",
         main = paste("Plot ",plotName,": Estimated total forb biomass (dry)",sep=""))
    lines(density(chm.forb.all.sim.sum.kg.ha,na.rm=TRUE), # density plot
          lwd = 2, # thickness of line
          col = "red")
    abline(v = mean(chm.forb.all.sim.sum.kg.ha, na.rm=TRUE),
           col = "black",
           lwd = 2)
    
    legend(x = "topright", # location of legend within plot area
           c("Density plot", 
             "Mean"
           ),
           col = c("red", "black"),
           lwd = c(2, 2))
    
    legend(x = "topleft", # location of legend within plot area
           c(paste(forb.prob.density.kg.ha[2],
                   " kg / ha (",
                   forb.prob.density.kg.ha[1],
                   " to ",
                   forb.prob.density.kg.ha[3],
                   ")",sep=""),
             paste(forb.prob.density.lb.ac[2],
                   " lbs / acre (",
                   forb.prob.density.lb.ac[1],
                   " to ",
                   forb.prob.density.lb.ac[3],
                   ")",sep=""))
    )
  }
  
  
  pdf("reports/biomass_est.pdf",width = 6.5,height = 10)
  par(mfrow=c(3,1))
  plot_biomass_total()
  plot_grass_total()
  plot_forb_total()
  plot_biomass_average()
  plot_grass_average()
  plot_forb_average()
  dev.off()
  
  #chm.herb.v <- rep(NA,length(aoi.v))
  #chm.herb.v[aoi.idx] <- chm.herb.v.test
  #chm.biomass <- raster(as.matrix(chm.herb.v, nrow=nrow(chm.herb.1m), ncol=ncol(chm.herb.1m), by.row=TRUE),template=chm.herb.1m)
  #writeRaster(chm.biomass,"products/biomass_est.tif",format = "GTiff",overwrite=TRUE)
  
  return(biomass.list)
}


# FUNCTION: Plot Sagebrush Heights -----------------------------------------------------------------
get_sagebrush_plot <- function(x) {
  
  sage_ht <- class.list$sage.ht
  
  setwd(x)
  pdf("reports/sagebrush_height.pdf",width = 6,height = 4)
  layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(3,8))
  
  par(mar=c(0, 5, 1.1, 5))
  boxplot(sage_ht , horizontal=TRUE , outline = FALSE , 
          ylim=c(0,max(sage_ht)), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F)
  
  par(mar=c(4, 5, 0, 5))
  h <- hist(sage_ht , breaks=40 , col="lightblue", border="black" , main="" , xlab="Sagebrush height (meters)", xlim=c(0,max(sage_ht)), axes = FALSE)
  axis(2)
  axis(1,seq(0,1,by=0.1))
  mtext("Sagebrush Height", side = 3, outer= TRUE, padj = 3)
  
  par(new = T)
  
  ec <- ecdf(sage_ht)
  plot(x = h$mids, y=ec(h$mids)*max(h$counts), col = rgb(0,0,0,alpha=0), axes=F, xlab=NA, ylab=NA)
  lines(x = h$mids, y=ec(h$mids)*max(h$counts), col ='red', lwd=1.5)
  axis(4, at=seq(from = 0, to = max(h$counts), length.out = 11), labels=seq(0, 1, 0.1), col = 'red', col.axis = 'red')
  mtext(side = 4, line = 3, 'Cumulative Density', col = 'red')
  
  dev.off()
}

# FUNCTION: Generate final output table ------------------------------------------------------------
gen_table_data <- function(x) {
  setwd(x)
  class.list$class.stats
  summary.table <- data.frame(
    Metric = c("Sagebrush Cover","Bareground Cover","Unclassified Cover","Herb. Biomass","Herb. Biomass","Grass Biomass","Grass Biomass","Forb Biomass","Forb Biomass","Mean Grass Height","Mean Grass Height"),
    Units = c("Percent","Percent","Percent","kg/ha","lbs/acres","kg/ha","lbs/acre","kg/ha","lbs/acre","meters","inches"),
    Result = c(class.list$class.stats$area_pct[3],class.list$class.stats$area_pct[2],class.list$class.stats$area_pct[1],
               biomass.stats[[2]][2],biomass.stats[[3]][2],biomass.stats[[4]][2],biomass.stats[[5]][2],biomass.stats[[6]][2],biomass.stats[[7]][2],grass.stats[4],grass.stats[4]* 39.3701)
  )
  
  summary.table[1:10,3] <- round(summary.table[1:10,3],0)
  summary.table[11:12,3] <- round(summary.table[11:12,3],2)
  
  names(summary.table)[3] <- biomass.stats[[1]]
  write.csv(summary.table,"reports/plot_summary.csv",row.names = FALSE)
}

# FUNCTION: Clean up scratch datasets ---------------------------------------------------------------
clean_scratch <- function(x, layers=TRUE, tiles=FALSE) {  # x is scratch directory
  setwd(x)
  scratch.files <- c("aoi.tif","bg_class.tif","bg_pred.tif","chm_align.tif","chm_update.tif",
                     "grass_ht.tif","sage_class.tif",
                     "aoi_proj.cpg","aoi_proj.dbf","aoi_proj.prj","aoi_proj.shp","aoi_proj.shx",
                     "tilePolys.RData")
  
  if (layers) {
    for (s.file in scratch.files) {
      if (file.exists(s.file)) {
        file.remove(s.file)
      }
    }
  }
    
  if (tiles) {
    tile.tifs <- list.files(file.path(x,"tiles"),pattern="*.tif")
    for (tif in tile.tifs) {
      if (file.exists(file.path(x,"tiles",tif))) {
        file.remove(file.path(x,"tiles",tif))
      }
    }
  }
}


# FUNCTION: remove left over temp files --------------------------------------------------------------
clean_temp <- function(x) { # x is tempdir()
 files <- list.files(x, pattern="*.tif")
 for (f in files) {
   if (file.exists(file.path(x,f))) {
     file.remove(file.path(x,f))
   }
 } 
}