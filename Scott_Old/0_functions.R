#Checks to see if all required lastools executibles exist.
check_lastools <- function(x) {
  tool.list<- c("lastile.exe",
                "lassort.exe",
                "lasthin.exe",
                "lasnoise.exe",
                "lasground.exe",
                "lasheight.exe",
                "las2las.exe",
                "lasmerge.exe",
                "las2shp.exe")
  tool <- paste(x,tool.list,sep="")
  result <- file.exists(tool)
  result <- data.frame(tool,result,stringsAsFactors = FALSE)
  result <- result[result$result == FALSE,1]
  return(result)
}

# AOI Polygon for clipping and alignment. Use camera coordinates to build rectangular 
# bounding box. 
generate_aoi <- function(x) { # X is project.scratch
  require(maptools)
  require(sp)
  require(rgeos)
  require(rgdal)
  require(shotGroups)
  
  setwd(x)
  
  #Target EPSG for Montana StatePlane CRS
  CRS.new <- CRS("+init=epsg:32100") 
  
  #move into photogrammetry directory and process
  setwd("photogrammetry/")
  
  #Load projections data. Provide Projection error (in pixels) and tiepoints.
  p.df <- read.table('projections.csv',sep='\t', stringsAsFactors = FALSE)
  
  #Load the reference data for camera locations.
  ref.df <- read.csv("reference.csv", skip=1, stringsAsFactors = FALSE)
  
  #remove file suffix to allow for table merge.
  ref.df$X.Merge <- gsub(".JPG","",ref.df$X.Label)
  
  #merge tables.
  t.df <- merge(ref.df,p.df,by.x = "X.Merge",by.y= "V1")
  
  #remove un-needed data 
  t.df <- subset(t.df,select=c(X.Label, #Image/Camera name
                               X.Longitude, #Original Longitude
                               Y.Latitude,  #Original Latitude
                               Z.Altitude,  #Original Alt (m)
                               X_est,  #Calculated Long for Camera 
                               Y_est,  #Calculated Lat for Camera
                               Z_est,  #Calculate Alt. for Camera
                               V2, #Projection Error (XY) in pixels
                               V3)) #Number of tie points per camera
  names(t.df) <- c("camera","oLong","oLat","oAlt","fLong","fLat","fAlt","projErr","TiePts")
  #Print summary projection errors
  print("Projection Errors")
  print(summary(t.df$projErr))
  
  #Pause if projection errors are high
  if(max(t.df$projErr > 0.35)){
    print("WARNING: projection error exceeds 0.35")
    #invisible(readline(prompt="Press [enter] to continue"))
  }
  
  #Print summary of tiepoints
  print("Camera Tiepoints")
  print(summary(t.df$TiePts))
  
  #Pause if tiepoints were less than 500
  if(min(t.df$TiePts < 500)){
    print("WARNING: Low Tiepoints")
    #invisible(readline(prompt="Press [enter] to continue"))
  }
  
  #Next we convert points into the output coordinate system
  #And calculate the minimum bounding box that will be used
  #later to define the AOI
  #The target coordinate system is defined as CRS.new and 
  #is defined at the top of this script.
  
  #Extract lat/long coordinates and assign to spatial object (d)
  d <- subset(t.df, select=c("fLong","fLat"))
  coordinates(d) <- c("fLong", "fLat")
  proj4string(d) <- CRS("+init=epsg:4326") # WGS 84
  
  #Reproject into desired coordiante system (CRS.new)
  d.new <- spTransform(d,CRS.new)
  
  #Add the coordinates of cameras in the final coordinate system to the 
  #Reference Table
  t.df$Xsp <- d.new$fLong
  t.df$Ysp <- d.new$fLat
  
  #Get minimum bounding box for later AOI selection. 
  #This creates a 4-side polygon of arbitrary rotation to clip 
  #the models to.
  
  xy <- as.matrix(data.frame(x = t.df$Xsp, y = t.df$Ysp))
  bb <- getMinBBox(xy) #from shotGroups package
  
  write.csv(t.df,"cameras.csv", row.names = FALSE) # Writes out table to Photgrammetry directory
  if(file.exists("cameras.csv")){
    file.rename("projections.csv","../scratch/projections.csv")
    file.rename("reference.csv","../scratch/reference.csv")
  }
  
  #Change directory to scratch to export spatial data for
  #1. camera locations 
  #2. AOI bounding box.
  setwd("../scratch")
  
  #Joint the cameras dataframe to the spatial data to export together 
  #as shapefile
  t.sd <- SpatialPointsDataFrame(d.new,t.df)
  writeOGR(t.sd,getwd(),"Cameras",driver = "ESRI Shapefile")
  
  #export AOI Poly as shapefile
  aoiPoly <- Polygon(bb$pts) #bb points to polygon class
  aoiPoly <- Polygons(list(aoiPoly),1) #wrap polygon in a list
  aoiPoly <- SpatialPolygons(list(aoiPoly), proj4string = CRS.new) #update class as Spatial
  bb[[1]] <- NULL #get rid of vertix coordinates
  bb <- as.data.frame(t(unlist(bb))) #force polygon stats into table
  aoiPoly <- SpatialPolygonsDataFrame(aoiPoly,bb,match.ID= FALSE)
  writeOGR(aoiPoly,getwd(),"aoiPoly", driver = "ESRI Shapefile",overwrite_layer = TRUE)
}

# Run a quick check to ensure that required files from photogrammetry 
# and pointcloud processing exist.
check_prelim_products <- function(x) { # x is project.scratch
  setwd(x)
  setwd("products")
  product.list <- c("ortho1.tif",
                "ortho2.tif",
                "dsm.tif",
                "dem.tif",
                "chm.tif",
                "aoi.tif",
                "pointcloud.laz",
                "pointcloud_ground.laz"
                )
  products <- paste(x,"/products/",product.list,sep="")
  result <- file.exists(products)
  result <- data.frame(products,result,stringsAsFactors = FALSE)
  result <- result[result$result == FALSE,1]
  return(result)
}

# Tile data to pass to U-Net CNN. 
tile_data <- function(x,y) { # Tiles RGB-CHM into 1024x1024x4 (RGBZ) tiles to pass to keras
  require(raster)
  require(gdalUtils)
  require(spatial.tools)
  require(TileManager)
  require(doParallel)  
  require(foreach) 
  require(rgdal)
  require(sp)
  
  setwd(x)
  nCores <- y
  
  ortho <- brick("products/ortho1.tif")
  ortho <- readAll(dropLayer(ortho,4))
  
  #TensorFlow does not like NA values.
  ortho[is.na(ortho[])] <- 0 
  
  #realign chm. Sometimes extent is a little wonky coming out of Arc.
  align_rasters("products/chm.tif","products/ortho1.tif","scratch/chm_align.tif", nThreads = "ALL_CPUS")
  
  chm <- readAll(raster("scratch/chm_align.tif"))
  
  #rescale chm 0 to 1 with max height (1) equal 2m.
  chm <- readAll(chm  / 2)
  chm[chm > 1] <- 1
  
  ## test add to address mislabelling at edge
  chm[is.na(chm[])] <- 0
  
  #combine rasters into a rgbz brick.
  combined_raster <- brick(ortho[[1]],ortho[[2]],ortho[[3]],chm)
  rm(chm,ortho)
  unlink("scratch/chm_align.tif", force = TRUE)
  gc()
  combined_raster <- readAll(combined_raster)
  gc()
  combined_raster <- modify_raster_margins(combined_raster, extent_delta = c(0,2000,2000,0), value = NA)
  gc()
  
  ts.ortho <- TileScheme(combined_raster,
                         dimByCell=c(1000,1000),
                         buffer = 12,
                         bufferspill = FALSE,
                         removeEmpty = TRUE)
  
  buff.poly <- ts.ortho$buffPolygons
  pred.poly <- ts.ortho$tilePolygons
  n.tiles <- dim(buff.poly)[1]
  poly.id <- as.character(buff.poly[[1]])
  
  if(dir.exists("scratch/tiles")) {
    unlink("scratch/tiles", recursive=TRUE, force=TRUE)
  }
  
  dir.create("scratch/tiles")
  save(buff.poly,pred.poly, file = "scratch/TilePolys.RData")
  
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  foreach(j=1:n.tiles) %dopar% {
    require(raster)
    
    t.extent <- extent(buff.poly[j,])
    
    t.rast <- crop(combined_raster,t.extent)
    writeRaster(t.rast, filename = paste("scratch/tiles/",poly.id[j],".tif",sep=""),
                format = "GTiff", overwrite=TRUE, datatype="FLT4S", options=c("PROFILE=BASELINE"))
    rm(t.extent,t.rast)
    gc()
  }
  
  stopCluster(cl)
  
  setwd("scratch/tiles")
  unlink(dir(pattern=".aux.xml$"))
}

pred_unet <- function(x,y,z) { #X is scratch diretory and Y is number of threads, z is keras batch size
  require(keras)
  require(raster)
  require(abind)
  require(parallel)
  require(doParallel)
  require(foreach)
  require(spatial.tools)
  require(TileManager)
  require(gdalUtils)
  require(rgdal)
  require(sp)
  
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
  
  # U-net 1024
  get_unet_1024 <- function(input_shape = c(1024, 1024, 4),  #input shape is 1024,1024,4 - RGB(CHM)
                            num_classes = 1) { 
    
    inputs <- layer_input(shape = input_shape)
    # 1024
    
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
    
    center <- down7_pool %>%
      layer_conv_2d(filters = 1024, kernel_size = c(5, 5), padding = "same") %>%
      layer_batch_normalization() %>%
      layer_activation("relu") %>%
      layer_conv_2d(filters = 1024, kernel_size = c(5, 5), padding = "same") %>%
      layer_batch_normalization() %>%
      layer_activation("relu") 
    # center
    
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
  
  cl <- makeCluster(y)
  
  clusterEvalQ(cl, {
    require(abind)
    require(raster)
  })
  
  registerDoParallel(cl)
  
  setwd(x)
  setwd("scratch/tiles")
  
  pred_samples <- length(dir(pattern=".tif$"))
  pred_names <- dir(pattern=".tif$")
  pred_id <- gsub(".tif$","",pred_names)

  pred.list <- list()
  
  pred.list <- foreach(j = seq_along(pred_names)) %dopar% {
    setwd(x)
    setwd("scratch/tiles")
    t.img <- brick(pred_names[j])
    rgb <- dropLayer(t.img,4)
    rgb <- rgb/255
    rgb <- as.array(rgb)
    
    chm <-t.img[[4]]
    chm <- as.matrix(chm)
    
    t.array <- abind(rgb,chm,along=3)
    return(t.array)
  }

  stopCluster(cl)
  gc()
  
  pred.array <- array(as.numeric(),dim=c(pred_samples,1024,1024,4))
  
  for(k in seq_along(pred.list)) {
    pred.array[k,,,] <- pred.list[[k]]
  }
  
  rm(pred.list)
  gc()
  
  #predict sagebrush
  model <- get_unet_1024()
  load_model_weights_hdf5(model, paste(binDir,"/KerasSagePred_Unet.hd5", sep=""))
  sage.out <- predict(model,pred.array,batch_size = z)
  gc()
  
  #predict bareground
  model <- get_unet_1024()
  load_model_weights_hdf5(model, paste(binDir,"/KerasBarePred_Unet.hd5", sep=""))
  bg.out <- predict(model,pred.array,batch_size = z)
  
  rm(pred.array)
  gc()
  setwd("../")
  
  load("TilePolys.RData")
  buff.names <- as.character(buff.poly$tileName)
  crop.names <- as.character(pred.poly$tileName)
  
  sage.list <- list()
  for(j in 1:pred_samples) {  
    buff.id <-which(pred_id[j] == buff.names)
    t.extent <- extent(buff.poly[buff.id,])
    t.crs <- crs(buff.poly[buff.id,])
    t.rast <- raster(sage.out[j,,,])
    extent(t.rast) <- t.extent
    crs(t.rast) <- t.crs
    crop.id <-which(pred_id[j] == crop.names)
    t.extent <- extent(pred.poly[crop.id,])
    t.rast <- crop(t.rast,t.extent)
    sage.list[[j]] <- t.rast
  }
  rm(sage.out)
  
  bg.list <- list()
  for(j in 1:pred_samples) {
    buff.id <-which(pred_id[j] == buff.names)
    t.extent <- extent(buff.poly[buff.id,])
    t.crs <- crs(buff.poly[buff.id,])
    t.rast <- raster(bg.out[j,,,])
    extent(t.rast) <- t.extent
    crs(t.rast) <- t.crs
    crop.id <-which(pred_id[j] == crop.names)
    t.extent <- extent(pred.poly[crop.id,])
    t.rast <- crop(t.rast,t.extent)
    bg.list[[j]] <- t.rast
  }
  rm(bg.out)

  names(sage.list)[1:2] <- c('x', 'y')
  sage.list$fun <- mean
  sage.list$na.rm <- TRUE
  mos2 <- do.call(mosaic,sage.list)
  writeRaster(mos2,filename = "sagepred.tif", format="GTiff",datatype = "INT1U", overwrite=TRUE)  
  rm(sage.list)
  
  names(bg.list)[1:2] <- c('x', 'y')
  bg.list$fun <- mean
  bg.list$na.rm <- TRUE
  mos2 <- do.call(mosaic,bg.list)
  writeRaster(mos2,filename = "bgpred.tif", format="GTiff",datatype = "INT1U", overwrite=TRUE)  
  rm(bg.list)
  gc()
}

#Merge sage and bareground classification rasters.
# 0 is no class, 
# 1 is bareground
# 2 is sagebrush

gen_cover <- function(x) { # x is project.scratch
  require(raster)
  require(gdalUtils)
  
  setwd(x)
  
  aoi <- raster("products/aoi.tif")
  bare <- raster("scratch/bg_class.tif")
  sage <- raster("scratch/sage_class.tif")
  
  sage <- reclassify(sage, 
                     rcl = matrix(c(0,NA,1,2), nrow=2,byrow = TRUE))
  rast.class <- raster::merge(sage,bare)
  writeRaster(rast.class,"scratch/prelim_class.tif", datatype = "INT1U", format="GTiff", overwrite=TRUE)
  align_rasters("scratch/prelim_class.tif","products/aoi.tif","scratch/prelim_class2.tif")
  rast.class <- raster("scratch/prelim_class2.tif")
  rast.class <- mask(rast.class,aoi)
  writeRaster(rast.class,"scratch/prelim_class.tif", datatype = "INT1U", format="GTiff", overwrite=TRUE)
  unlink("scratch/prelim_class2.tif")
}

copy_cover <- function(x) {
  setwd(x)
  flist1 <- dir("scratch",pattern = "prelim_class.",full.names = TRUE)
  flist2 <- gsub("^scratch/","products/",flist1)
  flist2 <- gsub("prelim_class.","CoverClass.",flist2)
  file.copy(flist1, flist2, overwrite = TRUE)
}

#Generate cover statistics.
get_class_stats <- function(x) {
  require(raster)
  setwd(x)
  rast.class <- raster("products/CoverClass.tif")
  class.stats <- as.data.frame(freq(rast.class, progress ="text"))
  class.stats <- class.stats[!is.na(class.stats$value),]
  class.stats$area_m2 <- class.stats$count * res(rast.class)[1] * res(rast.class)[2]
  total_area <- sum(class.stats$area_m2)
  class.stats$area_pct <- round(class.stats$area_m2 / total_area * 100,1)
  return(class.stats)
}

#Count number of sagebrush plants.
count_sage <- function(x) {
  require(rgdal)
  setwd(x)
  setwd("scratch")  
  shape <- readOGR(dsn = ".", "sage_points")
  shape <- shape[shape$RASTERVALU > 0.015,]
  names(shape)[3] <- "Height_m"
  shape <- shape[,-2]
  setwd(x)
  writeOGR(shape, "products","sage_plant",driver = "ESRI Shapefile", overwrite_layer = TRUE)
  sage.height <- shape$Height_m
  return(sage.height)
}

#Estimate standing biomass
plot_biomass <- function(x) {
  require(raster)
  require(gdalUtils)
  require(gam)
  setwd(x)  

  #rescale aoi and sage and bareground
  align_rasters("products/aoi.tif","scratch/chm2.tif","scratch/aoi.tif", overwrite=TRUE)
  align_rasters("scratch/sage_class.tif","scratch/chm2.tif","scratch/sage_class2.tif",overwrite=TRUE)
  
  #load rescaled aoi and sage and bareground
  aoi <- raster("scratch/aoi.tif")
  sage <- raster("scratch/sage_class2.tif")
  
  #Load cHM
  chm <- raster("scratch/chm2.tif")
  
  #reclassify sage for multiplication
  sage <- reclassify(sage,
                     matrix(c(0,1,1,0),byrow = TRUE,ncol=2))
  
  #Set max height at 0.5meter , deals with greasewood or other large shrubs/trees.
  chm.herb <- sage * chm
  chm.herb[chm.herb > 0.5] <- 0.5
  
  #Figure out how many cells to aggregate across
  agg.factor <- 1/mean(res(chm.herb))
  
  #aggregate rasters to 1m
  chm.herb.1m <- aggregate(chm.herb,agg.factor)
  aoi.1m <- aggregate(aoi,agg.factor)
  
  
  #Convert rasters to vectors for analysis
  chm.herb.v <- as.vector(chm.herb.1m)
  aoi.v <- as.vector(aoi.1m)
  aoi.idx <- which(aoi.v == 1)
  
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
  
  chm.herb.v <- rep(NA,length(aoi.v))
  chm.herb.v[aoi.idx] <- chm.herb.v.test
  chm.biomass <- raster(as.matrix(chm.herb.v, nrow=nrow(chm.herb.1m), ncol=ncol(chm.herb.1m), by.row=TRUE),template=chm.herb.1m)
  writeRaster(chm.biomass,"products/biomass_est.tif",format = "GTiff",overwrite=TRUE)
  
  return(biomass.list)
}

#Estimate grass height and get statistics
get_grass_stats <- function(x) {
  require(raster)

  setwd(x)
  setwd("scratch")
  grass <- raster("grassHt_Est.tif")
  grass.v <- as.vector(grass)
  grass.v <- grass.v[!is.na(grass.v)]
  agg.factor <- 1/res(grass)[1]
  grass.1m <- aggregate(grass,fact = agg.factor,fun = mean)
  grass.1v <- as.vector(grass.1m)
  grass.1v <- grass.1v[!is.na(grass.1v)]
  
  setwd(x)
  writeRaster(grass.1m, "products/grassht_est.tif", format="GTiff", overwrite=TRUE)
  
  setwd(x)
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

copy_cameras <- function(x) {
    setwd(x)
    flist1 <- dir("scratch",pattern = "^aoiPoly.",full.names = TRUE)
    flist2 <- gsub("^scratch/","products/",flist1)
    file.copy(flist1, flist2, overwrite = TRUE)
    
    flist1 <- dir("scratch",pattern = "^Cameras.",full.names = TRUE)
    flist2 <- gsub("^scratch/","products/",flist1)
    file.copy(flist1, flist2, overwrite = TRUE)
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

gen_table_data <- function(x) {
  setwd(x)
  class.list$class.stats
  summary.table <- data.frame(
    Metric = c("Sagebrush Cover","Bareground Cover","Unclassified Cover","Sagebrush Density","Herb. Biomass","Herb. Biomass","Grass Biomass","Grass Biomass","Forb Biomass","Forb Biomass","Mean Grass Height","Mean Grass Height"),
    Units = c("Percent","Percent","Percent","Plants per acre","kg/ha","lbs/acres","kg/ha","lbs/acre","kg/ha","lbs/acre","meters","inches"),
    Result = c(class.list$class.stats$area_pct[3],class.list$class.stats$area_pct[2],class.list$class.stats$area_pct[1],class.list$sage.density,
               biomass.stats[[2]][2],biomass.stats[[3]][2],biomass.stats[[4]][2],biomass.stats[[5]][2],biomass.stats[[6]][2],biomass.stats[[7]][2],grass.stats[4],grass.stats[4]* 39.3701)
  )
  
  summary.table[1:10,3] <- round(summary.table[1:10,3],0)
  summary.table[11:12,3] <- round(summary.table[11:12,3],2)
  
  names(summary.table)[3] <- biomass.stats[[1]]
  write.csv(summary.table,"reports/plot_summary.csv",row.names = FALSE)
}
  

  
  
    
  
