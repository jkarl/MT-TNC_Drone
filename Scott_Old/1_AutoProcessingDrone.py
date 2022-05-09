#compatibility Agisoft Metashape Professional 1.4.x
#use argument to specify the path to the "master folder"

import Metashape, os, sys, time, math

try:
    if not Metashape.app.activated:
        raise Exception
except:
    print("Cannot connect to the license server, or no activated license for Metashape.")

def process(path):

	# Enable GPU processing if possible...
	Metashape.app.gpu_mask = 2 ** len(Metashape.app.enumGPUDevices()) - 1 #setting GPU mask
	if Metashape.app.gpu_mask:
		Metashape.app.cpu_enable = False  
	else:
		Metashape.app.cpu_enable = True
		
	
	### processing parameters
	proj_crs = Metashape.CoordinateSystem("EPSG::32100") #output all data in State Plane (m)
	accuracy = 1  #align photos accuracy with High Accuracy (downscale=1)
	reference_preselection = True #use photo coordinates to help with pair selection
	generic_preselection = False #Reference preselection on, Generic Off
	use_adaptive_fitting = False
	keypoints = 40000 #align photos key point limit
	tiepoints = 4000 #align photos tie point limit
	threshold_RU = 10 #threshold value for Reconstruction uncertainity
	threshold1_PA = 3  #threshold value for projection accuracy , pass 1
	threshold2_PA = 2  #threshold value for projection accuracy , pass 2
	accuracy_TP = float(0.3) #update tiepoint_accuracy
	threshold1_RE = 0.5 #threshold value for reprojection accuracy, pass 1
	threshold2_RE = 0.4 #threshold value for reprojection accuracy, pass 2
	threshold3_RE = 0.35 #threshold value for reprojection accuracy, pass 3
	threshold4_RE = 0.3 #threshold value for reprojection accuracy, pass 4
	source_dsm = Metashape.DataSource.DenseCloudData #build mesh/DEM source
	#source_mesh = Metashape.DataSource.PointCloudData
	#surface = Metashape.SurfaceType.Arbitrary #build mesh surface type
	quality = 4 #build dense cloud at medium quality (downscale=4; high quality would be downscale of 2)
	filtering = Metashape.FilterMode.MildFiltering #depth filtering
	interpolation = Metashape.Interpolation.EnabledInterpolation #build mesh interpolation
	#blending = Metashape.BlendingMode.MosaicBlending #blending mode
	#face_num = Metashape.FaceCount.LowFaceCount #build mesh polygon count
	#mapping = Metashape.MappingMode.GenericMapping #build texture mapping
	colorBalance = True #color balance for ortho creation
	blendingMode = Metashape.BlendingMode.MosaicBlending #how to blend orthomosaic
	colorCorrection = True #use color correction in building orthomosaic
	atlas_size = 4096
	TYPES = ["jpg", "jpeg", "tif", "tiff"]
	
	###end of processing parameters definition

	print("Processing " + path)
	imagery = path + "/rawImagery"
	list_files = os.listdir(imagery)
	list_photos = list()
	for entry in list_files: #finding image files
		file = imagery + "/" + entry
		if os.path.isfile(file):
			if file[-3:].lower() in TYPES:
				list_photos.append(file)
	
	doc = Metashape.Document()	
	doc.save(path + "photogrammetry/projectPS.psx")
	chunk = doc.addChunk()
	chunk.label = path.rsplit("/", 1)[1]
	chunk.crs = proj_crs
	
	###align photos
	chunk.addPhotos(list_photos)
	chunk.matchPhotos(downscale=accuracy, generic_preselection=generic_preselection, reference_preselection=reference_preselection, filter_mask=False, keypoint_limit=keypoints, tiepoint_limit=tiepoints)
	chunk.alignCameras(adaptive_fitting = use_adaptive_fitting)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting) # Adaptive fitting may not work. Not in API reference
	chunk.resetRegion()
	doc.save()	
				
	###Gradual Selection
	f = Metashape.PointCloud.Filter()
	f.init(chunk, criterion = Metashape.PointCloud.Filter.ReconstructionUncertainty)
	f.removePoints(threshold_RU)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)

	##run second time for cleanup
	f.init(chunk, criterion = Metashape.PointCloud.Filter.ReconstructionUncertainty)
	f.removePoints(threshold_RU)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)
	doc.save()

	#Next Filter based on projection accuracy.
	f.init(chunk, criterion = Metashape.PointCloud.Filter.ProjectionAccuracy)
	f.removePoints(threshold1_PA)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)
	
	#run a second time for cleanup
	f.init(chunk, criterion = Metashape.PointCloud.Filter.ProjectionAccuracy)
	f.removePoints(threshold2_PA)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)
	doc.save()

	#Tighten tiepoint accuracy.
	chunk.tiepoint_accuracy = accuracy_TP
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)
	doc.save()

	#Filter Points based on reprojection error. Use 4 steps to avoid selection too many points
	f.init(chunk, criterion = Metashape.PointCloud.Filter.ReprojectionError)
	f.removePoints(threshold1_RE)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)
	f.init(chunk, criterion = Metashape.PointCloud.Filter.ReprojectionError)
	f.removePoints(threshold2_RE)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)
	f.init(chunk, criterion = Metashape.PointCloud.Filter.ReprojectionError)
	f.removePoints(threshold3_RE)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)
	f.init(chunk, criterion = Metashape.PointCloud.Filter.ReprojectionError)
	f.removePoints(threshold4_RE)
	chunk.optimizeCameras(adaptive_fitting = use_adaptive_fitting)
	doc.save()

	###building dense cloud
	chunk.buildDepthMaps(downscale = quality, filter = filtering)
	chunk.buildDenseCloud(point_colors = True, keep_depth = False)
	doc.save()

	#export pointcloud
	chunk.exportPoints(path = path + "/photogrammetry/pointcloud.laz",
		source = Metashape.DataSource.DenseCloudData,
		colors = True,
		projection = proj_crs)

	###building mesh for ortho2
	#chunk.buildModel(surface = surface, source = source_mesh, interpolation = interpolation, face_count = face_num)
	#doc.save()

	## build dem/dsm
	chunk.buildDem(source = source_dsm, 
		interpolation = interpolation,
		projection = proj_crs)
	doc.save()

	## export dsm
	chunk.exportDem(path = path + "/photogrammetry/dsm.tif",
		projection = proj_crs,
		tiff_big = True,
		tiff_overviews = True)	

	## build ortho1
	chunk.calibrateColors(source_data=Metashape.DataSource.ModelData, color_balance=colorBalance)
	chunk.buildOrthomosaic(surface=Metashape.ElevationData,
		blending=blendingMode,
		fill_holes = True, 
		projection = proj_crs)
	doc.save()
	
	## export ortho1
	chunk.exportOrthomosaic(path = path  + "/photogrammetry/ortho1.tif",
		projection = proj_crs,
		tiff_big = True,
		tiff_overviews = True,
		write_alpha = True)

	## build ortho2 
	#chunk.calibrateColors(source_data=Metashape.DataSource.ModelData, color_balance=colorBalance)
	#chunk.buildOrthomosaic(surface=Metashape.ModelData,
	#	blending=blendingMode,
	#	fill_holes = True,
	#	projection = proj_crs)
	#doc.save()

	## export ortho2
	#chunk.exportOrthomosaic(path = path  + "/photogrammetry/ortho2.tif",
	#	projection = proj_crs,
	#	tiff_big = True,
	#	tiff_overviews = True,
	#	write_alpha = True)

	chunk.exportReport(path = path + "/logs/Metashape.pdf")

	## Save camera reference data for later. Used for generating AOI.
	chunk.saveReference(path + "/photogrammetry/reference.csv", Metashape.ReferenceFormatCSV, items = Metashape.ReferenceItemsCameras, delimiter = ",")

	print("Processed " + chunk.label)
	return True

def main():

	t0 = time.time()
	print("Script started...")

	if len(sys.argv) < 2:
		print("No valid path input. Script aborted.")
		return False
	if os.path.isdir(sys.argv[1]):
		path = sys.argv[1]
	else:
		print("No valid path input. Script aborted.")
		return False	
	
	#Run process 
	process(path)

	t1 = time.time()
	t1 -= t0
	t1 = float(t1)	

	print("Image processing finished in " + "{:.2f}".format(t1) + " seconds.\n")

	return

main()
		
