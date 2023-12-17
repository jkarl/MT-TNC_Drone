# metashape_processing.py

# Metashape License Key
licenseKey = "A8RBJ-Y33UL-A6A9H-45H3P-ZZZN4"

# Load libraries
import sys
import os
import time
import Metashape
import pdal
from osgeo import gdal
from shapely.geometry import MultiPoint, mapping
from shapely.ops import transform
import fiona
#import pyproj

# Check and set metashape licensing
if not Metashape.app.activated:
    try:
        Metashape.License().activate(licenseKey)
    except:
        print("Cannot connect to the license server, or no activated license for Metashape.")


# set environment variables
path = 'C:/Users/jkarl/Downloads/2022/PL1' # switch to an argument
TYPES = ["jpg", "jpeg", "tif", "tiff"]
proj_crs = Metashape.CoordinateSystem("EPSG::32100") # output all data in State Plane (m)
projection = Metashape.OrthoProjection()
projection.crs = proj_crs
accuracy = 1  #align photos accuracy with High Accuracy (downscale=1)
reference_preselection = True #use photo coordinates to help with pair selection
generic_preselection = False #Reference preselection on, Generic Off
use_adaptive_fitting = False
keypoints = 40000 #align photos key point limit
tiepoints = 4000 #align photos tie point limit
threshold_RU = 10  # threshold value for Reconstruction uncertainity
threshold1_PA = 3  # threshold value for projection accuracy , pass 1
threshold2_PA = 2  # threshold value for projection accuracy , pass 2
accuracy_TP = float(0.3)  # update tiepoint_accuracy
threshold1_RE = 0.5  # threshold value for reprojection accuracy, pass 1
threshold2_RE = 0.4  # threshold value for reprojection accuracy, pass 2
threshold3_RE = 0.35  # threshold value for reprojection accuracy, pass 3
threshold4_RE = 0.3  # threshold value for reprojection accuracy, pass 4
source_dsm = Metashape.DataSource.DenseCloudData  # build mesh/DEM source
quality = 4  # build dense cloud at medium quality (downscale=4; high quality would be downscale of 2)
filtering = Metashape.MildFiltering  # depth filtering
interpolation = Metashape.Interpolation.EnabledInterpolation  # build mesh interpolation
colorBalance = True  # color balance for ortho creation
blendingMode = Metashape.BlendingMode.MosaicBlending  # how to blend orthomosaic
colorCorrection = True  # use color correction in building orthomosaic
source_mesh = Metashape.DataSource.PointCloudData
surface = Metashape.SurfaceType.Arbitrary #build mesh surface type
face_num = Metashape.FaceCount.LowFaceCount #build mesh polygon count
atlas_size = 4096

def process(path, outname):
    # stuff here...

    print("Processing " + path)

    # Get list of photos to add
    imagery = path + "/rawImagery"
    list_files = os.listdir(imagery)
    list_photos = list()
    for entry in list_files:  # finding image files
        file = imagery + "/" + entry
        if os.path.isfile(file):
            if file[-3:].lower() in TYPES:
                list_photos.append(file)

    # set up project
    doc = Metashape.Document()
    doc.save(path+"/photogrammetry/projectPS.psx")
    chunk = doc.addChunk()
    chunk.label = path.rsplit("/", 1)[1]
    chunk.crs = proj_crs


    ###align photos
    print("Aligning photos...")
    chunk.addPhotos(list_photos)
    chunk.matchPhotos(downscale=accuracy, generic_preselection=generic_preselection,
                      reference_preselection=reference_preselection, filter_mask=False, keypoint_limit=keypoints,
                      tiepoint_limit=tiepoints)
    chunk.alignCameras(adaptive_fitting=use_adaptive_fitting)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)  # Adaptive fitting may not work. Not in API reference
    chunk.resetRegion()
    doc.save()

    ## optimize sparse cloud
    ###Gradual Selection
    print("Optimizing sparse cloud...")
    f = Metashape.PointCloud.Filter()
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ReconstructionUncertainty)
    f.removePoints(threshold_RU)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)

    ##run second time for cleanup
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ReconstructionUncertainty)
    f.removePoints(threshold_RU)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)
    doc.save()

    # Next Filter based on projection accuracy.
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ProjectionAccuracy)
    f.removePoints(threshold1_PA)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)

    # run a second time for cleanup
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ProjectionAccuracy)
    f.removePoints(threshold2_PA)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)
    doc.save()

    # Tighten tiepoint accuracy.
    chunk.tiepoint_accuracy = accuracy_TP
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)
    doc.save()

    # Filter Points based on reprojection error. Use 4 steps to avoid selection too many points
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ReprojectionError)
    f.removePoints(threshold1_RE)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ReprojectionError)
    f.removePoints(threshold2_RE)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ReprojectionError)
    f.removePoints(threshold3_RE)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)
    f.init(chunk, criterion=Metashape.PointCloud.Filter.ReprojectionError)
    f.removePoints(threshold4_RE)
    chunk.optimizeCameras(adaptive_fitting=use_adaptive_fitting)
    doc.save()

    # Save camera reference data for later. Used for generating AOI.
    print("Saving camera reference data")
    chunk.exportReference(path + "/photogrammetry/reference.csv", Metashape.ReferenceFormatCSV,
                        items=Metashape.ReferenceItemsCameras, delimiter=",")

    camerasX = []
    camerasY = []
    for camera in chunk.cameras:
        camerasX.append(camera.reference.location[0])
        camerasY.append(camera.reference.location[1])
    cameraLocs = list(zip(camerasX,camerasY))

    # Create polygon object
    points = MultiPoint(cameraLocs)
    convexHull = points.convex_hull

    # Define a polygon feature geometry with one attribute
    schema = {
        'geometry': 'Polygon',
        'properties': {'id': 'int'},
    }

    # Write a new Shapefile
    with fiona.open(path+'/photogrammetry/aoi.shp', 'w', 'ESRI Shapefile', schema) as c:
        c.write({
            'geometry': mapping(convexHull),
            'properties': {'id': 1},
        })

    # Add the shape to the Metashape Project
    os.chdir(path+'/photogrammetry')
    chunk.importShapes('aoi.shp', boundary_type=Metashape.Shape.OuterBoundary,
                       crs=Metashape.CoordinateSystem("EPSG::4326"))

    # Build Dense Cloud
    print("Building dense cloud...")
    chunk.buildDepthMaps(downscale=quality, filter_mode=filtering)
    chunk.buildDenseCloud(point_colors=True, keep_depth=False)
    doc.save()

    # export pointcloud
    print("Exporting dense cloud...")
    chunk.exportPoints(path=path + "/photogrammetry/pointcloud.laz",
                       source_data=Metashape.DataSource.DenseCloudData,
                       save_colors=True,
                       crs=proj_crs,
                       clip_to_boundary=True)


    ## build dem/dsm
    print("Building DEM...")
    chunk.buildDem(source_data=Metashape.DataSource.DenseCloudData,
                   interpolation=interpolation,
                   projection=projection)
    doc.save()

    ## export dsm
    print("Exporting DEM...")
    chunk.exportRaster(source_data=Metashape.DataSource.ElevationData,
                    path=path + "/products/"+outname+"_dsm.tif",
                    projection=projection,
                    clip_to_boundary=True)

    ###building mesh for ortho2
    chunk.buildModel(surface_type = surface, source_data = source_mesh, interpolation = interpolation, face_count = face_num)
    doc.save()

    ## build ortho1
    print("Building Orthomosaic...")
    chunk.calibrateColors(source_data=Metashape.DataSource.ModelData)
    chunk.buildOrthomosaic(surface_data=Metashape.ElevationData,
                           blending_mode=blendingMode,
                           fill_holes=True,
                           projection=projection)
    doc.save()

    ## export ortho1
    print("Exporting orthomosaic...")
    compression = Metashape.ImageCompression()
    compression.tiff_compression = Metashape.ImageCompression.TiffCompressionNone
    chunk.exportRaster(source_data=Metashape.DataSource.OrthomosaicData,
                            path=path + "/products/"+outname+"_ortho.tif",
                            projection=projection,
                            save_alpha=True,
                            clip_to_boundary=True,
                            image_compression=compression)

    print("Writing processing report...")
    chunk.exportReport(path=path + "/reports/Metashape.pdf")



    print("Finished Processing " + chunk.label)
    return True

def setupDirectories(path):
    # Check for and create if necessary the directory structure for the outputs
    try:
        if not os.path.exists(path+"/photogrammetry"): os.mkdir(path+"/photogrammetry")
        if not os.path.exists(path + "/logs"): os.mkdir(path + "/logs")
        if not os.path.exists(path + "/products"): os.mkdir(path + "/products")
    except RuntimeError:
        print("Error: could not create directory structure")
        return False
    return True

def main():
    # Set timer
    t0 = time.time()
    print("Script started...")

    # Check for arguments
    if len(sys.argv) < 3:
        print("No valid path input. Script aborted.")
        return False
    if os.path.isdir(sys.argv[1]):
        path = sys.argv[1]
        outname = sys.argv[2]
    else:
        print("No valid path input. Script aborted.")
        return False

    # Check for Metashape active license. Quit if not found.
    try:
        if not Metashape.app.activated:
            raise Exception
    except:
        print("Cannot connect to the license server, or no activated license for Metashape.")
        return False

    # Set directory structure
    if not setupDirectories(path):
        print("Could not set up directory structure. Quitting.")
        return False

    # Run process
    process(path, outname)

    # Calculate time to completion
    t1 = time.time()
    t1 -= t0
    print("Image processing finished in " + "{:.2f}".format(float(t1)) + " seconds.\n")

    return

main()