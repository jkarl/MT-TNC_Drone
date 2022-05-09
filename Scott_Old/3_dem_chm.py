# -*- coding: utf-8 -*-

#SCRIPT TAKES ONE ARGUMENT AT THE COMMAND LINE: WORKING Directory

# Import arcpy module
import arcpy  

# Set working directory
arcpy.env.workspace = sys.argv[1]

arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")

# Local variables: Inputs
ortho1 = ".\\photogrammetry\\ortho1.tif"
ortho2 = ".\\photogrammetry\\ortho2.tif"
DSM = ".\\photogrammetry\\dsm.tif"
groundpoints = ".\\scratch\\points_ground.shp"
aoiPoly = ".\\scratch\\aoiPoly.shp"

#Local variables: outputs
DEM = ".\\scratch\\dem.tif"
tCHM = ".\\scratch\\tchm.tif"
tCHM2 = ".\\scratch\\tchm2.tif"

DEMout = ".\\products\\dem.tif"
DSMout = ".\\products\\dsm.tif"
CHMout = ".\\products\\chm.tif"
aoi = ".\\products\\aoi.tif"
o1out = ".\\products\\ortho1.tif"
o2out = ".\\products\\ortho2.tif"

#variables other
Input_false_raster_or_constant_value = "0"

#Set Env Variables for initial processing
arcpy.env.snapRaster = DSM
arcpy.env.extent = DSM
arcpy.env.cellSize = DSM
arcpy.env.overwriteOutput = True
field = ['value']

#First generate DEM from ground points
# Process: Natural Neighbor
print "Generating prelim DEM"
arcpy.NaturalNeighbor_3d(groundpoints, "Shape", DEM, "")

#Process: Minus 
#Generate preliminary CHM by subtracting DEM from DSM
print "Generating prelim CHM"
arcpy.Minus_3d(DSM, DEM, tCHM)

# Process: Remove negative values from CHM.
print "Remove negative CHM values"
#arcpy.gp.Con_sa(tCHM, tCHM, CHM, Input_false_raster_or_constant_value, "VALUE >= 0")
arcpy.gp.Con_sa(tCHM, tCHM, tCHM2, Input_false_raster_or_constant_value, "VALUE > 0")

#Delete intermediate data
print "Delete intermediate data"
arcpy.Delete_management(tCHM)

#covert AOI Poly to raster
arcpy.env.extent = aoiPoly
arcpy.AddField_management(aoiPoly,"value","SHORT")
with arcpy.da.UpdateCursor(aoiPoly,field) as cursor:
	for row in cursor:
		if row[0] == 0:
			row[0] = 1
			cursor.updateRow(row)

arcpy.PolygonToRaster_conversion(aoiPoly, "value", aoi)

#crop files for output
desc      = arcpy.Describe(aoi)
ExtObj    = desc.extent
clip      = "%d %d %d %d" % (ExtObj.XMin, ExtObj.YMin, ExtObj.XMax, ExtObj.YMax)

arcpy.Clip_management(DEM,clip,DEMout,aoiPoly,"#","ClippingGeometry","MAINTAIN_EXTENT")
arcpy.Clip_management(tCHM2,clip,CHMout,aoiPoly,"#","ClippingGeometry","MAINTAIN_EXTENT")
arcpy.Clip_management(DSM,clip,DSMout,aoiPoly,"#","ClippingGeometry","MAINTAIN_EXTENT")
arcpy.Clip_management(ortho1,clip,o1out,aoiPoly,"#","ClippingGeometry","MAINTAIN_EXTENT")
arcpy.Clip_management(ortho2,clip,o2out,aoiPoly,"#","ClippingGeometry","MAINTAIN_EXTENT")

#build statistics 

products = ".\\products"
arcpy.BuildPyramidsandStatistics_management(products, "NONE", "BUILD_PYRAMIDS", "CALCULATE_STATISTICS", "NONE", "", "NONE", "1", "1", "", "-1", "NONE", "NEAREST", "DEFAULT", "75", "SKIP_EXISTING", "")

#cleanup

arcpy.Delete_management(tCHM)
arcpy.Delete_management(tCHM2)
arcpy.Delete_management(DEM)

#EOF


