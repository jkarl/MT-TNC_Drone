# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy

arcpy.env.workspace = sys.argv[1]

# Local variables:
prelim_class_tif = "scratch/prelim_class.tif"
rat_class_dbf = "scratch/rat_class.dbf"

# Process: Build Raster Attribute Table
arcpy.BuildRasterAttributeTable_management(prelim_class_tif, "NONE")

# Process: Join Field
arcpy.JoinField_management(prelim_class_tif, "Value", rat_class_dbf, "Value", "Class")

