# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy

arcpy.env.overwriteOutput = True

arcpy.env.workspace = sys.argv[1]

arcpy.CheckOutExtension("Spatial")

# Local variables:
products = sys.argv[1]

chm = "chm2.tif"
hill = "chm2_hillshade.tif"
dem = "dem.tif"
dem_hill = "dem_hillshade.tif"

# Process: Hillshade
arcpy.gp.HillShade_sa(chm, hill, "315", "45", "NO_SHADOWS", "1")

# Process: Hillshade
arcpy.gp.HillShade_sa(dem, dem_hill, "315", "45", "NO_SHADOWS", "1")

# Process: Build Pyramids And Statistics
arcpy.BuildPyramidsandStatistics_management(products, "INCLUDE_SUBDIRECTORIES", "BUILD_PYRAMIDS", "CALCULATE_STATISTICS", "NONE", "", "NONE", "1", "1", "", "-1", "NONE", "NEAREST", "DEFAULT", "75", "SKIP_EXISTING", "")

