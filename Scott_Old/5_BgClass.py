# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("Spatial")

arcpy.env.workspace = sys.argv[1]


# Local variables:
bg_pred_tif = "scratch/bgpred.tif"
Input_true_raster_or_constant_value = "1"
Input_false_raster_or_constant_value = "0"
Con_tif1 = "scratch/con1.tif"
Con_tif2 = "scratch/con2.tif"
Majorit_Con_1 = "scratch/maj1.tif"
Boundar_Majo2 = "scratch/bound1.tif"
RegionG_Boun2 = "scratch/region1.tif"
Input_true_raster_or_constant_value__2_ = "0"
Con_RegionG_1 = "scratch/bg_class.tif"
aoi = "products/aoi.tif"

# Process: Con
arcpy.gp.Con_sa(bg_pred_tif, Input_true_raster_or_constant_value, Con_tif1, Input_false_raster_or_constant_value, "VALUE > 0.5")

# Process: Times
arcpy.gp.Times_sa(Con_tif1, aoi, Con_tif2)


# Process: Majority Filter
arcpy.gp.MajorityFilter_sa(Con_tif2, Majorit_Con_1, "EIGHT", "MAJORITY")

# Process: Boundary Clean
arcpy.gp.BoundaryClean_sa(Majorit_Con_1, Boundar_Majo2, "ASCEND", "TWO_WAY")

# Process: Region Group
arcpy.gp.RegionGroup_sa(Boundar_Majo2, RegionG_Boun2, "FOUR", "WITHIN", "ADD_LINK", "")

# Process: Con (2)
arcpy.gp.Con_sa(RegionG_Boun2, Input_true_raster_or_constant_value__2_, Con_RegionG_1, Boundar_Majo2, "LINK = 1 AND COUNT < 8")

# Process: Build Pyramids
arcpy.BuildPyramids_management(Con_RegionG_1, "-1", "NONE", "NEAREST", "DEFAULT", "75", "OVERWRITE")

# Process: Calculate Statistics
arcpy.CalculateStatistics_management(Con_RegionG_1, "1", "1", "", "OVERWRITE")

arcpy.Delete_management(Con_tif1)
arcpy.Delete_management(Con_tif2)
arcpy.Delete_management(Majorit_Con_1)
arcpy.Delete_management(Boundar_Majo2)
arcpy.Delete_management(RegionG_Boun2)
arcpy.Delete_management(bg_pred_tif)


