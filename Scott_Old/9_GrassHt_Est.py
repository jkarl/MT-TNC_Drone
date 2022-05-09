# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy

arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("Spatial")

arcpy.env.workspace = sys.argv[1]

# Local variables:
chm2_tif = "chm2.tif"
sage_class_tif = "sage_class.tif"
Con1_True = "0"
con1_tif = "con1.tif"
Expand_tif = "expand.tif"
Con2_True = "0"
Con2_False = "1"
Con2_tif = "con2.tif"
Times1_tif = "times1.tif"
Focal_Tif = "focal.tif"
Times2_Const = "1.5001"
GrassHt_Fine = "grassHt_Est.tif"

# Process: Con1
arcpy.gp.Con_sa(chm2_tif, Con1_True, con1_tif, chm2_tif, "VALUE > 0.5")

# Process: Expand
arcpy.gp.Expand_sa(sage_class_tif, Expand_tif, "4", "1")

# Process: con2
arcpy.gp.Con_sa(Expand_tif, Con2_True, Con2_tif, Con2_False, "VALUE = 1")

# Process: Times1
arcpy.gp.Times_sa(con1_tif, Con2_tif, Times1_tif)

# Process: FocalStat
arcpy.gp.FocalStatistics_sa(Times1_tif, Focal_Tif, "Circle 0.1524 MAP", "MAXIMUM", "DATA")

# Process: Times2
arcpy.gp.Times_sa(Focal_Tif, Times2_Const, GrassHt_Fine)


arcpy.Delete_management(con1_tif)
arcpy.Delete_management(Expand_tif)
arcpy.Delete_management(Con2_tif)
arcpy.Delete_management(Times1_tif)
arcpy.Delete_management(Focal_Tif)

