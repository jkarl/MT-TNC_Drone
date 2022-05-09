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
SageCHM = "Times_tif4.tif"
SageCHM2 = "SetNull_Time2.tif"
Neg_CHM = "Negate_SetNu2.tif"
Input_raster_or_constant_value_2 = "10"
inverse_chm = "Plus_Negate_2.tif"
FlowDir = "FlowDir_2.tif"
FlowAcc = "FlowAcc.tif"
Focal1 = "FocalSt_Flow3.tif"
Focal2 = "SetNull_x2.tif"
Local_Min = "Minus_SetNul1.tif"
Input_true_raster_or_constant_value = "1"
Sage_Apex1 = "Con_Minus_Se1.tif"
Sage_Apex2 = "RasterT_Con_Min2.shp"
Sage_Apex3 = "sage_points.shp"

# Process: Times
arcpy.gp.Times_sa(chm2_tif, sage_class_tif, SageCHM)

# Process: Set Null
arcpy.gp.SetNull_sa(SageCHM, SageCHM, SageCHM2, "VALUE = 0")

# Process: Negate
arcpy.gp.Negate_sa(SageCHM2, Neg_CHM)

# Process: Plus
arcpy.gp.Plus_sa(Neg_CHM, Input_raster_or_constant_value_2, inverse_chm)

# Process: Flow Direction
arcpy.gp.FlowDirection_sa(inverse_chm, FlowDir, "NORMAL","","D8")

# Process: Flow Accumulation
arcpy.gp.FlowAccumulation_sa(FlowDir, FlowAcc, "", "FLOAT", "D8")

# Process: Focal Statistics
arcpy.gp.FocalStatistics_sa(FlowAcc, Focal1, "Circle 0.5 MAP", "MAXIMUM", "DATA")

# Process: Set Null (2)
arcpy.gp.SetNull_sa(SageCHM, Focal1, Focal2, "VALUE = 0")

# Process: Minus
arcpy.gp.Minus_sa(Focal2, FlowAcc, Local_Min)

# Process: Con
arcpy.gp.Con_sa(Local_Min, Input_true_raster_or_constant_value, Sage_Apex1, "", "VALUE = 0")

# Process: Raster to Point
arcpy.RasterToPoint_conversion(Sage_Apex1, Sage_Apex2, "VALUE")

# Process: Extract Values to Points
arcpy.gp.ExtractValuesToPoints_sa(Sage_Apex2, chm2_tif, Sage_Apex3, "NONE", "VALUE_ONLY")

arcpy.Delete_management(SageCHM)
arcpy.Delete_management(SageCHM2)
arcpy.Delete_management(Neg_CHM)
arcpy.Delete_management(inverse_chm)
arcpy.Delete_management(FlowDir)
arcpy.Delete_management(FlowAcc)
arcpy.Delete_management(Focal1)
arcpy.Delete_management(Focal2)
arcpy.Delete_management(Local_Min)
arcpy.Delete_management(Sage_Apex1)
arcpy.Delete_management(Sage_Apex2)
