# canopy_height_model.py


# Load packages
import sys
import os
import time
import pdal

# Set environment variables
DTMradius = 0.25
resolution = 0.05
CHMradius = 0.05


def processCHM(path, outname):

    # Structure the PDAL pipeline objects
    dtmJSON = f'''{{
                "pipeline": [
                    "{path}/photogrammetry/pointcloud.laz",
                    {{
                        "type":"writers.gdal",
                        "output_type":"min",
                        "dimension":"Z",
                        "resolution":"{str(resolution)}",
                        "radius":"{str(DTMradius)}",
                        "filename":"{path}/products/{outname}_dtm.tif"
                    }}
                ]
            }}'''

    chmJSON = f'''{{
                "pipeline": [
                   "{path}/photogrammetry/pointcloud.laz", 
                    {{
                        "type":"filters.hag_dem",
                        "raster":"{path}/products/{outname}_dtm.tif"
                    }},
                    {{
                        "type":"filters.ferry",
                        "dimensions":"HeightAboveGround=>Z"
                    }},
                    {{
                        "type":"filters.delaunay"
                    }},
                    {{
                        "type":"filters.faceraster",
                        "resolution":"{str(CHMradius)}"
                    }},
                    {{
                        "type":"writers.raster",
                        "filename":"{path}/products/{outname}_chm.tif"
                    }}
                ]
            }}'''

    try:
        # Calculate the DTM and export to raster
        print("Calculating DTM.....")
        dtm = pdal.Pipeline(dtmJSON)
        dtmRun = dtm.execute()

        # Calculate the canopy height model and export raster
        print("Calculating CHM.....")
        chm = pdal.Pipeline(chmJSON)
        chmRun = chm.execute()
    except RuntimeError:
        print("Error calculating DTM/CHM. Exiting...")
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

    # Check for existence of subdirectories
    if not os.path.exists(path+'/photogrammetry'):
        print("Cannot find photogrammetry subdirectory.")
        return False
    if not os.path.exists(path+'/products'):
        print("Cannot find products subdirectory.")
        return False

    # Run process
    if processCHM(path, outname):
        print("CHM successfully calculated and written to products directory.")

    # Calculate time to completion
    t1 = time.time()
    t1 -= t0
    print("Processing finished in " + "{:.2f}".format(float(t1)) + " seconds.\n")

    return

main()
