# metashape_test.py
# Script for testing metashape API function in python and checking for valid license

try:
    import Metashape
except:
    print("Cannot load Metashape library")

try:
    if not Metashape.app.activated:
        raise Exception
    else:
        print("Metashape library loaded and license validated. Ready to go!")
except:
    print("Cannot connect to the license server, or no activated license for Metashape.")


