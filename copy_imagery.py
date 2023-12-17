# copy_imagery.py
# copies /rawimagery directories from the existing TNC file structure to a new directory

import os
import shutil

basepath = 'F:/UAS/MT-TNC/uavPlots'
outpath = 'F:/UAS/MT-TNC/rerunPlots'

# create output directory if it doesn't exist
if not os.path.exists(outpath):
    os.mkdir(outpath)

# get list of directories to copy raw imagery from
listcontents = os.listdir(basepath)
listdir = []
for i in listcontents:
    if os.path.isdir(os.path.join(basepath,i)):
        listdir.append(i)

# set up a list to hold directory names that do not have rawimagery subdirectories
faillist = []

# iterate over all the directories
for dir in listdir:

    # check for existence of rawimagery subdirectory
    if os.path.exists(os.path.join(basepath,dir,"rawimagery")):

        #print update
        print("Copying raw imagery from "+basepath+" to "+outpath)

        # create new directory in outpath
        if not os.path.exists(os.path.join(outpath,dir)):
            os.mkdir(os.path.join(outpath,dir))

        # copy rawimagery folder over to new directory in outpath
        shutil.copytree(os.path.join(basepath,dir,"rawimagery"),
                     os.path.join(outpath,dir,"rawimagery"))

    else:
        print("No raw imagery folder in "+basepath)
        faillist.append(dir)

print("Finished!")