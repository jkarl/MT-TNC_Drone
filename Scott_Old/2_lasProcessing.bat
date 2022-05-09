::
:: 2018-07-23
:: Processes raw pointclound data from photoscan DPC
:: This script does the following
:: 1. Sorts the points for faster reading and visualization and rescales to cm from mm
:: 2. Removes 'low noise' points to help with generating Bare-earth DEM
:: 3. Outputs a cleaned pointclound dataset 
:: 4. Runs a groundclassification to remove high points (non-ground)
:: 5. Creates a thinned pointclound dataset from 4 that is used to generate bareground  DEM
:: 6. Outputs a thinned pointcloud and a shapefile of thinned points
::
:: First Argument is the input file and second arguement is the number of cores(threads) to use.

::create temporary tiles directory

rmdir temp_tiles /s /q
mkdir temp_tiles

:: create a temporary tiling of the file into 7 by 7 meter tiles 
:: use a 10 meter buffer. Specify projection as EPSG 32100

lastile -i %1 ^
        -tile_size 7 ^
        -epsg 32100 ^
        -buffer 5 ^
        -flag_as_withheld ^
        -rescale 0.01 0.01 0.01 ^
        -odir temp_tiles -o split.laz

echo tiling complete
:: run spatial coherence for faster processing

rmdir temp_tiles_sorted /s /q
mkdir temp_tiles_sorted

lassort -v -i temp_tiles\split*.laz ^
        -odir temp_tiles_sorted -olaz ^
        -cores %2

echo lassort complete

rmdir temp_tiles /s /q

:: reclassify the highest point of every 2.5 by 2.5 meter grid cell with classification code 8. rmdir temp_tiles_thinned /s /q

rmdir temp_tiles_thinned /s /q
mkdir temp_tiles_thinned

lasthin -v -i temp_tiles_sorted\*.laz -step 0.5 ^
        -highest -classify_as 8 ^
        -odir temp_tiles_thinned -olaz ^
        -cores %2

echo lasthin complete 
rmdir temp_tiles_sorted /s /q

::Considering only those points classified as 8 in the last step we then run 
::lasnoise to find points that are highly isolated in wide and flat neighborhoods 
::that are then reclassified as 7.

rmdir temp_tiles_isolated /s /q
mkdir temp_tiles_isolated

lasnoise -v -i temp_tiles_thinned\split*.laz ^
        -ignore_class 0 ^
        -step_xy 1 -step_z 0.2 -isolated 4 ^
        -classify_as 7 ^
        -odir temp_tiles_isolated -olaz ^
        -cores %2
echo lasnoise complete
rmdir temp_tiles_thinned /s /q

::Now we run a temporary ground classification of only (!!!) on those points that are still 
::classified as 8 using the default parameters of lasground. Hence we only use the points 
::that were the highest points on the 2.5 by 2.5 meter grid and that were not classified 
::as noise in the previous step.

rmdir temp_tiles_ground /s /q
mkdir temp_tiles_ground

lasground -v -i temp_tiles_isolated\split*.laz ^
        -city -ultra_fine -ignore_class 0 7 ^
        -odir temp_tiles_ground -olaz ^
        -cores %2

echo lasground complete
rmdir temp_tiles_isolated /s /q

::The result of this temporary ground filtering is then merely used to mark all points that are 0.5 
::meter below the triangulated TIN of these temporary ground points with classification code 12 using lasheight.

rmdir temp_tiles_denoised /s /q         
mkdir temp_tiles_denoised

lasheight -v -i temp_tiles_ground\split*.laz ^
        -do_not_store_in_user_data ^
        -classify_below -0.5 12 ^
        -odir temp_tiles_denoised -olaz ^
        -cores %2

echo lasheight complete
rmdir temp_tiles_ground /s /q

:: Reclassify points back to zero and remove noise (class 12) 

rmdir temp_tiles_denoised2 /s /q         
mkdir temp_tiles_denoised2

las2las -v -i temp_tiles_denoised\split*.laz ^
        -change_classification_from_to 1 0 ^
        -change_classification_from_to 2 0 ^
        -change_classification_from_to 7 0 ^
        -drop_class 12 ^
        -odir temp_tiles_denoised2 -olaz ^
        -cores %2

echo last2las complete
rmdir temp_tiles_denoised /s /q      

:: Remove buffer from tiles and output a clean pointcloud (all points)

rmdir temp_out1 /s /q
mkdir temp_out1

lastile -v -i temp_tiles_denoised2\split*.laz ^
        -remove_buffer ^
        -odir temp_out1 -olaz ^
        -cores %2

echo lastile complete
lasmerge -v -i temp_out1\split*.laz ^
        -o pointcloud_clean.laz ^
        -olaz

echo lasmerge complete
rmdir temp_out1 /s /q

:: Run lasground to identify non-ground points in clean pointcloud data

rmdir temp_tiles_ground2 /s /q
mkdir temp_tiles_ground2

lasground -v -i temp_tiles_denoised2\split*.laz ^
        -odir temp_tiles_ground2 -olaz ^
        -archaeology ^
        -cores %2

echo lasground two complete
rmdir temp_tiles_denoised2 /s /q

::drop class 1 (not ground)

rmdir temp_tiles_drop /s /q
mkdir temp_tiles_drop

las2las -v -i temp_tiles_ground2\split*.laz ^
        -drop_class 1 ^
        -odir temp_tiles_drop -olaz ^
        -cores %2

echo las2las complete
rmdir temp_tiles_ground2 /s /q

:: Run lasthin to select lowest points in 2x2 meter step

rmdir temp_tiles_thin2 /s /q   
mkdir temp_tiles_thin2

lasthin -v -i temp_tiles_drop\split*.laz ^
        -step 1 ^
        -odir temp_tiles_thin2 ^
        -olaz ^
        -cores %2

echo lasthin two complete
rmdir temp_tiles_drop /s /q

:: Remove buffer on thinned points

rmdir temp_tiles_out2 /s /q
mkdir temp_tiles_out2

lastile -v -i temp_tiles_thin2\split*.laz ^
        -remove_buffer ^
        -odir temp_tiles_out2 -olaz ^
        -cores %2

echo lastile remove buffers complete
rmdir temp_tiles_thin2 /s /q

:: merge and write out thinned pointcloud
lasmerge -v -i temp_tiles_out2\split*.laz ^
        -o pointcloud_ground.laz ^
        -olaz
echo final merge complete
rmdir temp_tiles_out2 /s /q

:: create shapefile of thinned points
las2shp -v -i pointcloud_ground.laz -o points_ground.shp -single_points
echo shapefile created