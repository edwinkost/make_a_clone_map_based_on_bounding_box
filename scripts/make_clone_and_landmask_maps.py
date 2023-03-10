#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from __future__ import print_function

import os
import sys
import glob
import subprocess
import time
import datetime
import shutil

import numpy as np
import math

import pcraster as pcr
import virtualOS as vos


def boundingBox(pcrmap):
    ''' derive the bounding box for a map, return xmin,ymin,xmax,ymax '''
    bb = []
    xcoor = pcr.xcoordinate(pcrmap)
    ycoor = pcr.ycoordinate(pcrmap)
    xmin  = pcr.cellvalue(pcr.mapminimum(xcoor), 1, 1)[0]
    xmax  = pcr.cellvalue(pcr.mapmaximum(xcoor), 1, 1)[0]
    ymin  = pcr.cellvalue(pcr.mapminimum(ycoor), 1, 1)[0]
    ymax  = pcr.cellvalue(pcr.mapmaximum(ycoor), 1, 1)[0]
    return [math.floor(xmin), math.floor(ymin), math.ceil(xmax), math.ceil(ymax)]
    
def spatialInterpolation2PCR(fieldArray, pcrType, MV):
	#-interpolates the field array to the full extent
	field = pcr.numpy2pcr(pcrType, fieldArray, MV)
	cellID = pcr.nominal(pcr.uniqueid(pcr.defined(field)))
	zoneID = pcr.spreadzone(cellID,0,1)
	if pcrType == pcr.Scalar:
		field = pcr.areaaverage(field,zoneID)
	else:
		field = pcr.areamajority(field,zoneID)
	return field

def define_landmask(input_file, clone_map_file, output_map_file):

    # define the landmask based on the input     
    cmd = "gdalwarp -tr 0.5 0.5 -te -180 -90 180 90 -r max " + str(input_file) + " " + output_map_file + ".tif"
    print(cmd); os.system(cmd)
    cmd = "gdal_translate -of PCRaster " + output_map_file + ".tif " + output_map_file
    print(cmd); os.system(cmd)
    cmd = "mapattr -c " + clone_map_file + " " + output_map_file
    print(cmd); os.system(cmd)
    cmd = "rm " + output_map_file + ".*"
    print(cmd); os.system(cmd)
    
    landmask = pcr.defined(pcr.readmap(output_map_file))
    landmask = pcr.ifthen(landmask, landmask)
    # ~ pcr.aguila(landmask)
    
    return landmask


# bounding box extent (-te xmin ymin xmax ymax)
# - thailand: latitudes: 4-22 and longitudes: 95-107
xmin =  97. #  95. #  95. #  95. #  97.
ymin =   5. #   4. #   4. #   4. #   5.
xmax = 106. # 107. # 110. # 107. # 106.
ymax =  21. #  22. #  22. #  22. #  21.

# ldd file
global_ldd_inp_file = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/routing/surface_water_bodies/version_2020-05-XX/lddsound_30sec_version_202005XX.map"

# output folder
out_folder = "/scratch/depfg/sutan101/thailand_30sec/clone_maps/version_20220310/"

# output files
out_clone_file = "clone_thailand.map"
out_mask_file  =  "mask_thailand.map"

# output files in a full path
out_clone_file = out_folder + "/" + out_clone_file
out_mask_file  = out_folder + "/" + out_mask_file

def main():

    # output folder
    clean_out_folder = True
    if os.path.exists(out_folder): 
        if clean_out_folder:
            shutil.rmtree(out_folder)
            os.makedirs(out_folder)
    else:
        os.makedirs(out_folder)
    os.chdir(out_folder)    
    os.system("pwd")

    # tmp folder
    tmp_folder = out_folder + "/tmp/"
    if os.path.exists(tmp_folder): shutil.rmtree(tmp_folder)
    os.makedirs(tmp_folder)
    
    # set the clone map
    print("set the clone based on the ldd input") 
    pcr.setclone(global_ldd_inp_file)

    # read ldd
    print("read the ldd") 
    ldd_map = pcr.readmap(global_ldd_inp_file)

    # ~ # - extend ldd (not needed)
    # ~ ldd_map = pcr.ifthen(landmask, pcr.cover(ldd_map, pcr.ldd(5)))

    # get the cell size/resolution based on the ldd map
    cellsize = pcr.clone().cellSize()
    
    # create the bounding box
    print("create the bounding box") 
    bounding_box = pcr.boolean(1.0)
    bounding_box = pcr.ifthen(pcr.xcoordinate(bounding_box) >= xmin, bounding_box)
    bounding_box = pcr.ifthen(pcr.ycoordinate(bounding_box) >= ymin, bounding_box)
    bounding_box = pcr.ifthen(pcr.xcoordinate(bounding_box) <= xmax, bounding_box)
    bounding_box = pcr.ifthen(pcr.ycoordinate(bounding_box) <= ymax, bounding_box)
    
    # include the catchment map of the bounding box
    print("include upstream areas of the bounding box") 
    bounding_box_catchment = pcr.catchment(ldd_map, bounding_box)
    bounding_box_catchment = pcr.ifthen(bounding_box_catchment, bounding_box_catchment)
    
    # option to use only cells that have either 'complete' catchments or 'complete/correct' upstream areas (TODO: make this optional)
    print("using only cells that have either 'complete' catchments or 'complete/correct' upstream areas") 
    #
    # - define ldd at bounding box only
    ldd_map_at_bounding_box = pcr.lddrepair(pcr.lddmask(ldd_map, bounding_box))
    #
    # - catchment/upstream area sizes (at bounding box only)
    ldd_map_at_bounding_box_catchment_size = pcr.catchmenttotal(pcr.scalar(1.0), ldd_map_at_bounding_box)
    # -- if complete catchments are preferred
    # --- catchment classes (at bounding box only)
    ldd_map_at_bounding_box_catchment = pcr.catchment(ldd_map_at_bounding_box, pcr.pit(ldd_map_at_bounding_box))
    ldd_map_at_bounding_box_catchment = pcr.ifthen(pcr.scalar(ldd_map_at_bounding_box_catchment) > 0.0, ldd_map_at_bounding_box_catchment)
    # --- complete catchment size (at bounding box only)
    ldd_map_at_bounding_box_catchment_size = pcr.areamaximum(ldd_map_at_bounding_box_catchment_size, ldd_map_at_bounding_box_catchment)
    #
    # - correct catchment sizes (based on the global ldd) 
    ldd_map_catchment_size = pcr.catchmenttotal(pcr.scalar(1.0), ldd_map)
    # -- if complete catchments are preferred
    ldd_map_catchment_size = pcr.areamaximum(ldd_map_catchment_size, ldd_map_at_bounding_box_catchment)
    # ~ #
    # - bounding box used 
    bounding_box_catchment = ldd_map_at_bounding_box_catchment_size == ldd_map_catchment_size
    bounding_box_catchment = pcr.ifthen(bounding_box_catchment, bounding_box_catchment)

    # option to use only the largest clump
    print("using only the largest clump") 
    clump_ids   = pcr.clump(bounding_box_catchment)
    clump_areas = pcr.areatotal(pcr.scalar(bounding_box_catchment), clump_ids)
    bounding_box_catchment = pcr.ifthen(clump_areas == pcr.mapmaximum(clump_areas), pcr.boolean(1.0))
    bounding_box_catchment = pcr.ifthen(bounding_box_catchment, bounding_box_catchment)
    
    # report the bounding_box_catchment as a scalar map
    pcr.report(pcr.scalar(bounding_box_catchment), out_mask_file)

    # checking using aguila
    cmd = "aguila " +  out_mask_file + " + " + global_ldd_inp_file
    print(cmd)
    os.system(cmd) 
    
    # bounding box coordinates for the clone
    print("get the bounding box coordinates for the clone map") 
    xmin_clone, ymin_clone, xmax_clone, ymax_clone = boundingBox(bounding_box_catchment)
    # - number of rows and columns
    num_rows = int( round( (ymax_clone - ymin_clone) / cellsize ))
    num_cols = int( round( (xmax_clone - xmin_clone) / cellsize ))
    
    # make the clone map using mapattr 
    clonemap_mask_file = out_clone_file
    cmd = "mapattr -s -R %s -C %s -B -P yb2t -x %s -y %s -l %s %s" %(str(num_rows), str(num_cols), str(xmin_clone), str(ymax_clone), str(cellsize), clonemap_mask_file)
    print(cmd); os.system(cmd)

    # set the local landmask for the clump
    pcr.setclone(clonemap_mask_file)
    local_mask = vos.readPCRmapClone(v = out_mask_file, \
                                         cloneMapFileName = clonemap_mask_file, 
                                         tmpDir = tmp_folder, \
                                         absolutePath = None, isLddMap = False, cover = None, isNomMap = False)
    local_mask_boolean = pcr.defined(local_mask)
    local_mask_boolean = pcr.ifthen(local_mask_boolean, local_mask_boolean)
    pcr.report(local_mask_boolean, out_mask_file)

        
if __name__ == '__main__':
    sys.exit(main())
