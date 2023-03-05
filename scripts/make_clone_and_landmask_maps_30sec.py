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



# original ldd and clone
global_ldd_inp_file = "/scratch/depfg/sutan101/data/global_ldd_reservoirs_and_lakes_before_data_lost_on_eejit/lddsound_30sec_version_202005XX.map"

# global subdomain file
global_subdomain_file = "/scratch/depfg/sutan101/make_global_subdomains/version_2021-02-17/general_subdomains_using_threshold_of_25_cells/global_subdomains_30min_final_filled.map"

# output_folder
out_folder            = "/scratch/depfg/sutan101/make_global_subdomains/version_2021-02-17/general_subdomains_using_threshold_of_25_cells/30sec/"

# cell size in arcmin
# - 0.5 arcmin = 30 arcsec
cellsize_in_arcmin = 0.5


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
    print("set the clone based on the ldd inpu") 
    pcr.setclone(global_ldd_inp_file)
    

    # define the landmask
    print("define the landmask based on the ldd inpup") 
    landmask = pcr.defined(pcr.readmap(global_ldd_inp_file))
    landmask = pcr.ifthen(landmask, landmask)
    

    # read ldd
    print("define the ldd") 
    ldd_map = pcr.readmap(global_ldd_inp_file)
    # ~ # - extend ldd (not needed)
    # ~ ldd_map = pcr.ifthen(landmask, pcr.cover(ldd_map, pcr.ldd(5)))
    
    # copy ldd file
    cmd = "cp " + str(global_ldd_inp_file) + " ."
    print(cmd); os.system(cmd)


    # make catchment map
    print("make catchment map") 
    catchment_map = pcr.catchment(ldd_map, pcr.pit(ldd_map))


    # read global subdomain file
    print("read global subdomain file") 
    global_subdomain_map = vos.readPCRmapClone(v = global_subdomain_file, cloneMapFileName = global_ldd_inp_file, tmpDir = tmp_folder, absolutePath = None, isLddMap = False, cover = None, isNomMap = True)


    # set initial subdomain
    print("assign subdomains to all catchments") 
    subdomains_initial = pcr.areamajority(global_subdomain_map, catchment_map)
    subdomains_initial = pcr.ifthen(landmask, subdomains_initial)

    pcr.aguila(subdomains_initial)

    pcr.report(subdomains_initial, "global_subdomains_initial.map")

    print(str(int(vos.getMinMaxMean(pcr.scalar(subdomains_initial))[0])))
    print(str(int(vos.getMinMaxMean(pcr.scalar(subdomains_initial))[1])))


    
    print("Checking all subdomains, avoid too large subdomains") 

    num_of_masks = int(vos.getMinMaxMean(pcr.scalar(subdomains_initial))[1])

    # clone code that will be assigned
    assigned_number = 0
    
    subdomains_final = pcr.ifthen(pcr.scalar(subdomains_initial) < -7777, pcr.nominal(0))
    
    for nr in range(1, num_of_masks + 1, 1):

        msg = "Processing the landmask %s" %(str(nr))
        print(msg)

        mask_selected_boolean = pcr.ifthen(subdomains_initial == nr, pcr.boolean(1.0))
        
        process_this_clone = False
        if pcr.cellvalue(pcr.mapmaximum(pcr.scalar(mask_selected_boolean)), 1, 1)[0] > 0: process_this_clone = True
        
        # ~ if nr == 1: pcr.aguila(mask_selected_boolean)
        
        # - initial check value
        check_ok = True

        if process_this_clone:
            xmin, ymin, xmax, ymax = boundingBox(mask_selected_boolean)
            area_in_degree2 = (xmax - xmin) * (ymax - ymin)
            
            # ~ print(str(area_in_degree2))
            
            # check whether the size of bounding box is ok
            reference_area_in_degree2 = 2500.
            if area_in_degree2 > 1.50 * reference_area_in_degree2: check_ok = False
            if (xmax - xmin) > 10* (ymax - ymin): check_ok = False
        
        if check_ok == True and process_this_clone == True:

            msg = "Clump is not needed."
            msg = "\n\n" +str(msg) + "\n\n"
            print(msg)

            # assign the clone code
            assigned_number = assigned_number + 1

            # update global landmask for river and land
            mask_selected_nominal = pcr.ifthen(mask_selected_boolean, pcr.nominal(assigned_number))
            subdomains_final = pcr.cover(subdomains_final, mask_selected_nominal) 
        
        if check_ok == False and process_this_clone == True:
			
            msg = "Clump is needed."
            msg = "\n\n" +str(msg) + "\n\n"
            print(msg)

            # make clump
            clump_ids = pcr.nominal(pcr.clump(mask_selected_boolean))
            
            # merge clumps that are close together 
            clump_ids_window_majority = pcr.windowmajority(clump_ids, 10.0)
            clump_ids = pcr.areamajority(clump_ids_window_majority, clump_ids) 
            # ~ pcr.aguila(clump_ids)
            
            # minimimum and maximum values
            min_clump_id = int(pcr.cellvalue(pcr.mapminimum(pcr.scalar(clump_ids)),1)[0])
            max_clump_id = int(pcr.cellvalue(pcr.mapmaximum(pcr.scalar(clump_ids)),1)[0])

            for clump_id in range(min_clump_id, max_clump_id + 1, 1):
            
                msg = "Processing the clump %s of %s from the landmask %s" %(str(clump_id), str(max_clump_id), str(nr))
                msg = "\n\n" +str(msg) + "\n\n"
                print(msg)

                # identify mask based on the clump
                mask_selected_boolean_from_clump = pcr.ifthen(clump_ids == pcr.nominal(clump_id), mask_selected_boolean)
                mask_selected_boolean_from_clump = pcr.ifthen(mask_selected_boolean_from_clump, mask_selected_boolean_from_clump)

                # check whether the clump is empty
                check_mask_selected_boolean_from_clump = pcr.ifthen(mask_selected_boolean, mask_selected_boolean_from_clump)
                check_if_empty = float(pcr.cellvalue(pcr.mapmaximum(pcr.scalar(pcr.defined(check_mask_selected_boolean_from_clump))),1)[0])
                
                if check_if_empty == 0.0: 
                
                    msg = "Map is empty !"
                    msg = "\n\n" +str(msg) + "\n\n"
                    print(msg)

                else:
                
                    msg = "Map is NOT empty !"
                    msg = "\n\n" +str(msg) + "\n\n"
                    print(msg)

                    # assign the clone code
                    assigned_number = assigned_number + 1
                    
                    # update global landmask for river and land
                    mask_selected_nominal = pcr.ifthen(mask_selected_boolean_from_clump, pcr.nominal(assigned_number))
                    subdomains_final = pcr.cover(subdomains_final, mask_selected_nominal)
                    
    # ~ # kill all aguila processes if exist
    # ~ os.system('killall aguila')
    
    pcr.aguila(subdomains_final)

    print("")
    print("")
    print("")

    print("The subdomain map is READY.") 

    pcr.report(subdomains_final, "global_subdomains_final.map")

    num_of_masks = int(vos.getMinMaxMean(pcr.scalar(subdomains_final))[1])
    print(num_of_masks)

    print("")
    print("")
    print("")


    print("Making the clone and landmask maps for all subdomains") 

    num_of_masks = int(vos.getMinMaxMean(pcr.scalar(subdomains_final))[1])

    # clone and mask folders
    clone_folder = out_folder + "/clone/"
    if os.path.exists(clone_folder): shutil.rmtree(clone_folder)
    os.makedirs(clone_folder)
    mask_folder = out_folder + "/mask/"
    if os.path.exists(mask_folder): shutil.rmtree(mask_folder)
    os.makedirs(mask_folder)


    print("")
    print("")

    for nr in range(1, num_of_masks + 1, 1):

        msg = "Processing the subdomain %s" %(str(nr))
        print(msg)

        # set the global clone
        pcr.setclone(global_ldd_inp_file)
        
        mask_selected_boolean = pcr.ifthen(subdomains_final == nr, pcr.boolean(1.0))
        
        mask_selected_nominal = pcr.ifthen(subdomains_final == nr, pcr.nominal(nr))
        mask_file = "mask/mask_%s.map" %(str(nr))
        pcr.report(mask_selected_nominal, mask_file)

        xmin, ymin, xmax, ymax = boundingBox(mask_selected_boolean)
        area_in_degree2 = (xmax - xmin) * (ymax - ymin)
        
        print(str(nr) + " ; " + str(area_in_degree2) + " ; " + str((xmax - xmin)) + " ; " + str((ymax - ymin)))

        # cellsize in arcdegree 
        cellsize = cellsize_in_arcmin / 60.
        
        # number of rows and cols
        num_rows = int(round(ymax - ymin) / cellsize)
        num_cols = int(round(xmax - xmin) / cellsize)
        
        # make the clone map using mapattr 
        clonemap_mask_file = "clone/clonemap_mask_%s.map" %(str(nr))
        cmd = "mapattr -s -R %s -C %s -B -P yb2t -x %s -y %s -l %s %s" %(str(num_rows), str(num_cols), str(xmin), str(ymax), str(cellsize), clonemap_mask_file)
        print(cmd); os.system(cmd)
        
        # set the local landmask for the clump
        pcr.setclone(clonemap_mask_file)
        local_mask = vos.readPCRmapClone(v = mask_file, \
                                         cloneMapFileName = clonemap_mask_file, 
                                         tmpDir = tmp_folder, \
                                         absolutePath = None, isLddMap = False, cover = None, isNomMap = True)
        local_mask_boolean = pcr.defined(local_mask)
        local_mask_boolean = pcr.ifthen(local_mask_boolean, local_mask_boolean)
        pcr.report(local_mask_boolean, mask_file)
        


    print("")
    print("")
    print("")

    print(num_of_masks)

        
if __name__ == '__main__':
    sys.exit(main())
