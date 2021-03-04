#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 17:12:35 2019
Create level1 (l1) netCDF files out of level0 (l0) netCDF files
    - l1 data is corrected for not perfect alignment of the lidars to north (bearing)
    - l1 data is corrected for wrong assignment of the range gate centre (gc_offset). 
    - l1 data give infromation about lidar location (latitude, longitude, altitude)
    - l1 filenames give information about performed scan type (scan type name is user specific!!!)    
    - l0 input data was created using marenha/PIANO/netCDF_l0.py

To run this code, the module hpl2NetCDF.py and piano_paramters.py are required. 
They can be found in the GITHub repositories:
    marenha/doppler_wind_lidar_toolbox/2NetCDF
    marenha/PIANO/piano_parameters

Used directories:
    path_2NetCDF    - contains the module hpl2netCDFpy
    path_l0_in      - input directory of *_l0.nc files; structure: path_l0_in/[lidar_str]/yyyy/yyyymm/yyyymmdd/*_l0.nc
    path_l1_out     - output directory of *_l1.nc files; structure: path_l1_out/[lidar_str]_l1/yyyymmdd/*_l1.nc

Stare scan files are saved in an extra direcotry, because these hourly files will
be combined to daily files in a follow up step (marenha/PIANO/combine_stare_l1.py)

During PIANO, the SL_88 lidar performed continously 6beam scans. Therefore, scan 
files are available every 30 sec. These scans are also combined to daily files
in a follow up step (marenha/PIANO/combine_6beam_l1)
@author: maren
"""
import sys,os
import matplotlib.dates as mdates
import xarray as xr
import numpy as np

# Import modules
path_parent = os.path.abspath('..')
path_2NetCDF = os.path.join(path_parent,'doppler_wind_lidar_toolbox','2NetCDF')
sys.path.append(path_2NetCDF)
import piano_parameters as pp
import hpl2NetCDF as h2n

#%% path of input level0 files and level1 output directory
'''
The directory structure and the files names of level0 data matchs the 
directory structure and files names created by the Halo Photonics StreamLine software; 
For level1 data a modified directory structure and other files names are used
'''
path_l0_in = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public','level0/')
path_l1_out = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public','level1/')

#%% system id of lidars corresponding to the folders the raw data is stored
lidars_str=['SLXR_142','SL_75','SL_74','SL_88']
lidars_info=pp.load_lidars_info()


#%% loop through lidars
for lidar_str in lidars_str:
    lidar_info = lidars_info[lidar_str]

    path_l0_lidar=os.path.join(path_l0_in,lidar_str)
    
    # loop through all availale level0 sub-/directories 
    path_l0_lidar_yyyy = sorted(os.listdir(path_l0_lidar))

    for yyyy_temp in path_l0_lidar_yyyy:
        path_l0_lidar_yyyymm = sorted(os.listdir(os.path.join(path_l0_lidar,yyyy_temp)))

        for yyyymm_temp in path_l0_lidar_yyyymm:
            path_l0_lidar_yyyymmdd = sorted(os.listdir(os.path.join(path_l0_lidar,yyyy_temp,yyyymm_temp)))
    
            for yyyymmdd_temp in path_l0_lidar_yyyymmdd:
                date_num = mdates.datestr2num(yyyymmdd_temp)
                
                # data is only converted to level1 data if metadata is valid for this time period
                if (date_num<lidar_info.start_num) | (date_num>lidar_info.end_num): continue
                
                
                path_l0_lidar_date = os.path.join(path_l0_lidar,yyyy_temp,yyyymm_temp,yyyymmdd_temp)
                files_all = os.listdir(path_l0_lidar_date)
                
                
                #%% Convert Stare files to level1
                '''
                In this step, stare data is saved in an extra folder, this data
                will be combined into daily files (combine_stare_l1.py) in the 
                next step and stored in the same directory as the other level1 
                scan files
                '''
                folder_name = '%s_stare_l1_hourly' % lidar_info.name
                files_stare = [file for file in files_all if 'Stare' in file]
                path_out = os.path.join(path_l1_out,folder_name,yyyymmdd_temp)
                
                for file_temp in files_stare:
                    hour_temp = file_temp.split('_')[3]
                    file_name_out = '%s_stare_l1_%s_%s.nc' %(lidar_info.name,yyyymmdd_temp,hour_temp)
                    file_path = os.path.join(path_l0_lidar_date,file_temp)
                    h2n.to_netcdf_l1(file_path,file_name_out,lidar_info,path_out)
                
                #%% Convert User files to level1 data
                '''
                This part depends very much on the performed scan scenarios and
                scan pattern; the naming and methods used here cannot directly
                transfered to other campaign and use cases;
                SL88 performed continously 6beam scans during PIANO;
                As for the stare data, the 6beam scans are first converted into
                l1 single netCDF files and later combined for a whole day with 
                combine_6beam_l1.py
                '''
                files_user = [file for file in files_all if 'User' in file]
                if lidar_str in ['SL_75','SLXR_142','SL_74']:
                    '''
                    Information about scan scenarios during PIANO is imported; 
                    the scan types included in these scenarios are part of the 
                    level1 file names; scans which are not part of any scenario 
                    are not converted into level1 data 
                    '''
                    scan_types_s0split=pp.load_scaninfo('scenario0_split')
                    scan_types_s0=pp.load_scaninfo('scenario0')
                    scan_types_s1=pp.load_scaninfo('scenario1')
                    scan_types_sc=pp.load_scaninfo('scenario_christmans')
                    
                    #loop through user files
                    for file_temp in files_user:
                        file_path = os.path.join(path_l0_lidar_date,file_temp)
                        # import scan type name for each user file
                        with xr.open_dataset(file_path) as ds:
                            scan_file_temp=ds.scan_type.split()[0]
                    
                        # check if scan type is part of any scan scenario
                        scan_type_temp_s0split_list = [scan_type for scan_type in scan_types_s0split if scan_file_temp in scan_type.scan_files]
                        scan_type_temp_s0_list = [scan_type for scan_type in scan_types_s0 if scan_file_temp in scan_type.scan_files]
                        scan_type_temp_s1_list = [scan_type for scan_type in scan_types_s1 if scan_file_temp in scan_type.scan_files]
                        scan_type_temp_sc_list = [scan_type for scan_type in scan_types_sc if scan_file_temp in scan_type.scan_files]
                        
                        scan_types_all = [scan_type_temp_s0split_list,scan_type_temp_s0_list,scan_type_temp_s1_list,scan_type_temp_sc_list]
                        scan_types_all_inside = [len(scan_type) == 1 for scan_type in scan_types_all]
                        
                        # only if file can be assigned to a scan scenario, level1 data is created
                        if any(scan_types_all): 
                            
                            #if a scan_type is included in several scenarios, the first one is taken
                            # TODO: find more elegant solution for this selection
                            scan_types_all_inside_int=np.where(scan_types_all_inside)[0][0]
                            scan_type_temp=scan_types_all[scan_types_all_inside_int][0]
                            
                            # output directory
                            folder_name=folder_name = '%s_l1' % (lidar_info.name)
                            path_out = os.path.join(path_l1_out,folder_name,yyyymmdd_temp)
                            
                            # file name of level1 netCDF data file
                            file_name_out = '%s_%s_l1_%s_%s.nc' %(lidar_info.name,scan_type_temp.name,file_temp.split('_')[2],file_temp.split('_')[3])
                            
                            h2n.to_netcdf_l1(file_path,file_name_out,lidar_info,path_out)
                            
                        else:
                            print('%s: %s (%s) cannot assigned to scan type' %(file_temp,yyyymmdd_temp,scan_file_temp))
                            continue
                
                # create l1 netCDF files for SL88 lidar
                elif lidar_str == 'SL_88':
                    for file_temp in files_user:
                        
                        file_path = os.path.join(path_l0_lidar_date,file_temp)
                        folder_name=folder_name = '%s_l1' % (lidar_info.name)
                        path_out = os.path.join(path_l1_out,folder_name,yyyymmdd_temp)
                        
                        file_name_out = '%s_6beam_l1_%s_%s.nc' %(lidar_info.name,file_temp.split('_')[2],file_temp.split('_')[3])

                        h2n.to_netcdf_l1(file_path,file_name_out,lidar_info,path_out)

                    
                