#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 11:00:38 2021
Convert .hpl files collected during the PIANO measurement campaign in 
autumn/early winter 2017 into .nc files 
    - data is neither added, nor removed, nor changed
    - .nc filenames correspond the input .hpl files
    - information about "institution" and "contact" are added to the .nc files

To run this code, the module hpl2NetCDF.py is required. It can be found in the 
GITHub repository:
    marenha/doppler_wind_lidar_toolbox/2NetCDF
    
Used directories:
    path_2NetCDF    - contains the module hpl2netCDFpy
    path_hpl        - input directory of .hpl files; sturcutre: path_hpl/[lidar_str]/Proc/yyyy/yyyymm/yyyymmdd/*.hpl
    path_l0         - output directory of .nc files; structure: path_l0/[lidar_str]/yyyy/yyyymm/yyyymmdd/*.nc
@author: maren
"""
import numpy as np
import sys,os
import matplotlib.dates as mdates
 
# Import modules
path_parent = os.path.abspath('..')
path_2NetCDF = os.path.join(path_parent,'doppler_wind_lidar_toolbox','2NetCDF')
sys.path.append(path_2NetCDF)
import hpl2NetCDF as h2n

#%% path of .hpl files and path level0 data is stored
path_hpl = os.path.join('/mnt','PIANO_campaign','lidars','raw_data')
path_l0 = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public','level0/')

#%% system id of lidars operated during PIANO
''' 
the system id lidar_str matches the directory name containing the .hpl data
'''
lidars_str =  ['SLXR_98','SLXR_142','SL_75','SL_74','SL_88']

#%% define PIANO measurement campaign period
'''
here, start_str to end_str covers the period, for which the SL_88 was installed 
at the rooftop of the HTL Anichstraße; all data available within this period 
is converted into level0 netCDF data
'''

start_str = '20170719'
end_str = '20180308'

# convert to matplotlib datenum and create day array
start_num,end_num = mdates.datestr2num(start_str), mdates.datestr2num(end_str)
days_num = np.arange(start_num,end_num+1,1)

#%% loop through lidars and days
for lidar_str in lidars_str:
    # within path_hpl_lidar, the folder structure equals the one the StreaLine software creates
    path_hpl_lidar = os.path.join(path_hpl,lidar_str,'Proc')
    
    for day_num in days_num:
        # mdates.num2datestr is not autmatically included in matplotlib.dates and has to be added manually
        date_str = mdates.num2datestr(day_num,'%Y%m%d')
        path_hpl_lidar_date = os.path.join(path_hpl_lidar,date_str[0:4],date_str[0:6],date_str)
        
        # check if data is available for this day
        if not os.path.isdir(path_hpl_lidar_date): continue
        
        # list of all data files, excluding Background_dmHH-110503.txt files
        files_all = os.listdir(path_hpl_lidar_date)
        files_all_no_bg = [file for file in files_all if '.hpl' in file and '.swp' not in file]
        
        for file_temp in files_all_no_bg:
            file_path = os.path.join(path_hpl_lidar_date,file_temp)
            
            #create .nc file and save to path_l0_lidar_date
            path_out = os.path.join(path_l0,lidar_str)
            h2n.hpl_to_netcdf(file_path,path_out,\
                        institution='University of Innsbruck, Department of Atmospheric and Cryospheric sciences (ACINN), AUSTRIA',\
                        contact='Alexander Gohm, Maren Haid (maren.haid@uibk.ac.at)',\
                        overwrite=True)
            
    



