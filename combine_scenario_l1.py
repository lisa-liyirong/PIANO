#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 11:58:45 2019
Prepare level1 coordinated scan files in order to calc 2D retrievals
here, only scenario0_split is considered
!!! this code can only be used for scenario0_split; for the other scenarios
manuel chagnes are necessary!!!
@author: maren
"""
import numpy as np
import xarray as xr
import os,sys
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

import halo_data_operations as hdo
import piano_parameters as pp
#import piano_utilities as pu
z_ref=570 
dtn=24*60*60

#%% path of level1 files and output folder
path_l1 = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public','level1/')


lidars_info = pp.load_lidars_info()
scan_scenario = 'scenario0_split'
scan_types = pp.load_scaninfo(scan_scenario)

path_out = os.path.join(path_l1,scan_scenario)

'''
in the netCDF file created in this code, the lidars location and the measurement points
are defined in dependent of on certain point. Here, the SLXR142 lidar is chosen 
since it is part of each coordinated scan pattern
'''
origin=lidars_info['SLXR_142']

for scan_type in scan_types:
    if scan_type.name == 'vad24': continue

    lidars_temp=scan_type.lidars
    
    for lidar_temp in lidars_temp:
        lidar_info = lidars_info[lidar_temp]
        path_l1_lidar_scan = os.path.join(path_l1,'%s_l1' %lidar_info.name)
        
        for yyyymmdd_temp in sorted(os.listdir(path_l1_lidar_scan)):
            files_path = os.path.join(path_l1_lidar_scan,yyyymmdd_temp)
            
            files_all = os.listdir(files_path)
            files_scan = sorted([file for file in files_all if scan_type.name in file])
            
            fn=len(files_scan)
            fi=0
            while fi < fn:
                print(fi)
                file_temp=files_scan[fi]
                file_path = os.path.join(files_path,file_temp)
                ds_list = []
                with xr.open_dataset(file_path,decode_times=False) as ds:
                    scan_file_temp = ds.scan_type.split()[0]
                    start_file_temp = mdates.datestr2num(ds.start_time)
                    ds_list.append(ds)
                    
                # only do this operation for scan_types included in the chosen scan scenario
                if scan_file_temp in scan_type.scan_files:
                    #%% Combine splitted scans to one pre prepaired data file
                    if (fi+1) < fn:
                        file_path_next = os.path.join(files_path,files_scan[fi+1])
                        with xr.open_dataset(file_path_next,decode_times=False) as ds:
                            scan_file_temp_next = ds.scan_type.split()[0]
                            start_file_temp_next = mdates.datestr2num(ds.start_time)

                        #check if next scan_file is part of the actual scan pattern
                        if (np.abs(start_file_temp_next-start_file_temp)*dtn < (18*60)) and (scan_file_temp_next in scan_type.scan_files):
                            ds_list.append(ds)
                            
                            fi+=1
                        

                    data=xr.concat(ds_list,dim='NUMBER_OF_RAYS')
            
                    data['center_of_gate']=xr.DataArray(data.gate_centers.values[0,:],dims=['NUMBER_OF_GATES'], \
                          attrs={'units':'m','long_name':'center of range gates 1d'})
            
                    data.attrs['dlx']=lidar_info.x-origin.x
                    data.attrs['dly']=lidar_info.y-origin.y
                    data.attrs['dlz']=lidar_info.zsl
                    
                    #shift angle to middle of averaging interval 
                    #angles of ppi with passes 0 degree is shifted
                    data = hdo.correct_angles(lidar_temp,data,correct_bearing=False)
                    
                    data = hdo.separate_scans_l1(data)    

                    path_netcdf_out_temp=os.path.join(path_out,yyyymmdd_temp)
                    if not os.path.isdir(path_netcdf_out_temp):  os.makedirs(path_netcdf_out_temp)
                    
                    '''
                    the file_name of the output data contains the date and hour, the lidar scan is assigned to. E.g, a scan
                    performed between 11:20 to 11:40 UTC belongs to the 12th hour of the day --> hour=12
                    this makes it easier to find coodinated scans from two lidars to estimate retrievals
                    '''
                    hour_num = np.ceil(start_file_temp*24)/24
                    file_name = '%s_%s_l1_%s.nc' % (lidar_info.name,scan_type.name,mdates.num2datestr(hour_num,'%Y%m%d_%H'))
                    
                    print('%s saved' %file_name)
                    
                    data.to_netcdf(os.path.join(path_netcdf_out_temp,file_name))
                    data.close()
                fi+=1

