#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 10:45:29 2018
Calculation of two-dimensional wind fields from level1 netCDF data and export
as .nc files. This code only works for the PIANO scan scenario scenario1b (details
see piano_parameters.py). The inout data files were prepared with combine_scenario_l1.py.

To run this code, the module hpl2NetCDF.py is required. It can be found in the 
GITHub repository:
    marenha/doppler_wind_lidar_toolbox/coplanar_retrieval
@author: maren
"""
import numpy as np
import xarray as xr
import os,sys
import matplotlib.dates as mdates
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd

# Import modules
path_parent = os.path.abspath('..')
path_cr = os.path.join(path_parent,'doppler_wind_lidar_toolbox','coplanar_retrieval')
sys.path.append(path_cr)
import hpl2NetCDF as h2n
import halo_data_operations as hdo
import calc_retrieval as cr
import piano_parameters as pp


#%% input directory (path_l1) and output directory (path_l2)
path_l1 = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public','level1/')
path_l2 = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public','level2/')
# path_plot_out=os.path.join('/mnt','PIANO_campaign','lidars','quicklooks','vads')

#%% import information about the instruments and the scan strategy
lidars_info=pp.load_lidars_info()


scan_scenario = 'scenario1b'
scan_types_s0split=pp.load_scaninfo(scan_scenario)


weight='lidar'
delta_l=30
ground_level_ibk=570 #approx. height of Inn River near UIBK 

path_l1_scenario = os.path.join(path_l1,scan_scenario)
for scan_type in scan_types_s0split:
    path_l2_scan = os.path.join(path_l2,'%s_l2' %scan_type.name)
    #here, no ret are estimated for conical scans
    if scan_type.orientation in ['conical','rhisn','rhiew']:
        continue

    path_l1_yyyymmdd = sorted(os.listdir(path_l1_scenario))
        
    for yyyymmdd_temp in path_l1_yyyymmdd:
        date_num = mdates.datestr2num(yyyymmdd_temp)
        
        
        path_l2_scan_date = os.path.join(path_l2_scan,yyyymmdd_temp)
        if not os.path.exists(path_l2_scan_date):
            os.makedirs(path_l2_scan_date)
        
        path_l1_scenario_temp = os.path.join(path_l1_scenario,yyyymmdd_temp)
        files_all = os.listdir(path_l1_scenario_temp)

        files_scan=[file for file in files_all if scan_type.name in file]

        #group files by scan time; does only work for this type of naming: lidarName_scanType_l1_yyyymmdd_HH.nc
        files_times_start=['%s_%s' %(file.split('_')[3],file.split('_')[4].split('.')[0]) for file in files_scan]
        scan_times=np.unique(files_times_start)
        
        
        files_scan_sort=[sorted([file for file in files_scan if scan_time in file]) for scan_time in scan_times]
        sn=len(files_scan_sort)
        for scan_time,files_scan_temp in zip(scan_times,files_scan_sort):   
            
            
            if len(files_scan_temp) != len(scan_type.lidars): continue
          
            # number of scans has to be the same
            sn_all=np.unique([xr.open_dataset(os.path.join(path_l1_scenario_temp,file),decode_times=False).NUMBER_OF_SCANS.size for file in files_scan_temp])
            if sn_all.size>1: continue
        
            calc_ret=True
            si=0
            retrievals_episode=[]
            scan_start_episode,scan_end_episode=[],[]
            retrieval_start_episode,retrieval_end_episode=[],[]
            scan_episode=[]
            while calc_ret:  
                #write file in 
                scan_list,dl_loc_list,dl_id_list=[],[],[]
                scan_start_list,scan_end_list=[],[]
                
                
                for lidar_temp in scan_type.lidars:
                    # probably not necessary step to sort by lidar
                    lidar_info = lidars_info[lidar_temp]
                    file_temp = [file for file in files_scan_temp if lidar_info.name in file][0]
                    with xr.open_dataset(os.path.join(path_l1_scenario_temp,file_temp),decode_times=False) as data_temp:
                        sn=data_temp.NUMBER_OF_SCANS.size
                        if si>=sn: continue
                        el_temp=data_temp.el_scan.sel(NUMBER_OF_SCANS=si).values
                        az_temp=data_temp.az_scan.sel(NUMBER_OF_SCANS=si).values
                        vr_temp=data_temp.radial_velocity_scan.sel(NUMBER_OF_SCANS=si).values
                        intensity_temp=data_temp.intensity_scan.sel(NUMBER_OF_SCANS=si).values
                        snr_temp=intensity_temp-1
                        delta_g=int(data_temp.range_gate_length.split()[0])
                        r=data_temp.center_of_gate.values
                        scan_start_list.append(data_temp.start_scan.sel(NUMBER_OF_SCANS=si).values/24+date_num)
                        scan_end_list.append(data_temp.end_scan.sel(NUMBER_OF_SCANS=si).values/24+date_num)
                    
                    vr_filtered_temp=hdo.noise_reduce(lidar_info.name,vr_temp,intensity_temp,False)
    #                vr_filtered_temp_th=hdo.noise_reduce(lidar_id_temp,vr_temp,intensity_temp,True)
                    
                    dl_loc=[data_temp.dlx,data_temp.dly,data_temp.dlz]
                    dl_loc_list.append(dl_loc)
                    dl_id_list.append(lidar_info.name)
                    scan_list.append(cr.scan(el_temp,az_temp,vr_filtered_temp,snr_temp,dl_loc,r))

                    
                    
                grid=pp.load_grid(scan_type.name,delta_l)      
                [scan.to_grid(grid) for scan in scan_list]
                scan_start,scan_end=np.min(scan_start_list),np.max(scan_end_list)
                
                retrieval_temp=cr.calc_retrieval(scan_list,grid,weight)
    
                scan_episode.append(scan_list)
                retrievals_episode.append(retrieval_temp)
                retrieval_start_episode.append(scan_start)
                retrieval_end_episode.append(scan_end)
                scan_start_episode.append(scan_start_list)
                scan_end_episode.append(scan_end_list)

                si+=1
                if si==sn:
                    calc_ret=False 
            
    
    
            #dimensions
            ln=len(scan_list)
            sn=len(scan_episode)
            (gyn,gxn)=grid.xx.shape
            
            vr_grid_array=np.full([ln,sn,gyn,gxn],np.nan)
            for si,scan_list in enumerate(scan_episode):
                for li,scan in enumerate(scan_list):
                    vr_grid_array[li,si,:,:]=scan.vr_grid
                
            #combine retrievals to multi dimensional arrays
            ws_array=np.array([retrieval.ws for retrieval in retrievals_episode])
            u_array=np.array([retrieval.u for retrieval in retrievals_episode])
            v_array=np.array([retrieval.v for retrieval in retrievals_episode])
            
            scan_start_array=np.array(retrieval_start_episode)
            scan_end_array=np.array(retrieval_end_episode)
            
            dlx=np.array([scan[0] for scan in dl_loc_list])
            dly=np.array([scan[1] for scan in dl_loc_list])
            dlz=np.array([scan[2] for scan in dl_loc_list])
            
    
            file_name='%s_l2_%s.nc' %(scan_type.name,scan_time)
            file_path_out=os.path.join(path_l2_scan_date,file_name)
            if os.path.isfile(file_path_out):
                os.remove(file_path_out)
            dataset=Dataset(file_path_out,'w',format='NETCDF4_CLASSIC')
            
            NUMBER_OF_LIDARS=dataset.createDimension('NUMER_OF_LIDARS',ln)
            NUMBER_OF_SCANS= dataset.createDimension('NUMBER_OF_SCANS',sn)
            X_DIM=dataset.createDimension('X_DIM',gxn)
            Y_DIM=dataset.createDimension('Y_DIM',gyn)
           
            x_var=dataset.createVariable('x',np.float64,('X_DIM')) 
            y_var=dataset.createVariable('y',np.float64,('Y_DIM')) 
            
            dlx_var=dataset.createVariable('dlx',np.float64,('NUMER_OF_LIDARS')) 
            dly_var=dataset.createVariable('dly',np.float64,('NUMER_OF_LIDARS')) 
            dlz_var=dataset.createVariable('dlz',np.float64,('NUMER_OF_LIDARS')) 
            
            ws_var=dataset.createVariable('ws',np.float64,\
                                          ('NUMBER_OF_SCANS','Y_DIM','X_DIM')) 
            u_var=dataset.createVariable('u',np.float64,\
                                          ('NUMBER_OF_SCANS','Y_DIM','X_DIM')) 
            v_var=dataset.createVariable('v',np.float64,\
                                          ('NUMBER_OF_SCANS','Y_DIM','X_DIM')) 
            vr_grid_var=dataset.createVariable('vr_grid',np.float64,\
                                          ('NUMER_OF_LIDARS','NUMBER_OF_SCANS','Y_DIM','X_DIM')) 
            
            scan_start_var=dataset.createVariable('scan_start',np.float64,('NUMBER_OF_SCANS')) 
            scan_end_var=dataset.createVariable('scan_end',np.float64,('NUMBER_OF_SCANS')) 
            
            vr_grid_var[:]=vr_grid_array
            dlx_var[:]=dlx
            dly_var[:]=dly
            dlz_var[:]=dlz
            if scan_type.name in ['rhisn','rhiew']:
                x_var[:]=grid.xy
                y_var[:]=grid.z
            else:
                x_var[:]=grid.x
                y_var[:]=grid.y
            ws_var[:]=ws_array
            u_var[:]=u_array
            v_var[:]=v_array
            scan_start_var[:]=scan_start_array
            scan_end_var[:]=scan_end_array
            
            dataset.description='Retrieved two dimensional wind field from coplanar Doppler wind lidar scans'
            dataset.weight=weight
            dataset.lattice_length=delta_l
            
            dataset.close()
        
    
    
            
            
    
