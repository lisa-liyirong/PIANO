#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 16:44:42 2018
Estimate vertical profiles of horizontal wind from corrected data (*_l1.nc) 
measured during the PIANO measurement campaign. During the campaign regular  
conical scans (el = 70 deg) in ss mode with 24 beams distributed regularly on 
the cone were performed with the SL74, SL75 and SLXR142 lidar. This scan
pattern was performed once about every 20 min (three profiles per hour). 
The SL88 lidar performed continously the 6beam scan pattern (five beams in a 
cone with el = 70 deg and regularly distributed azimuth angles and a sixth beam
vertical, el = 90 deg). 
For each VAD24 scan, one vertical profile of horizontal wind is estimated. For
the SL88 lidar, 10 min averaged profiles are retrieved. The wind profiles are 
collected for one day and stored into daily .nc files.
Used directories:
    path_lidar_in           - input directory of *_l1.nc files; structure: path_lidar_in/yyyy/yyyymm/yyyymmdd/*_l1.nc
    path_lidar_out          - output directory of retrieved wind profiles: path_lidar_out/yyyy/yyyymmdd/*_vad.nc

Corrected .nc data files are created with netCDF_l1.py

To run this code, the module calc_vad.py is required. 
It can be found in the GITHub repositories:
    marenha/VAD_retrieval
@author: maren
"""
import os,sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from netCDF4 import Dataset
import time as time_pack

# Import own modules
path_parent = os.path.abspath('..')
path_VAD = os.path.join(path_parent,'doppler_wind_lidar_toolbox','VAD_retrieval')
path_2NetCDF = os.path.join(path_parent,'doppler_wind_lidar_toolbox','2NetCDF')
sys.path.append(path_VAD)
sys.path.append(path_2NetCDF)
import calc_vad as calc_vad
import piano_parameters as pp
import vad2NetCDF as v2n

#paramters
dtn = 24*60*60 # seconds per day

#%% path of level0 files and output folder
# path_data = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public')
path_data = os.path.join('P:\\','ZENODO_data_sets','Lidars')
path_l1 = os.path.join(path_data,'level1')
path_l2 = os.path.join(path_data,'level2')
# path_plot_out=os.path.join('/mnt','PIANO_campaign','lidars','quicklooks','vads')

#%% system id of lidars corresponding to the folders the raw data is stored
lidars_str=['SLXR_142','SL_75','SL_74','SL_88']

lidars_info=pp.load_lidars_info()

#%% loop through lidars

for lidar_str in lidars_str:
    lidar_info = lidars_info[lidar_str]
    
    # create ouput directory
    path_l2_lidar_out = os.path.join(path_l2,'%s_vad_l2' %lidar_info.name)  
    if not os.path.exists(path_l2_lidar_out): os.makedirs(path_l2_lidar_out)    

    # input directory of level1 data
    path_l1_lidar_in = os.path.join(path_l1,'%s_l1' %lidar_info.name)
    
    #%% loop through all available level1 sub-/directories 
    if lidar_str in ['SLXR_142','SL_75','SL_74']:
        
        path_l1_lidar_yyyymmdd = sorted(os.listdir(path_l1_lidar_in))
        for yyyymmdd_temp in path_l1_lidar_yyyymmdd:
            path_temp = os.path.join(path_l1_lidar_in,yyyymmdd_temp)
            files_all_sorted = sorted(os.listdir(path_temp))   
            files_vad = [file for file in files_all_sorted if 'vad24' in file]
            
            #%% estimate vertical profile of horizontal wind for each single scan
            # retrievals are collected in a list and combined to daily files
            datenum_list,snr_list=[],[]
            ws_list,wd_list,u_list,v_list,w_list=[],[],[],[],[]
            rv_fluc_list = []
            an_list = []
            # loop through each vad file
            for file_temp in files_vad:
                path_temp = os.path.join(path_l1_lidar_in,yyyymmdd_temp,file_temp)
                
                with xr.open_dataset(path_temp,decode_times=False) as data_temp:
                    az_rad_temp = np.deg2rad(data_temp.azimuth.values)
                    el_temp = data_temp.elevation.values
                    el_rad_temp = np.deg2rad(el_temp)
                    rv_temp = data_temp.radial_velocity.values
                    snr_temp = data_temp.intensity.values-1
                    gn = rv_temp.shape[0]
                    an_list.append(az_rad_temp.shape[0])
                    
                    date_num_temp = data_temp.datenum_time.values
                    
                    el_deg = np.unique(np.round(el_temp))
                    
                    range_gate_length = np.unique(np.diff(data_temp.gate_centers.values))[0]
                    gz_temp = data_temp.gate_centers.values*np.sin(np.deg2rad(el_deg))
                    snr_ray = np.mean(snr_temp,axis=1)
    
                # estimate horizontal wind for each range gate
                ws_ray,wd_ray = np.full(gn,np.nan),np.full(gn,np.nan)
                u_ray,v_ray,w_ray = np.full(gn,np.nan),np.full(gn,np.nan),np.full(gn,np.nan)
                rv_mean_ray,rv_fluc_ray = np.full(gn,np.nan), np.full(gn,np.nan)  
                for gi in range(0,gn):  
                    rv_tempp = rv_temp[gi,:]
                    ind_nan = np.where(~np.isnan(rv_tempp))[0]
    
                    rn=len(ind_nan)
                    if rn/rv_tempp.size<.8: continue
                    
                    u_ray[gi],v_ray[gi],w_ray[gi],ws_ray[gi],wd_ray[gi],rv_fluc_ray[gi] = calc_vad.calc_vad_3d(rv_tempp[ind_nan],el_rad_temp[ind_nan],az_rad_temp[ind_nan])
                    
                    rv_mean_ray[gi] = np.nanmean(rv_tempp)
                    
                # colltect retrievals in daily lists
                datenum_list.append(date_num_temp[0])
                ws_list.append(ws_ray)
                wd_list.append(wd_ray)
                u_list.append(u_ray)
                v_list.append(v_ray)
                w_list.append(w_ray)
                snr_list.append(snr_ray)
                rv_fluc_list.append(rv_fluc_ray)
                
            # transform to multi dimensional numpy.ndarray
            ws_array = np.vstack(ws_list).T
            wd_array = np.vstack(wd_list).T
            u_array = np.vstack(u_list).T
            v_array = np.vstack(v_list).T
            w_array = np.vstack(w_list).T
            snr_array = np.vstack(snr_list).T
            an_array = np.array(an_list)
            dn_array = np.array(datenum_list)
            rv_fluc_array = np.vstack(rv_fluc_list).T
            
            #%% filter data
            '''
            data is filtered with an SNR threshold of 18.2 dB and the first 
            range gates are invalid
            '''
            snr_threshold = -18.2
            snr_db_array = 10*np.log10(snr_array)
            ws_f_array,wd_f_array,u_f_array,v_f_array,w_f_array=np.copy(ws_array),np.copy(wd_array),np.copy(u_array),np.copy(v_array),np.copy(w_array)
            ind_th=np.where((snr_db_array<snr_threshold) | (np.isnan(snr_db_array)))
            for var in [ws_f_array,wd_f_array,u_f_array,v_f_array,w_f_array]:
                if lidar_info.name=='SLXR142': 
                    var[:5,:]=np.nan
                else:
                    var[:3,:]=np.nan
                var[ind_th]=np.nan
        
            #%% create netCDF file
            vad_var_temp = v2n.vad(dn_array,gz_temp,u_f_array,v_f_array,w_f_array,ws_f_array,wd_f_array,rv_fluc_array,snr_array,range_gate_length,snr_threshold,el_deg,an_array,u_array,v_array)
            v2n.to_netcdf(lidar_info, vad_var_temp, path_l2_lidar_out)
            

    #%% Retrievals for 6beam scans
    elif lidar_str in ['SL_88']:
        files_all_sorted = sorted(os.listdir(path_l1_lidar_in))
        files_vad = [file for file in files_all_sorted if '6beam' in file]

        for file_temp in files_vad:
            path_temp = os.path.join(path_l1_lidar_in,file_temp)
            yyyymmdd_temp = file_temp.split('_')[-1].split('.')[0]
            date_num = mdates.datestr2num(yyyymmdd_temp)
            
            with xr.open_dataset(path_temp,decode_times=False) as data_temp:
                az_temp = data_temp.azimuth.values
                el_temp = data_temp.elevation.values
                
                rv_temp = data_temp.radial_velocity.values
                snr_temp = data_temp.intensity.values-1
                date_num_temp = data_temp.datenum_time.values
                
                range_gate_length = np.unique(np.diff(data_temp.gate_centers.values))[0]
                
                gc_temp = data_temp.gate_centers.values
                
            # use only data with el = 90 deg
            ind_el = np.where(np.abs(el_temp - 70 ) <1)[0]
            az_temp = az_temp[ind_el]
            el_temp = el_temp[ind_el]
            date_num_temp = date_num_temp[ind_el]
            rv_temp = rv_temp[:,ind_el]
            snr_temp = snr_temp[:,ind_el]
            
            az_rad_temp = np.deg2rad(az_temp)
            el_rad_temp = np.deg2rad(el_temp)
            
            gn = rv_temp.shape[0]

            el_deg = np.unique(np.round(el_temp))

            gz_temp = gc_temp*np.sin(np.deg2rad(el_deg))

            dt_min = 10 # minutes per profile
            time_array = np.arange(date_num*dtn,(date_num+1)*dtn,dt_min*60)/dtn
            tn = time_array.size

            ws_array, wd_array = np.full([gn,tn],np.nan), np.full([gn,tn],np.nan)
            u_array, v_array, w_array = np.full([gn,tn],np.nan),np.full([gn,tn],np.nan), np.full([gn,tn],np.nan)
            snr_array = np.full([gn,tn],np.nan)
            rv_fluc_array, rv_mean_array = np.full([gn,tn],np.nan), np.full([gn,tn],np.nan)
            an_temp = np.full(tn, np.nan)
            for ti,tps in enumerate(time_array):  
                tpe = tps + 10/(24*60)
                ind_t = np.where((date_num_temp>tps) &(date_num_temp<=tpe))[0]
                
                # at least 10 beams have to be within the interval
                if len(ind_t) < 10: continue 
                
                snr_tempp = snr_temp[:,ind_t]
                snr_array[:,ti] = np.nanmean(snr_tempp,axis = 1)
                
                rv_tempp = rv_temp[:,ind_t]
                el_rad_tempp = el_rad_temp[ind_t]
                az_rad_tempp = az_rad_temp[ind_t]
                for gi in range(gn):  
                    rv_temppp = rv_tempp[gi,:]
                    ind_nan = np.where(~np.isnan(rv_temppp))[0]
    
                    rn=len(ind_nan)
                    if rn/rv_temppp.size<.8: continue
                    
                    u_array[gi,ti],v_array[gi,ti],w_array[gi,ti],ws_array[gi,ti],wd_array[gi,ti],rv_fluc_array[gi,ti] = calc_vad.calc_vad_3d(rv_temppp[ind_nan],el_rad_tempp[ind_nan],az_rad_tempp[ind_nan])
                    
                    rv_mean_array[gi,ti] = np.nanmean(rv_temppp)
                    
                an_temp[ti] = ind_t.size

            #%% filter data
            '''
            data is filtered with an SNR threshold of 18.2 dB and the first 
            range gates are invalid
            '''
            snr_threshold = -22
            snr_db_array = 10*np.log10(snr_array)
            ws_f_array,wd_f_array,u_f_array,v_f_array,w_f_array=np.copy(ws_array),np.copy(wd_array),np.copy(u_array),np.copy(v_array),np.copy(w_array)
            ind_th=np.where((snr_db_array<snr_threshold) | (np.isnan(snr_db_array)))
            for var in [ws_f_array,wd_f_array,u_f_array,v_f_array,w_f_array,rv_fluc_array]:
                var[:3,:]=np.nan
                var[ind_th]=np.nan
                
            #%% create netCDF file
            vad_var_temp = v2n.vad(time_array,gz_temp,u_f_array,v_f_array,w_f_array,ws_f_array,wd_f_array,rv_fluc_array,snr_array,range_gate_length,snr_threshold,el_deg,an_temp,u_array,v_array)
            
            v2n.to_netcdf(lidar_info, vad_var_temp, path_l2_lidar_out)

            
    