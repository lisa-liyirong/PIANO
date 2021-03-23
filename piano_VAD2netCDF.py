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
sys.path.append(path_VAD)
import calc_vad as calc_vad
import piano_parameters as pp

#paramters
dtn = 24*60*60 # seconds per day

#%% path of level0 files and output folder
# path_data = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public')
path_data = os.path.join('D:\\','ACINN_DL','DL_PIANO_public','Data')
path_l1 = os.path.join(path_data,'level1')
path_l2 = os.path.join(path_data,'level2')
# path_plot_out=os.path.join('/mnt','PIANO_campaign','lidars','quicklooks','vads')

#%% system id of lidars corresponding to the folders the raw data is stored
lidars_str=['SLXR_142','SL_75','SL_74','SL_88']
lidars_str = ['SL_88']

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
            datenum_list,snr_list,z_list=[],[],[]
            ws_list,wd_list,u_list,v_list,w_list=[],[],[],[],[]

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
                    an = az_rad_temp.shape[0]
                    
                    date_num_temp = data_temp.datenum_time.values
                    
                    el_deg = np.unique(np.round(el_temp))
                    
                    range_gate_length = np.unique(np.diff(data_temp.gate_centers.values))[0]
                    z = data_temp.gate_centers.values*np.sin(np.deg2rad(el_deg))
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
                z_list.append(z)
        # transform to multi dimensional numpy.ndarray
        ws_array = np.vstack(ws_list).T
        wd_array = np.vstack(wd_list).T
        u_array = np.vstack(u_list).T
        v_array = np.vstack(v_list).T
        w_array = np.vstack(w_list).T
        snr_array = np.vstack(snr_list).T
        
    #%% Retrievals for 6beam scans
    elif lidar_str in ['SL_88']:
        files_all_sorted = sorted(os.listdir(path_l1_lidar_in))
        files_vad = [file for file in files_all_sorted if '6beam' in file]

        for file_temp in files_vad:
            path_temp = os.path.join(path_l1_lidar_in,file_temp)
            date_str = file_temp.split('_')[-1].split('.')[0]
            date_num = mdates.datestr2num(date_str)
            
            with xr.open_dataset(path_temp,decode_times=False) as data_temp:
                az_temp = data_temp.azimuth.values
                el_temp = data_temp.elevation.values
                
                rv_temp = data_temp.radial_velocity.values
                snr_temp = data_temp.intensity.values-1
                date_num_temp = data_temp.datenum_time.values
                
                range_gate_length = np.unique(np.diff(data_temp.gate_centers.values))[0]
                
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
            an = az_rad_temp.shape[0]
            
            el_deg = np.unique(np.round(el_temp))

            z = data_temp.gate_centers.values*np.sin(np.deg2rad(el_deg))
            snr_ray = np.mean(snr_temp,axis=1)
            
            
            dt_min = 10 # minutes per profile
            time_array = np.arange(date_num*dtn,(date_num+1)*dtn,dt_min)/dtn
            tn = time_array.size
            
             
            
            ws_array, wd_array = np.full([gn,tn],np.nan), np.full([gn,tn],np.nan)
            u_array, v_array, w_array = np.full([gn,tn],np.nan),np.full([gn,tn],np.nan), np.full([gn,tn],np.nan)
            snr_array = np.full([gn,tn],np.nan)
            rv_fluc_array, rv_mean_array = np.full([gn,tn],np.nan), np.full([gn,tn],np.nan)
            for ti,tps in enumerate(time_array):  
                tpe = tps + 10/(24*60)
                ind_t = np.where((date_num_temp>tps) &(date_num_temp<=tpe))[0]
                
                # at least 10 beams have to be within the interval
                if len(ind_t) < 10: continue 
                
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
                
    
        #%% filter data
        '''
        data is filtered with an SNR threshold of 18.2 dB and the first 
        range gates are invalid
        '''
        snr_db_array = 10*np.log10(snr_array)
        ws_f_array,wd_f_array,u_f_array,v_f_array,w_f_array=np.copy(ws_array),np.copy(wd_array),np.copy(u_array),np.copy(v_array),np.copy(w_array)
        ind_th=np.where((snr_db_array<-18.2) | (np.isnan(snr_db_array)))
        for var in [ws_f_array,wd_f_array,u_f_array,v_f_array,w_f_array]:
            if lidar_info.name=='SLXR142': 
                var[:5,:]=np.nan
            else:
                var[:3,:]=np.nan
            var[ind_th]=np.nan
                
            
        #%% create netCDF file
        datenum_array=np.array(datenum_list)
        tn=datenum_array.shape[0]
            
            
        file_name = '%s_vad_%s.nc' %(lidar_info.name,yyyymmdd_temp)
        file_path = os.path.join(path_l2_lidar_out,file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)
            
        dataset_temp = Dataset(file_path,'w',format ='NETCDF4')
        # define dimensions
        dataset_temp.createDimension('NUMBER_OF_GATES',gn)
        dataset_temp.createDimension('NUMBER_OF_SCANS',tn)
        dataset_temp.createDimension('STATION_KEY',1)
    #    dataset_temp.createDimension('TEXT',1)
    
        # Metadata
        dataset_temp.description = "Profiles of horizontal wind speed and direction"
        dataset_temp.institution = "Department of Atmospheric and Cryospheric sciences (ACINN), University of Innsbruck, AUSTRIA",
        dataset_temp.contact = "Alexander Gohm (alexander.gohm@uibk.ac.at), Maren Haid, Lukas Lehner"
        dataset_temp.range_gate_length = '%i m' % range_gate_length
        dataset_temp.system_id = lidar_info.lidar_id
        dataset_temp.history = 'File created %s by M. Haid' %time_pack.strftime('%d %b %Y')
        dataset_temp.lat = lidar_info.lat
        dataset_temp.lon = lidar_info.lon
        dataset_temp.alt = lidar_info.zsl
        
        elevation=dataset_temp.createVariable('elevation',np.int64, ('STATION_KEY'))
        elevation.units = 'degrees'
        elevation.long_name = 'elevation'
        elevation.description = 'elevation of rays (number of rays specified in variable: rays) for VAD scan'
        elevation[:] = el_deg
        
        rays=dataset_temp.createVariable('rays',np.int64, ('STATION_KEY'))
        rays.units = 'unitless'
        rays.long_name = 'number of rays'
        rays.description = 'number of rays used to calculate mean wind profile (VAD algorithm) within the interval'
        rays[:] = an
        
        datenum = dataset_temp.createVariable('datenum',np.float64, ('NUMBER_OF_SCANS'))
        datenum.units = 'Number of days from January 0, 0000'
        datenum.long_name = 'start time of each conical scan'
        datenum.description = 'datenum (matlab) timestamp'
        datenum[:] = datenum_array
        
        unixtime=(datenum_array-mdates.datestr2num('19700101'))*(24*60*60)
        
        time = dataset_temp.createVariable('time',np.int64, ('NUMBER_OF_SCANS'))
        time.units = 'Seconds since 01-01-1970 00:00:00'
        time.long_name = 'start time of each conical scan'
        time.description = 'UNIX timestamp'
        time[:] = unixtime
        
        ff = dataset_temp.createVariable('ff', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ff.units = 'm s-1'
        ff.long_name  ='mean horizontal wind speed'
        ff.description = 'wind speed filtered with snr threshold of 0.015 (-18.2 dB)'
        ff[:,:] = ws_f_array
        
        dd = dataset_temp.createVariable('dd', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        dd.units = 'degrees'
        dd.long_name = 'mean horizontal wind direction'
        dd.description = 'wind direction filtered with snr threshold of 0.015 (-18.2 dB)'
        dd[:,:] = wd_f_array
        
        ucomp = dataset_temp.createVariable('ucomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ucomp.units = 'm s-1'
        ucomp.long_name = 'u component of horizontal wind vector'
        ucomp.description = 'u component filtered with snr threshold of 0.015 (-18.2 dB)'
        ucomp[:,:] = u_f_array
        
        vcomp = dataset_temp.createVariable('vcomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        vcomp.units = 'm s-1'
        vcomp.long_name = 'v component of horizontal wind vector'
        vcomp.description = 'v component filtered with snr threshold of 0.015 (-18.2 dB)'
        vcomp[:,:] = v_f_array
        
        wcomp = dataset_temp.createVariable('wcomp', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        wcomp.units = 'm s-1'
        wcomp.long_name = 'vertical velocity'
        wcomp.description = 'vertical compnent of three-dimensional wind vector filtered with snr threshold of 0.015 (-18.2 dB)'
        wcomp[:,:] = w_f_array
        
        ucomp_unfiltered = dataset_temp.createVariable('ucomp_unfiltered', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        ucomp_unfiltered.units = 'm s-1'
        ucomp_unfiltered.long_name = 'u component of horizontal wind vector'
        ucomp_unfiltered.description = 'non filtered u component'
        ucomp_unfiltered[:,:] = u_array
        
        vcomp_unfiltered = dataset_temp.createVariable('vcomp_unfiltered', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'))
        vcomp_unfiltered.units = 'm s-1'
        vcomp_unfiltered.long_name = 'v component of horizontal wind vector'
        vcomp_unfiltered.description = 'non filtered v component'
        vcomp_unfiltered[:,:] = v_array

        snr = dataset_temp.createVariable('snr', np.float32,('NUMBER_OF_GATES','NUMBER_OF_SCANS'),)
        snr.units = 'unitless'
        snr.long_name = 'signal to noise ratio (SNR)'
        snr.description = 'averaged profiles of snr; translation into dB: 10*log10(snr_array)'
        snr[:,:] = snr_array
        
        height = dataset_temp.createVariable('height',np.float32,('NUMBER_OF_GATES'))
        height.units = 'm'
        height.long_name = 'heigth of range gate centers'
        height.description = 'height of range gate centers above ground: gate_centers = (range_gate + 0.5) * range_gate_length * sin(elevation)'
        height[:] = z
            
        dataset_temp.close()
        