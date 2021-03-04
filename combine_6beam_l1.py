#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 09:55:14 2018
Combine 6beam data measured by SL88 to daily files
@author: maren
"""
import numpy as np
import xarray as xr
import os,sys 
from netCDF4 import Dataset
import piano_parameters as pp
import matplotlib.dates as mdates
dtn=24*60*60

#%% path of level0 files and output folder
path_l1 = os.path.join('/mnt','PIANO_campaign','lidars','netCDF_public','level1/')


#%% system id of lidars corresponding to the folders the raw data is stored
lidar_str='SL_88'

lidars_info=pp.load_lidars_info()
lidar_info = lidars_info[lidar_str]

path_l1_lidar_in=os.path.join(path_l1,'%s_l1' %lidar_info.name)
path_l1_lidar_out=os.path.join(path_l1,'%s_6beam_l1' %lidar_info.name)

# loop through all availale level0 sub-/directories 
path_l0_lidar_yyyymmdd = sorted(os.listdir(path_l1_lidar_in))

# loop through days --> one stare file per day
for yyyymmdd_temp in path_l0_lidar_yyyymmdd:
    path_temp=os.path.join(path_l1_lidar_in,yyyymmdd_temp)
    files_stare_sorted=sorted(os.listdir(path_temp))
    
    
    if mdates.datestr2num(yyyymmdd_temp) < mdates.datestr2num('20171126'): continue

    #%% Import and collect available data
    dwell_time_list,ds_list=[],[]
    for file in files_stare_sorted:
        file_path=os.path.join(path_temp,file)
        with xr.open_dataset(file_path,decode_times=False) as ds_temp:
             dwell_time_list.append(ds_temp.pulses_per_ray/lidar_info.pulse_frequency) # necessary to create uniform time array
             ds_list.append(ds_temp)                  

    #%% combine data to 2D arrays time x range_gate
    decimal_time_mat = np.hstack([ds_temp.decimal_time.values for ds_temp in ds_list])
    datenum_time_mat = np.hstack([ds_temp.datenum_time.values for ds_temp in ds_list])
    time_mat = np.hstack([ds_temp.time.values for ds_temp in ds_list])
    rv_mat = np.hstack([ds_temp.radial_velocity.values for ds_temp in ds_list])
    intensity_mat = np.hstack([ds_temp.intensity.values for ds_temp in ds_list])
    beta_mat = np.hstack([ds_temp.beta.values for ds_temp in ds_list])
    azimuth_mat = np.hstack([ds_temp.azimuth.values for ds_temp in ds_list])
    elevation_mat = np.hstack([ds_temp.elevation.values for ds_temp in ds_list])
    pitch_mat = np.hstack([ds_temp.pitch_angle.values for ds_temp in ds_list])
    roll_mat = np.hstack([ds_temp.roll_angle.values for ds_temp in ds_list])
    gc_mat = np.unique([ds_temp.gate_centers.values for ds_temp in ds_list])
    
    
    tn = decimal_time_mat.size
    gn = gc_mat.size
    
    #%% create netCDF files which have the same structure as other level1 data files
    file_out = '%s_6beam_l1_%s.nc' %(lidar_info.name,yyyymmdd_temp)
    file_out_path = os.path.join(path_l1_lidar_out,file_out)

    if os.path.isfile(file_out_path):
        os.remove(file_out_path)
            
    dataset_temp=Dataset(file_out_path,'w',format ='NETCDF4')
            
    # define dimensions
    dataset_temp.createDimension('NUMBER_OF_GATES',gn)
    dataset_temp.createDimension('NUMBER_OF_RAYS',tn)
    dataset_temp.createDimension('LOCATION',1)
    
    # Metadata: it is assumed that all scan files of the day have the same metadata
    dataset_temp.description = ds_temp.description + '; 6beam files combined for one day'
    dataset_temp.focus_range = ds_temp.focus_range
    dataset_temp.range_gate_length = ds_temp.range_gate_length
    dataset_temp.pulses_per_ray = ds_temp.pulses_per_ray
    dataset_temp.system_id = ds_temp.system_id
    dataset_temp.scan_type = ds_temp.scan_type
    dataset_temp.resolution = ds_temp.resolution
    
    
    latitude=dataset_temp.createVariable('lat', np.float32,('LOCATION'))
    latitude.units = 'decimal degree north'
    latitude.long_name = 'latitude'
    latitude.description = 'latitude, north is positive'
    latitude.missing_value = '-999.'
    latitude[:] = lidar_info.lat
    
    longitude=dataset_temp.createVariable('lon', np.float32,('LOCATION'))
    longitude.units = 'decimal degree east'
    longitude.long_name = 'latitude'
    longitude.description = 'latitude, east is positive'
    longitude.missing_value = '-999.'
    longitude[:] = lidar_info.lon  
    
    altitude=dataset_temp.createVariable('zsl', np.float32,('LOCATION'))
    altitude.units='m'
    altitude.long_name='altitude'
    altitude.missing_value = -999.
    altitude.difference_to_geoid = lidar_info.diff_geoid
    altitude.difference_to_geoid_descr = 'Difference between zsl and height above geoid for each STATION_KEY, e.g. zsl[1] + zsl.difference_to_geoid[1] = Height above Geoid in m'
    altitude.difference_to_bessel =  lidar_info.diff_bessel
    altitude.difference_to_bessel_descr = 'Difference between zsl and height above Bessel 1841 reference ellipsoid for geodetic datum MGI for each STATION_KEY; e.g. zsl[1] + zsl.difference_to_bessel[1] = Height above Bessel ellipsoid in m'
    altitude.difference_to_WGS84 = lidar_info.diff_WGS84
    altitude.difference_to_WGS_descr = 'Difference between zsl and height above WGS84 reference ellipsoid  for each STATION_KEY, e.g., zsl[1] + zsl.difference_to_WGS84[1] = Height above WGS84 ellispoid'
    altitude[:] = lidar_info.zsl
    
    bearing=dataset_temp.createVariable('bearing', np.float32,('LOCATION'))
    bearing.units = 'degrees'
    bearing.description = 'angle between the theoretical position of North direction in projected plane (epsg:31254) and North alignment of the instrument: bearing = North_lidar - North_projected '
    bearing.missing_value = -999.
    bearing[:] = lidar_info.bearing
    
    datenum_time = dataset_temp.createVariable('datenum_time', np.float64, ('NUMBER_OF_RAYS'))
    datenum_time.units = 'Days since 01-01-0001 00:00:00'
    datenum_time.long_name = 'start time number of each ray'
    datenum_time[:] = datenum_time_mat
    
    decimal_time=dataset_temp.createVariable('decimal_time', np.float64, ('NUMBER_OF_RAYS'))
    decimal_time.units = 'decimal time (hours) UTC'
    decimal_time.long_name = 'start time number of each ray'
    decimal_time[:] = decimal_time_mat
    
    time = dataset_temp.createVariable('time', np.float64, ('NUMBER_OF_RAYS'))
    time.units = 'Seconds since 01-01-1970 00:00:00'
    time.long_name = 'start time of each ray'
    time.description = 'UNIX timestamp'
    time[:] = time_mat
    
    azi = dataset_temp.createVariable('azimuth', np.float32, ('NUMBER_OF_RAYS'))
    azi.units = 'degrees, meteorologically'
    azi.long_name = 'azimuth angle'
    azi.comment = ds_temp.azimuth.comment
    azi[:] = azimuth_mat
    
    ele = dataset_temp.createVariable('elevation', np.float32, ('NUMBER_OF_RAYS'))
    ele.units = 'degrees, meteorologically'
    ele.long_name = 'elevation angle'
    ele[:] = elevation_mat
    
    pitch = dataset_temp.createVariable('pitch_angle', np.float32, ('NUMBER_OF_RAYS'))
    pitch.units = 'degrees'
    pitch.long_name = 'pitch angle'
    pitch[:] = pitch_mat
    
    roll = dataset_temp.createVariable('roll_angle', np.float32, ('NUMBER_OF_RAYS'))
    roll.units = 'degrees'
    roll.long_name = 'roll angle'
    roll[:] = roll_mat
    
    rv = dataset_temp.createVariable('radial_velocity', np.float32,('NUMBER_OF_GATES','NUMBER_OF_RAYS'))
    rv.units = 'm s-1'
    rv.long_name = 'Doppler velocity along line of sight'
    rv[:,:] = rv_mat

    intensity = dataset_temp.createVariable('intensity', np.float32,('NUMBER_OF_GATES','NUMBER_OF_RAYS'))
    intensity.units = 'unitless'
    intensity.long_name = 'SNR + 1'
    intensity[:,:] = intensity_mat
    
    beta = dataset_temp.createVariable('beta', np.float32,('NUMBER_OF_GATES','NUMBER_OF_RAYS'))
    beta.units = 'm-1 sr-1'
    beta.long_name = 'attenuated backscatter'
    beta[:,:] = beta_mat
    
    gate_centers = dataset_temp.createVariable('gate_centers',np.float32,('NUMBER_OF_GATES'))
    gate_centers.units = 'm'
    gate_centers.long_name = 'center of range gates'
    gate_centers[:] = gc_mat
    
    dataset_temp.close()