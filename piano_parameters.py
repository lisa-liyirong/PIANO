#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 17:24:45 2019
PIANO paramters
classes:
    lidar_location(lat,lon,zsl,name,lidar_id,bearing,gc_corr,pulse_frequency,start_str,end_str,diff_WGS84=np.nan,diff_geoid=np.nan,diff_bessel=np.nan)
    scan_type(name,lidars,scan_files,orientation)
functions: 
    load_lidars_info() - returns lidar_location() for SLXR_142, SL_88, SL_75 and SL_74
    load_scaninfo(scan_scenario) - scan_scenario = {'scenario1a', 'scenario1b', 'scenario2','scenario3'}
    load_grid(scan,delta_l) - scan = {'rhiew','rhisn','ppi3'}
@author: maren
"""
import sys,os
import numpy as np
import pyproj as ppj
import matplotlib.dates as mdates

path_parent = os.path.abspath('..')
path_cr = os.path.join(path_parent,'doppler_wind_lidar_toolbox','coplanar_retrieval')
sys.path.append(path_cr)
import calc_retrieval as cr

class lidar_location:     
    def __init__(self,lat,lon,zsl,name,lidar_id,bearing,gc_corr,pulse_frequency,start_str,end_str,diff_WGS84=np.nan,diff_geoid=np.nan,diff_bessel=np.nan):
        self.lat = lat
        self.lon = lon
        self.zsl = zsl #height in austrian levelled heights
        self.name = name
        self.lidar_id = lidar_id
        self.bearing = bearing # deviation from true geogrpahic north
        self.gc_corr = gc_corr
        self.pulse_frequency = pulse_frequency
        
        # define time period for which these settings are valid
        self.start_str = start_str
        self.end_num = end_str
        self.start_num = mdates.datestr2num(start_str)
        self.end_num = mdates.datestr2num(end_str)
        
        self.diff_WGS84 = diff_WGS84
        self.diff_geoid = diff_geoid
        self.diff_bessel = diff_bessel
        
        proj_tiris = ppj.Proj(init='epsg:31254')
        proj_wgs84 = ppj.Proj(init='epsg:4326')
        
        self.x,self.y,alt_dummy = ppj.transform(proj_wgs84,proj_tiris,lon,lat,zsl)
    
class scan_type():
    def __init__(self,name,lidars,scan_files,orientation):
        self.name = name
        self.lidars = lidars
        self.scan_files = scan_files
        self.orientation = orientation
        
def load_lidars_info():
    lidars=dict()
    
    lidars['SLXR_142'] = lidar_location(47.2639112,11.3846693,619.64 ,'SLXR142',142,-0.23,-18,10000,'20171007','20180308', diff_WGS84=49.17, diff_geoid=-0.18, diff_bessel=0.92)
    lidars['SL_74'] = lidar_location(47.2660760,11.4011756,629.13,'SL74',74,-0.56,0,15000,'20171007','20171218', diff_WGS84=49.12, diff_geoid=-0.18, diff_bessel=0.88)
    lidars['SL_75']  = lidar_location(47.2533684,11.3928383,623.79,'SL75',75,2.3,0,15000,'20171007','20171218', diff_WGS84=49.19, diff_geoid=-0.18, diff_bessel=0.934)
    lidars['SL_88'] = lidar_location(47.2649805,11.3893394,584.88,'SL88',88,0,0,15000,'20170719','20180308', diff_WGS84=49.15,diff_geoid=-0.18, diff_bessel=0.91)
    
    return lidars


def load_scaninfo(scan_scenario):
    if scan_scenario == 'scenario1b':
        rhisn = scan_type('rhiew',['SLXR_142','SL_74'],['csm_rhi_78.30_0.00-45.00_16x_s1_w500','csm_rhi_258.30_10.00-100.00_16x_s3_w500','csm_rhi_258.30_0.00-90.00_16x_s3_w500'],'vertical')
        rhiew = scan_type('rhisn',['SLXR_142','SL_75'],['csm_rhi_151.40_0.00-45.00_16x_s1_w500','csm_rhi_331.40_10.00-100.00_16x_s3_w500','csm_rhi_331.40_0.00-90.00_16x_s3_w500'],'vertical')
        ppi3 = scan_type('ppi3',['SLXR_142','SL_75','SL_74'],['csm_ppi_1.60_70.00-160.00_16x_s3_w500','csm_ppi_1.00_320.00-410.00_16x_s3_w500','csm_ppi_0.50_180.00-270.00_16x_s3_w500'],'horizontal')
        vad24 = scan_type('vad24',['SLXR_142','SL_75','SL_74'],['ss_vad_70.00_24rays_1x'],'conical')
        
        scan_types=[rhisn,rhiew,ppi3,vad24]
        
    elif scan_scenario == 'scenario3':
        rhi = scan_type('rhi',['SLXR_142'],['csm_rhi_258.30_0.00-45.00_32x_s1_w500','csm_rhi_151.40_0.00-45.00_32x_s1_w500','csm_rhi_78.30_0.00-45.00_32x_s1_w500'],'vertical')
        vad24 = scan_type('vad24',['SLXR_142'],['ss_vad_70.00_24rays_1x'],'conical')
        
        scan_types=[rhi,vad24]
        
    elif scan_scenario == 'scenario1a':
        rhisn = scan_type('rhiew',['SLXR_142','SL_74'],['csm_rhi_78.30_0.00-40.00_24x_s1_w500','csm_rhi_258.30_0.00-90.00_30x_s3_w500','csm_rhi_78.30_0.00-45.00_33x_s1_w1000','csm_rhi_258.30_0.00-90.00_33x_s3_w1000'],'vertical')
        rhiew = scan_type('rhisn',['SLXR_142','SL_75'],['csm_rhi_151.40_0.00-60.00_16x_s1_w500','csm_rhi_331.40_0.00-90.00_33x_s3_w1000','csm_rhi_151.40_0.00-45.00_33x_s1_w1000','csm_rhi_331.40_0.00-90.00_33x_s3_w500'],'vertical')
        ppi3 = scan_type('ppi3',['SLXR_142','SL_75','SL_74'],['csm_ppi_1.60_25.00-205.00_10x_s2_w500','csm_ppi_0.50_180.00-270.00_33x_s3_w1000','csm_ppi_1.00_320.00-410.00_33x_s3_w1000','csm_ppi_1.60_70.00-160.00_33x_s3_w1000','csm_ppi_1.00_320.00-410.00_33x_s3_w500'],'horizontal')
        vad24 = scan_type('vad24',['SLXR_142','SL_75','SL_74'],['ss_vad_70.00_24rays_1x'],'conical')
        
        scan_types=[rhisn,rhiew,ppi3,vad24]
        
    elif scan_scenario == 'scenario2':
        rhisn = scan_type('rhiew3',['SLXR_142','SL_75','SL_74'],['csm_rhi_78.30_0.00-40.00_24x_s1_w500','csm_rhi_348.30_0.00-120.00_24x_s3_w500','csm_rhi_258.30_0.00-120.00_24x_s3_w500'],'vertical')
        rhiew = scan_type('rhisn3',['SLXR_142','SL_75','SL_74'],['csm_rhi_151.40_0.00-60.00_16x_s1_w500','csm_rhi_331.40_0.00-180.00_16x_s3_w500','csm_rhi_241.40_0.00-180.00_16x_s3_w500'],'vertical')
        ppi3 = scan_type('ppi3',['SLXR_142','SL_75','SL_74'],['csm_ppi_1.60_25.00-205.00_10x_s2_w500','csm_ppi_1.00_270.00-450.00_10x_s2_w500','csm_ppi_0.50_150.00-330.00_10x_s2_w500'],'horizontal')
        vad24 = scan_type('vad24',['SLXR_142','SL_75','SL_74'],['ss_vad_70.00_24rays_1x'],'conical')
        
        scan_types=[rhisn,rhiew,ppi3,vad24]
        
    return scan_types


'''
Define grid for specific scan pattern (here: rhisn, rhiew and ppi3) 
and lattice length delta_l
input:
    scan (string)   - identifier of scan pattern
    delta_l (float) - lattice length of grid
return:
    calc_retrieval.grid
The grid is defined between lidar locations. Therefore, one lidar is put in 
the horizontal centre (x=0, y=0) of a metric spherical coordinate system (x,y,z) 
and the other lidar locations are defined relative to this centre. 
'''
def load_grid(scan,delta_l):
    lidars_info = load_lidars_info()
    # central lidar
    origin_x, origin_y = lidars_info['SLXR_142'].x, lidars_info['SLXR_142'].y
    
    delta_z = 2000 # vertical extent of the vertically orientated planes 
    

    if scan == 'rhisn': # vertical grid 
        lidars = ['SLXR_142','SL_75']
        for lidar in lidars:
            loc = lidars_info[lidar]
            loc.x, loc.y = loc.x-origin_x, loc.y-origin_y
            
        # distance between the two lidars
        delta_x = np.diff([lidars_info[lidar].x for lidar in lidars])
        delta_y = np.diff([lidars_info[lidar].y for lidar in lidars])
        distance = np.sqrt(delta_x**2+delta_y**2)+500
        
        # lower edge of vertical grid
        z_min = np.floor(np.min([lidars_info[lidar].zsl for lidar in lidars])/10)*10
        
        # angle between the two lidars
        alpha = np.arctan2(delta_y,delta_x) 
        
        # define 2D grid in 3D space
        xy = np.arange(0,distance,delta_l)
        x = xy*np.cos(alpha)
        y = xy*np.sin(alpha)
        z = np.arange(z_min,z_min+delta_z,delta_l)
      
     
    elif scan == 'rhiew': # vertical grid
        lidars = ['SLXR_142','SL_74']
        for lidar in lidars:
            loc = lidars_info[lidar]
            loc.x, loc.y = loc.x-origin_x, loc.y-origin_y
          
        # distance between the two lidars
        delta_x = np.diff([lidars_info[lidar].x for lidar in lidars])
        delta_y = np.diff([lidars_info[lidar].y for lidar in lidars])
        distance = np.sqrt(delta_x**2+delta_y**2)+500
        
        # lower edge of vertical grid
        z_min = np.floor(np.min([lidars_info[lidar].zsl for lidar in lidars])/10)*10
        
        # angle between the two lidars
        alpha = np.arctan2(delta_y,delta_x) # angle between x axis and silo location
        
        # define 2D grid in 3D space
        xy = np.arange(0,distance,delta_l)
        x = xy*np.cos(alpha)
        y = xy*np.sin(alpha)
        z = np.arange(z_min,z_min+delta_z,delta_l)
        
    
    elif scan == 'ppi3': # horizontal grid
        lidars = ['SLXR_142','SL_75','SL_74']
        for lidar in lidars:
            loc = lidars_info[lidar]
            loc.x,loc.y = loc.x-origin_x,loc.y-origin_y
            
        # define edges of horizontal grid
        x_min = np.floor(np.min([lidars_info[lidar].x for lidar in lidars])/10)*10
        x_max = np.ceil(np.max([lidars_info[lidar].x for lidar in lidars])/10)*10+200
        y_min = np.floor(np.min([lidars_info[lidar].y for lidar in lidars])/10)*10-200
        y_max = np.ceil(np.max([lidars_info[lidar].y for lidar in lidars])/10)*10+200
        
        x = np.arange(x_min,x_max,delta_l)
        y = np.arange(y_min,y_max,delta_l)
        z = np.array([0])
        
    # create grid for calc_retrieval.calc_retrieval(scan_list,grid,weight=None)
    grid=cr.grid(x,y,z,delta_l)
    
    return grid