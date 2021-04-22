#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 11:07:02 2017
operations for measured halo data
- 
@author: maren
"""
from netCDF4 import Dataset
import xarray as xr
import numpy as np
import datetime as dt
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
import os
##
#import plot_dd_retrieval as pdd

'''
Algorithm finds start and end indixes of single Scans 
+ can be applied for both rhi and ppi scans
+ sort them into new attribute of xarray data frame (_scan ending)

'''
def separate_scans_l1(data):
    
    el=data.elevation.values
    az=data.azimuth.values

    if (az.max()-az.min())<0.1: #rhi
        angles_=el
    elif (el.max()-el.min())<0.1: #ppi
        angles_=az
    else:
        data=data.expand_dims('NUMBER_OF_SCANS',1)
        return data
    
    #find minimum and maximum angles as waypoints
    ag_min=np.round(np.min(angles_),2)
    ag_max=np.round(np.max(angles_),2)
    
    # mean angle schange
    ag_step=np.round(np.median(np.abs(np.diff(angles_))),2)
    
    scan_start_ind,scan_end_ind=np.array([],'int'),np.array([],'int')
    ag_act=angles_[0] #actual angle --> changes after new scan starts
    agi=0 #counter
    while agi<(angles_.shape[0]-5): #assumption: one scan is longer then 5 rays
        #start with one anlge
        scan_start_ind=np.append(scan_start_ind,agi)
        temp=[]
        # two options: ag_act is close to a)ag_max or b)ag_min
        # find next angle which has a distance to a)ag_min/b)ag_max smaller then 
        # ag_step/2
        # exception: if there is no other a)ag_max/b)ag_min following--> end of loop
        if np.abs(ag_act-ag_max)<np.abs(ag_act-ag_min): #a)closer to ag_max
            if ~np.any(np.abs(angles_[agi:]-ag_max)>ag_step):
                break
            temp=scan_start_ind[-1]+np.argwhere(np.abs(angles_[agi:]-ag_min)<ag_step/2)
        elif np.abs(ag_act-ag_max)>np.abs(ag_act-ag_min): #b)closer to ag_min
            if ~np.any(np.abs(angles_[agi:]-ag_min)>ag_step):
                break
            temp=scan_start_ind[-1]+np.argwhere(np.abs(angles_[agi:]-ag_max)<ag_step/2)
        
        if len(temp)==0:
            break
        scan_end_ind=np.append(scan_end_ind,temp[0])
        agi=scan_end_ind[-1]+1
        ag_act=angles_[agi]
    if len(temp)==0:
        scan_start_ind=scan_start_ind[:-1]
    number_of_scans=scan_start_ind.shape[0]
    number_of_rays_per_scan=int(np.max(scan_end_ind-scan_start_ind))
    
    # todo: über zwei Schleifen programmieren. Sieht hier super umständlich aus.
    #create names of coordiante information...
    coors=['x','y','z','xy']
    extends=['','_bottom','_top','_left','_right']
    systems=['_global','_local']
    coordinate_names=[]
    for coor in coors:
        for system in systems:
            for extend in extends:
                coordinate_names.append('%s%s%s' %(coor,system,extend))
    
    var_names=['radial_velocity','rv_filtered','intensity','beta']+coordinate_names
    for var_name in var_names:
        if not hasattr(data,var_name): #e.g. global coordinates are not available...
            continue
        var_temp=getattr(data,var_name).values
        temp=np.full([var_temp.shape[0],number_of_rays_per_scan,number_of_scans],np.nan)
        for si in range(0,number_of_scans):
            scan_start_temp=scan_start_ind[si]
            scan_end_temp=scan_end_ind[si]
            pulses_per_scan_temp=scan_end_temp-scan_start_temp
            
            temp[:,0:pulses_per_scan_temp,si]=var_temp[:,scan_start_temp:scan_end_temp]
            data['%s_scan' %var_name]=xr.DataArray(temp,dims=['NUMBER_OF_GATES','NUMBER_OF_RAYS_PER_SCAN','NUMBER_OF_SCANS'])
    
    #kann noch schöner gelöst werden, z.B. in obere Schleife integrieren
    az_scan=np.full([number_of_rays_per_scan,number_of_scans],np.nan)
    el_scan=np.full([number_of_rays_per_scan,number_of_scans],np.nan)
    start_scan=np.full(number_of_scans,np.nan)
    end_scan=np.full(number_of_scans,np.nan)
    for si in range(0,number_of_scans):
        scan_start_temp=scan_start_ind[si]
        scan_end_temp=scan_end_ind[si]
        pulses_per_scan_temp=scan_end_temp-scan_start_temp

        az_scan[0:pulses_per_scan_temp,si]=data.azimuth.values[scan_start_temp:scan_end_temp]
        el_scan[0:pulses_per_scan_temp,si]=data.elevation.values[scan_start_temp:scan_end_temp]
        
        start_scan[si]=data.decimal_time.values[scan_start_temp]
        end_scan[si]=data.decimal_time.values[scan_end_temp]
        

    data['el_scan']=xr.DataArray(el_scan,dims=['NUMBER_OF_RAYS_PER_SCAN','NUMBER_OF_SCANS'])
    data['az_scan']=xr.DataArray(az_scan,dims=['NUMBER_OF_RAYS_PER_SCAN','NUMBER_OF_SCANS'])
    data['start_scan']=xr.DataArray(start_scan,dims=['NUMBER_OF_SCANS'])
    data['end_scan']=xr.DataArray(end_scan,dims=['NUMBER_OF_SCANS'])
    data['start_scan_ind']=xr.DataArray(scan_start_ind,dims=['NUMBER_OF_SCANS'])
    data['end_scan_ind']=xr.DataArray(scan_end_ind,dims=['NUMBER_OF_SCANS'])
    
    return data

'''
Import processed data for time period times_scan_start 
'''
def import_data_combined(path_date_temp,times_scan_start,times_scan_end):
    files_all=os.listdir(path_date_temp)
   
    # scans of one scan pattern should be combined --> makes further work easier
    data_list=[]
    t_n=len(times_scan_start)
    for ti in range(0,t_n):
        if ti==0:
            file_temp=[file for file in files_all if (file.split('_')[0]=='User5')\
                       &(file.split('_')[3][0:4]>times_scan_start[ti])\
                       &(file.split('_')[3][0:4]<times_scan_end[ti])]
        else:
            file_temp=[file for file in files_all if (file.split('_')[0]=='User5')\
                       &(file.split('_')[3][0:4]>=times_scan_start[ti])\
                       &(file.split('_')[3][0:4]<times_scan_end[ti])]   
        if len(file_temp)==0:   
            continue
            
        file_temp=file_temp[0]
        path_temp=os.path.join(path_date_temp,file_temp)
        data_temp=xr.open_dataset(path_temp)
        data_list.append(data_temp)
        
        
    #combine both datasets along the number_of_rays variable
    data=xr.concat(data_list,dim='NUMBER_OF_RAYS')
    data['center_of_gate']=xr.DataArray(data.gate_centers.values[0,:],dims=['NUMBER_OF_GATES'], \
             attrs={'units':'m','long_name':'center of range gates 1d'})
    
    return data

'''
Center of measuremetn is shifted in the middle of the area for which radial velocity estimate is valid
+ borders of this area are given in local coordinate system
+ bearing can be corrected in case that hasn't be done before
'''
def correct_angles(lidar,data,correct_bearing=False):
    if np.max(np.abs(np.diff(data.azimuth.values)))>350:
        data.azimuth.values[data.azimuth.values<180]+=360
    
    if correct_bearing==True:
        data.azimuth.values=data.azimuth.values-lidar.bearing
    
    #angle is stored at the beginning of scanning segment
    #should represent the mass center of the averaging area
    angles=['elevation','azimuth']
    for angle in angles:
        angle_temp=np.copy(getattr(data,angle).values)
        angle_diff=np.diff(angle_temp)
        angle_diff=np.append(angle_diff,angle_diff[-1]) #assumption: scanner moving direction is the same for the last and forelast beam
        data[angle].values=angle_temp+angle_diff/2
        
    return data

def noise_reduce(lidar_id,rv,intensity,use_th=True):
    rv_=np.copy(rv)
    snr_=intensity-1
    snr_db=10*np.log10(snr_)
    rv_[np.isnan(snr_db)]=np.nan
    rv_all=np.full(rv_.shape,np.nan)
    
    gn=rv_.shape[0]
    tn=rv_.shape[1]

    if use_th==False:
        diff_rl=np.diff(rv_)[1:-1,:-1]
        diff_lr=np.fliplr(np.diff(np.fliplr(rv_)))[1:-1,1:]
        diff_ud=np.diff(rv_,axis=0)[:-1,1:-1]
        diff_du=np.flipud(np.diff(np.flipud(rv_),axis=0))[1:,1:-1]
        
        #borders will be ignored
        diff_mat=np.array([diff_rl,diff_lr,diff_ud,diff_du])
        diff_max_ind=np.argmax(np.abs(diff_mat),axis=0)
        diff_ind=np.indices((diff_mat.shape[1],diff_mat.shape[2]))
        #diff_mat(diff_max,np.indices((5, 5)))
        diff_mat[diff_max_ind.flatten(),diff_ind[0].flatten(),diff_ind[1].flatten()]=np.nan
    

        diff_sum=np.nansum(diff_mat,axis=0)
#        rv_var=np.nanmean(diff_mat,axis=0)
        
        diff_sum_median=np.median(diff_sum.flatten())
        diff_sum_percentile=np.percentile(diff_sum.flatten(),75)
        
        rv_all=np.full(rv_.shape,np.nan)
        rv_=rv_[1:-1,1:-1]
        
        #here threshold!!!!!! search for objective value
        rv_[np.abs(diff_sum)>(diff_sum_median+diff_sum_percentile)]=np.nan

        rv_all[1:-1,1:-1]=rv_
        
    elif use_th==True:
#        snr_threshold=10**(-17/10)+1
        snr_threshold=-20
        
        rv_[snr_db<snr_threshold]=np.nan
        rv_all=rv_
    
    #remove values which are surrounded by NaNs (these values are supposed to be NaNs as well)
    
    # for li in range(0,5):
    #     rv_nan=np.isnan(rv_all)
    #     for ri in range(1,rv_nan.shape[0]-1):
    #         for ci in range(1,rv_nan.shape[1]-1):
    #             sum_temp=np.sum([rv_nan[ri-1,ci],rv_nan[ri+1,ci],rv_nan[ri,ci-1],rv_nan[ri,ci+1],\
    #                         rv_nan[ri+1,ci-1],rv_nan[ri+1,ci+1],rv_nan[ri-1,ci-1],rv_nan[ri-1,ci+1]])
    #             if sum_temp>=5:
    #                 rv_all[ri,ci]=np.nan
                    
    rv_extend=np.full([gn-2,tn-2,9],np.nan)
    for li in range(5): #repeat procedure as often as required
        rv_extend[:,:,:]=np.nan
        rv_filtered_nan=np.isnan(rv_all)
        count=0
        for ci in range(3):
            for ri in range(3):
                re=gn-2+ri
                ce=tn-2+ci
                rv_extend[:,:,count]=rv_filtered_nan[ri:re,ci:ce]
                count+=1
           
        sum_rv_extend=np.sum(rv_extend,axis=2)
        #if measurement point is surroundet by more than 5 nan it's set on nan, too
        rv_all[1:-1,1:-1][sum_rv_extend>5]=np.nan 
                        
    if lidar_id=='SLXR142':
        rv_all[:5,:]=np.nan   
    elif lidar_id in ['SL75','SL74','SL88']:
        rv_all[:3,:]=np.nan 

    return rv_all

def decimal2datetime(date_str,decimal_time):
    
    if (type(decimal_time) is np.float64) or (type(decimal_time) is float):
        date_time=dt.datetime(int(date_str[0:4]),int(date_str[4:6]),int(date_str[6:8]),int(decimal_time),(np.floor((decimal_time*60) % 60)).astype(int),(np.floor((decimal_time*3600) % 60)).astype(int))
    else:   
        date_time=[dt.datetime(int(date_str[0:4]),int(date_str[4:6]),int(date_str[6:8]),x.astype(int),(np.floor((x*60) % 60)).astype(int),(np.floor((x*3600) % 60)).astype(int)) for x in decimal_time]

    return mdates.date2num(date_time)

def decimal2unix(date_str,decimal_time):
    date_num=mdates.datestr2num(date_str)
    if (type(decimal_time) is np.float64) or (type(decimal_time) is float):
        date_time=dt.datetime(int(date_str[0:4]),int(date_str[4:6]),int(date_str[6:8]),int(decimal_time),(np.floor((decimal_time*60) % 60)).astype(int),(np.floor((decimal_time*3600) % 60)).astype(int))
        unixtime=(date_time-dt.datetime(1970,1,1)).total_seconds()
    else: 
#        date_time=np.full(len(decimal_time),0,dtype='object')
#        for ti,x in enumerate(decimal_time):
#            x_num=date_num+(x/24)
#            hour=int(mdates.num2datestr(x_num,'%H'))
#            minute=int(mdates.num2datestr(x_num,'%M'))
#            second=int(mdates.num2datestr(x_num,'%S'))
#            date_time[ti]=dt.datetime(int(date_str[0:4]),int(date_str[4:6]),int(date_str[6:8]),hour,minute,second)
        date_time=[]
        for di,x in enumerate(decimal_time):
            if x==24.0:
                temp=dt.datetime(int(mdates.num2datestr(date_num+1,'%Y')),int(mdates.num2datestr(date_num+1,'%m')),int(mdates.num2datestr(date_num+1,'%d')))
            else:
                temp=dt.datetime(int(date_str[0:4]),int(date_str[4:6]),int(date_str[6:8]),x.astype(int),(np.floor((x*60) % 60)).astype(int),(np.floor((x*3600) % 60)).astype(int))
            
           
            date_time.append(temp)
                


#        date_time=[dt.datetime(int(date_str[0:4]),int(date_str[4:6]),int(date_str[6:8]),x.astype(int),(np.floor((x*60) % 60)).astype(int),(np.floor((x*3600) % 60)).astype(int)) for x in decimal_time]
        unixtime=[(d-dt.datetime(1970,1,1)).total_seconds() for d in date_time]
    
    
    return unixtime

'''
Corrdinates of measurements are shifted in a way that they will lay in an Cartesian 
coordinate system valid for all measurements
'''
def shift_coordinate_system(lidar_temp,data_temp,origin):
    
    name_extends=['_bottom','','_top','_left','_right'] 
    for name_extend in name_extends:
        x_temp=getattr(data_temp,'x_local%s' %name_extend).values+lidar_temp.x-origin.x
        y_temp=getattr(data_temp,'y_local%s' %name_extend).values+lidar_temp.y-origin.y
        z_temp=getattr(data_temp,'z_local%s' %name_extend).values+lidar_temp.zs
        
        xy_temp=np.sqrt(x_temp**2+y_temp**2)
        xy_temp[(x_temp<0)]=-xy_temp[(x_temp<0)]
#        xy_temp[(x_temp<0)|(y_temp<0)]=-xy_temp[(x_temp<0)|(y_temp<0)]
        
        data_temp['x_global%s' %name_extend]=xr.DataArray(x_temp,dims=['NUMBER_OF_GATES','NUMBER_OF_RAYS'], \
                 attrs={'units':'m','long_name':'global coordinate in east direction'})
        data_temp['y_global%s' %name_extend]=xr.DataArray(y_temp,dims=['NUMBER_OF_GATES','NUMBER_OF_RAYS'], \
                 attrs={'units':'m','long_name':'global coordinate in north direction'})
        data_temp['z_global%s' %name_extend]=xr.DataArray(z_temp,dims=['NUMBER_OF_GATES','NUMBER_OF_RAYS'], \
                 attrs={'units':'m','long_name':'measurement height in levelled height'})
        data_temp['xy_global%s' %name_extend]=xr.DataArray(xy_temp,dims=['NUMBER_OF_GATES','NUMBER_OF_RAYS'], \
                 attrs={'units':'m','long_name':'distance to origin'})
    
    return data_temp
