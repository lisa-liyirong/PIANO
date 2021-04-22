# PIANO
## Abstract 
Script collection for data processing of [HALO Photonics](https://halo-photonics.com/) Doppler wind lidars measurements operated during the [PIANO](https://www.uibk.ac.at/projects/piano/index.html.en) (Penetration and Interruption of Alpine Foehn) measurement campaign. From 18 Sep 2017 to 18 Dec 2017 (individual Lidar installation periods differ), four Doppler wind lidars, model Stream Line (SL_74, SL_75 and SL_88) and Stream Line XR (SLXR_142), were operated by the [ACINN](https://www.uibk.ac.at/acinn/index.html.en) (Department of Atmospheric and Cryospheric Sciences, University of Innsbruck) in the city of Innsbruck. The lidars have a rotating scanner head which allows for a variety of scan pattern. The Stream Line software provides profiles of radial velocity and Signal-to-Noise ratio stored in txt files.

## Level0 data
**Conversion of Halo Photonics StreamLine .hpl files into .nc data files.**
### `netCDF_l0.py`
The StreamLine software generates .hpl files (text files) containing time-distance resolved data of radial velocity, Signal-to-Noise (SNR) and attenuated backscatter (+ spectral width for newer software version). Additionally, for each time stamp, the scanner position is given as azimuth and elevation angle. This data is converted into matrices (time x distance) and converted into netCDF (.nc) level0 (l0) data (e.g., `hpl_to_netcdf()` in marenha/doppler_wind_lidar_toolbox/2NetCDF/hpl2NetCDF.py). Level 0 data contains the same data as the original .hpl files.

## Level1 data
**Generation of .nc files of corrected data.**
### `netCDF_l1.py`
In `netCDF_l1.py`, the data of the l0 .nc files is corrected and stored as level 1 (l1) netCDF data. The corrected variables are the azimuth angle and the range gate centre. The azimuth angle needs to be corrected due to a misalignment of the lidar scanner to geographical North in the field. Further, the range gates centres of the SLXR_142 show to have an offset which needs to be corrected. Additionally, information about the lidar location is added to the l1 netCDF files.

### `combine_stare_l1.py`
The StreamLine software stored hourly files of stare data. In `combine_stare_l1.py`, the hourly l1 .nc files are combined to daily files. 

### `combine_6beam_l1.py`
The SL_88 lidar performed 6beam scans during the complete campaign period. For each 6beam scan (about 26 seconds) one .hpl file is created. In `combine_6beam_l1.py`, the l1 .nc files are combined into daily l1 .nc files.

### `combine_scenario_l1.py`
For the estimation of two-dimensional wind fields with `netCDF_ret_l2.py` from multiple-Doppler lidar coplanar scans, the data is prepared. This code only works for scan patterns performed within the scan scenario *scenario1b* (details in `piano_paramters.py`). 

## Level2 data
### `piano_VAD2netCDF.py`
**Generation of daily .nc files of vertical profiles of horizontal wind.**

During the campaign, the SL_74, SL_75 and SLXR_142 performed regularly VAD scans at an elevation angle of 70 deg. The data of these scans is used to derive vertical profiles of horizontal wind using the VAD algorithm. For the SL_88, the 6beam scans are used to derive 10-minute averages of vertical profiles of horizontal wind. 
For each day one data file is created containing the retrieved wind profiles.

### `netCDF_ret_l2.py`
**Generation of .nc files of two-dimensional wind fields of differente scan planes.**

During the campaign, the SL_74, SL_75 and SLXR_142 performed synchronized coplanar scans for different planes. For these scans, the two-dimensional wind on the plane is derived. The core algorithm can be found in marenha/doppler_wind_lidar_toolbox/coplanar_retrieval/calc_retrieval.py. 

## Paramter files and help functions 
### `piano_paramters.py`
Parameter file containing information about lidar locations, lidar parameters, scan scenarios and scan pattern. Additionally, the function `load_grid()` creates a grid for the performed coplanar scan pattern. 

### `halo_data_operations.py`
Functions of HaloPhonics Doppler wind lidar treatment. The functions are PIANO campaign specific and cannot be used directly for other scan scenarios and pattern. This includes a filter function `noise_reduce()` and a angle correction function `correct_angles()`.
