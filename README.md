# BTS 2048 data post-processing

[![License: CC BY-SA 4.0](https://licensebuttons.net/l/by-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-sa/4.0/)

The repository contains post-processing scripts to read and visualize solar radiation data exported by the  "Solarscan Software System" (hereinafter referred as Solarscan).


## Description
### General issue
The spectroradiometer BTS2048 is currently operated via a Solarscan. Solarscan is responsible for the a continuous scheduling and the measurement and the radiometric calibration of the BTS2048. Solarscan exports it's data in a specific ASCII-Format, see also [Solarscan User Manual](doc/Solarscan_BTS2048.pd).

The used file format of the Solarscan export is not well documented. However, the data should be archived in a well documented file format, which will be used of a majority of the Earth science community. The scripts of that repository will take account of this.

### Objectives
The software package is written to post-process observed solar radiation data of the array spectroradiometer BTS2018.
The software package includes
* a module to read OR0-Files from UV observations with the array spectroradiometer BTS2048 and
* a module to visualize the numerical data and
* a module to store data as netcdf file(s).


## Requirements

* <img src="https://www.python.org/static/community_logos/python-logo-generic.svg" alt="Python logo" style="width:130px;"> version??
 * argparse
 * calendar
 * datetime
 * matplotlib
 * NetCDF
 * numpy
 * os
 * pygal

## Usage
```
python BTS_main_process.py -s 20180913 -e 20181231
```

## Authors
* Rico Hengst
* Lionel Doppler
* Nicolas Bayer

## Cooperation
<img src="doc/TROPOS-Logo_ENG.svg" alt="TROPOS" style="width:100px;">  <img src="doc/Deutscherwetterdienst-logo.svg" alt="DWD" style="width:100px;">


## ToDo
* improve script to read, plot auxiliary data
