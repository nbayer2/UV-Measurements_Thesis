#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 13:30:44 2019

@author: bayer
"""
import netCDF4 as nc4
from netCDF4 import date2num,num2date
import numpy as np
from datetime import datetime
# Hartwigs sunpos routine
from trosat import sunpos as sp
import logging

# create logger, name important
module_logger = logging.getLogger('uv-processing.BTS2NetCDF')

def netCDF_file(d_bts1day,nc_file, cfjson ):

    # """checking if the file does already exist, and delete if so."""
    # if os.path.isfile(netCDF_path+i8date+'.nc'):
    #     os.remove(netCDF_path+i8date+'.nc')
    
    ## calculate cosine of zenith angle
    szen, sazi = sp.sun_angles(d_bts1day['datetime'],cfjson['variables']['lat']['data'], cfjson['variables']['lon']['data'])
    esd=sp.earth_sun_distance(d_bts1day['datetime'])

    """ write dict data to numpy array https://www.quora.com/Is-it-possible-to-convert-a-Python-set-and-dictionary-to-a-NumPy-array """
    np_datetime = np.array(list(d_bts1day["datetime"]))
    
    """ prepare time variable to store seconds since ... in netcdf file"""
    second_since = date2num(np_datetime, cfjson['variables']['time']['attributes']['units'])
    
    
    """creating the netCDF file"""    
    f = nc4.Dataset(nc_file,'w', format='NETCDF4') #'w' stands for write 

    """Creating the dimensions in the NetCDF file from the JSON file"""
    for name,shape in cfjson['dimensions'].items(): 
        f.createDimension(name,len(d_bts1day[name]))
        # print(f.dimensions)
        
    """Building variables    time = uv_grp.createVariable('Time', 'i4', 'time')"""
    for variable in cfjson['variables']:
        f.createVariable(variable,cfjson['variables'][variable]['type'],(cfjson['variables'][variable]['shape'])) #,zlib=True,least_significant_digit=3)
        for attr in cfjson['variables'][variable]['attributes']:
            setattr(f[variable],attr,cfjson['variables'][variable]['attributes'][attr])
            # f[variable].attr=cfjson['variables'][variable]['attributes'][attr]
        # f[variable].units=cfjson['variables'][variable]['attributes']['units']
        """Adding units to the variables"""
        if variable=='time':
            f[variable][:]= second_since
        elif variable=='lon' or variable=='lat' or variable=='elev':
            f[variable][:]= cfjson['variables'][variable]['data']
        elif variable=='uva_irrad':
            f[variable][:]=d_bts1day['uva']
        elif variable=='uvb_irrad':
            f[variable][:]=d_bts1day['uvb']        
        elif variable=='uv_index':
            f[variable][:]=d_bts1day['uvind']        
        elif variable=='uv_irrad':
            f[variable][:]=d_bts1day['uvint']
        elif variable=='spec_irrad':
            f[variable][:,:]=d_bts1day['spect']
        elif variable=='wvl':
            f[variable][:]=d_bts1day['wvl']
        elif variable=='esd':
            f[variable][:]=esd
        elif variable=='szen':
            f[variable][:]=szen
        elif variable=='sazi':
            f[variable][:]=sazi
        else:
            print('This variable was not saved in the netCDF-File:',variable)
            continue
        # # f[x].setncatts(cfjson['variables'][x]['attributes'])
        # if variable=='time':
        #     f[variable][:]= second_since
        # elif cfjson['variables'][variable]['shape']==['time']:
        #     f[variable][:]=d_bts1day[variable]
        # elif len(cfjson['variables'][variable]['shape'])==2:
        #     f[variable][:,:]=d_bts1day[variable]
        # elif cfjson['variables'][variable]['shape']==['wvl']:
        #     f[variable][:]=d_bts1day[variable]
    today = datetime.today()
    f.history = "Created " + today.strftime("%d/%m/%y")

    """adding attributer"""
    for atributes in cfjson['attributes']:
        setattr(f,atributes,cfjson['attributes'][atributes])
    f.close()

    module_logger.info( 'Write data to netcdf: ' + nc_file )  
 
    
# """read netCDF"""
# nc = nc4.Dataset('/vols/satellite/home/bayer/uv/netCDF/2019/08/bts2048_uv_20190826_mp_l1_c1.nc','r')
# for i in nc.variables:
#     print(i,nc.variables[i] )#, nc.variables[i].units, nc.variables[i].shape)
#     print(nc.variables)
#     i=nc.variables[i][:]
#     #print(i)
# print(nc.history)
# print(nc.variables["uva_irrad"][:])
# print(nc)
    
# d_bts1day=bts.read_oro_bts('/vols/satellite/datasets/surfrad/uv_radiometer/raw/2019/08/26/MP190826.OR0','globals',20190826)
# nc = nc4.Dataset('/vols/satellite/home/bayer/uv/netCDF/2019/08/bts2048_uv_20190826_mp_l1_c1.nc','r')
# with open('/vols/satellite/home/bayer/uv/src/config/templates/uv_js_meta.json') as f:
#     cfjson= json.load(f)
