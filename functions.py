#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 13:10:19 2020

@author: nbayer
"""

import load_cams_data as ld
import numpy as np
import pandas as pd

def create_aerosol_file(lon,lat,date,saving_path,wvls=[469,550,670,865,1240]):
    """
    This function creates the aerosol file from CAMS for LibRadtran. It creates 
    the file for a specific latitude and longitude
    """
    
    """Define path+file where the CAMS.nc are located"""
    fsfc = "/vols/satellite/home/jonas/documents/paper/2020_clearsky_aerosoleffect/scripts/ecmwf/data/nc/cams-ra_"+date.strftime('%Y-%m-%d')+"_sfc.nc" 
    fml = "/vols/satellite/home/jonas/documents/paper/2020_clearsky_aerosoleffect/scripts/ecmwf/data/nc/cams-ra_"+date.strftime('%Y-%m-%d')+"_ml.nc"
    
    """load the CAMS data from the chosen netCDF""" 
    c = ld.CAMS(fsfc,fml)
    
    """Choose a Latitude to be used for the input_file for LibRadtran"""
    lats=c.lats.reshape((8,40,40))  # print(lats[0,:,0]) to see all lats
    lat_grid=[]
    lat_grid.extend(abs(lats[0,:,0]-lat))
    x=lat_grid.index(min(lat_grid))     # x is the position of the chosen lat

    """Choose a Longitude to be used for the input_file for LibRadtran"""
    lons=c.lons.reshape((8,40,40))  # print(lons[0,0,:]) to see all lons
    lon_grid=[]
    lon_grid.extend(abs(lons[0,0,:]-lon))
    y=lon_grid.index(min(lon_grid))     # x is the position of the chosen lon
    
    """load the aerosol data for the chosen wavelengths"""
    AP_sfc,AP_ml = c.aerosol_optprop(wvls)
    
    """crating DataFrame for the aerosol_file"""
    aerosol_file=pd.DataFrame({'z(km)':[],'aer_layer':[]})
        
    for z in range(0,60):
        layer_data=pd.DataFrame({'wavelength':[],'extintion coeffient [km-1]':[],'single scattering albedo':[],
                                 '0':[],'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],})
            
        for n in range(0,len(wvls)):
            layer_data=layer_data.append({'wavelength':wvls[n],
                                          'extintion coeffient [km-1]':('{0:.7}'.format(np.array(AP_ml.ext).reshape((8,40,40,60,len(wvls)))[4,x,y,z,n])),
                                          'single scattering albedo':('{0:.7}'.format(np.array(AP_ml.ssa).reshape((8,40,40,60,len(wvls)))[4,x,y,z,n])),
                                          '0':('{0:.4f}'.format(AP_ml.g.values.reshape((8,40,40,60,len(wvls)))[4,x,y,z,n]**0)),
                                          '1':('{0:.4f}'.format(AP_ml.g.values.reshape((8,40,40,60,len(wvls)))[4,x,y,z,n]**1)),
                                          '2':('{0:.4f}'.format(AP_ml.g.values.reshape((8,40,40,60,len(wvls)))[4,x,y,z,n]**2)),
                                          '3':('{0:.4f}'.format(AP_ml.g.values.reshape((8,40,40,60,len(wvls)))[4,x,y,z,n]**3)),
                                          '4':('{0:.4f}'.format(AP_ml.g.values.reshape((8,40,40,60,len(wvls)))[4,x,y,z,n]**4)),
                                          '5':('{0:.4f}'.format(AP_ml.g.values.reshape((8,40,40,60,len(wvls)))[4,x,y,z,n]**5)),
                                          '6':('{0:.4f}'.format(AP_ml.g.values.reshape((8,40,40,60,len(wvls)))[4,x,y,z,n]**6))
                                          },ignore_index=True)
        """Defining path + name for saving layer file """
        layer_file=saving_path+'Aerosol/layers/'+date.strftime('%y%m%d_')+'lat_lon-'+str(x)+'_'+str(y)+'-z'+str(z)+'.LAYER'    
        
        """saving layer file """
        layer_data.to_csv(layer_file, columns=['wavelength',
                                               'extintion coeffient [km-1]',
                                               'single scattering albedo',
                                               '0',
                                               '1',
                                               '2',
                                               '3',
                                               '4',
                                               '5',
                                               '6'],sep=' ', encoding='utf-8', header=False,index=False)
        """creating aerosol file"""
        aerosol_file=aerosol_file.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,40,40,60))[4,x,y,z]/1000)),
                                          'aer_layer': layer_file
                                          },ignore_index=True)
        """saving aerosol file """
        aerosol_file.to_csv(saving_path+'Aerosol/'+date.strftime('%y%m%d_')+'lat_lon-'+str(x)+'_'+str(y), 
                            columns=['z(km)','aer_layer'],
                            sep=' ', encoding='utf-8', header=False,index=False)


def create_atm_file(lon,lat,date,atmo_path):
    
    """Define path+file where the CAMS.nc are located"""
    fsfc = "/vols/satellite/home/jonas/documents/paper/2020_clearsky_aerosoleffect/scripts/ecmwf/data/nc/cams-ra_"+date.strftime('%Y-%m-%d')+"_sfc.nc" 
    fml = "/vols/satellite/home/jonas/documents/paper/2020_clearsky_aerosoleffect/scripts/ecmwf/data/nc/cams-ra_"+date.strftime('%Y-%m-%d')+"_ml.nc"
    
    """load the CAMS data from the chosen netCDF""" 
    c = ld.CAMS(fsfc,fml)
    SI=ld.SI

    """Choose a Latitude to be used for the input_file for LibRadtran"""
    lats=c.lats.reshape((8,40,40))  # print(lats[0,:,0]) to see all lats
    lat_grid=[]
    lat_grid.extend(abs(lats[0,:,0]-lat))
    x=lat_grid.index(min(lat_grid))     # x is the position of the chosen lat

    """Choose a Longitude to be used for the input_file for LibRadtran"""
    lons=c.lons.reshape((8,40,40))  # print(lons[0,0,:]) to see all lons
    lon_grid=[]
    lon_grid.extend(abs(lons[0,0,:]-lon))
    y=lon_grid.index(min(lon_grid))     # x is the position of the chosen lon

    """crating DataFrame for the atmospheric_file and gases_files"""
    data=pd.DataFrame({'z(km)':[],'p(mb)':[],'T(K)':[],'air(# * cm-3)':[]})
    ozone=pd.DataFrame({'z(km)':[],'O3 (mass mixing ratio) [kg kg^-1]':[]})
    no2=pd.DataFrame({'z(km)':[],'NO2 (mass mixing ratio) [kg kg^-1]':[]})
    h2o=pd.DataFrame({'z(km)':[],'H2O (mass mixing ratio) [kg kg^-1]':[]})
    
    for z in range(0,60):
        data=data.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[4,x,y,z]/1000)),
                  'p(mb)':('{0:.5f}'.format(c.P_mlvl.reshape((8,14,14,60))[4,x,y,z]/100)),
                  'T(K)':('{0:.3f}'.format(c.cams_ml.t[4,z,x,y].values)),
                  'air(# * cm-3)':('{0:.7}'.format((c.P_mlvl.reshape((8,14,14,60))[4,x,y,z]/1.e+6/c.cams_ml.t[4,z,x,y].values/SI.k)))
                  },ignore_index=True)
    
        ozone=ozone.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[4,x,y,z]/1000)),
                            'O3 (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.go3[4,z,x,y].values)
                            },ignore_index=True)
        
        no2=no2.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[4,x,y,z]/1000)),
                        'NO2 (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.no2[4,z,x,y].values)
                        },ignore_index=True)
        
        h2o=h2o.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[4,x,y,z]/1000)),
                        'H2O (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.q[4,z,x,y].values/(1-c.cams_ml.q[4,z,x,y].values))
                        },ignore_index=True)

    """ create path for saving the variables files"""
    atmo_file=atmo_path+date.strftime('%Y%m%d')+'_lat_lon:'+str(x)+'_'+str(y)+'.dat'
    mol_file_o3=atmo_path+date.strftime('%Y%m%d')+'_lat_lon:'+str(x)+'_'+str(y)+'-ozone.dat'
    mol_file_no2=atmo_path+date.strftime('%Y%m%d')+':_lat_lon:'+str(x)+'_'+str(y)+'-no2.dat'
    mol_file_h2o=atmo_path+date.strftime('%Y%m%d')+'_lat_lon:'+str(x)+'_'+str(y)+'-h2o.dat'
    
    """save the atmospheric data file as .dat"""
    data.to_csv(atmo_file, columns=['z(km)','p(mb)','T(K)','air(# * cm-3)'],
                sep=' ', encoding='utf-8', header=False,index=False)
    ozone.to_csv(mol_file_o3, columns=['z(km)','O3 (mass mixing ratio) [kg kg^-1]'],
                 sep=' ', encoding='utf-8', header=False,index=False)
    no2.to_csv(mol_file_no2, columns=['z(km)','NO2 (mass mixing ratio) [kg kg^-1]'],
               sep=' ', encoding='utf-8', header=False,index=False)
    h2o.to_csv(mol_file_h2o, columns=['z(km)','H2O (mass mixing ratio) [kg kg^-1]'],
               sep=' ', encoding='utf-8', header=False,index=False)

    
    
    
    
    
    
    
    