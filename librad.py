#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:43:23 2020

@author: nbayer
"""
import load_cams_data as ld
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import netCDF4 as nc4
from netCDF4 import date2num,num2date
from datetime import datetime
import xarray as xr

"""
#Process for all latitudes and longitudes
python librad.py -d 2015-06-30 -wvls 280 310 340 370 400 469 500

#Process for a specific latitude and longitude
python librad.py -d 2015-06-30 -wvls 280 310 340 370 400 469 500 -lat 51.25 -lon 12.9

"""

"""for calling the function from the terminal"""
parser = argparse.ArgumentParser(description='Process to create and running input_files for LibRadTran from the CAMS_ra for a specific date')
parser.add_argument('-d', type=str, required=True, dest='date', # the variable is saved in args.date as a string
                    help='Insert the date as 2019-01-07(YYYY-MM-DD)')
parser.add_argument('-wvls', nargs='+',type=float, required=True, dest='wvls', # the variable is saved in args.wvls as a string
                    help='Insert the wavelengths for the calculations as wvls1 wvls2 ... wvlsn')
parser.add_argument('-lat', type=float, dest='lat',
                    help="for using a specific latitude")
parser.add_argument('-lon',  type=float, dest='lon', 
                    help="for using a specific longitude")
args = parser.parse_args()

"""Define path+file where the CAMS.nc are located"""
fsfc = "/vols/satellite/home/jonas/documents/paper/2020_clearsky_aerosoleffect/scripts/ecmwf/data/nc/cams-ra_"+args.date+"_sfc.nc" 
fml = "/vols/satellite/home/jonas/documents/paper/2020_clearsky_aerosoleffect/scripts/ecmwf/data/nc/cams-ra_"+args.date+"_ml.nc"                  

"""load the CAMS data from the chosen netCDF""" 
c = ld.CAMS(fsfc,fml)

"""load atmospheric and gases Constants from load_cams_data.py""" 
SI=ld.SI
constants=ld.CONSTANTS

#t_step=14*14 # every time step starts after t_step (grid points)
t_step=41*41
"""Choose a Latitude to be used for the input_file for LibRadtran"""
if args.lat:
    lat=51.5
    lats=c.lats.reshape((8,14,14))  # print(lats[0,:,0]) to see all lats
    lat_grid=[]
    lat_grid.extend(abs(lats[0,:,0]-args.lat))
    x=lat_grid.index(min(lat_grid))     # x is the position of the chosen lat

"""Choose a Longitude to be used for the input_file for LibRadtran"""
if args.lon:
    lon=12.93
    lons=c.lons.reshape((8,14,14))  # print(lons[0,0,:]) to see all lons
    lon_grid=[]
    lon_grid.extend(abs(lons[0,0,:]-args.lon))
    y=lon_grid.index(min(lon_grid))     # x is the position of the chosen lon

"""Choose the Wavelengths """
wvls=args.wvls

times=[]        #contains all the time steps ##times[1].strftime("%m/%d/%y")
for ti in range(0, len(c.times),t_step):
    times.append(pd.to_datetime(c.times[ti]))

"""Define the path to save the atmospheric_file"""
atmo_path='/vols/satellite/home/bayer/libradtran/Libradtran-files/atmospheric_files/'


"""load the aerosol data for the chosen wavelengths"""
AP_sfc,AP_ml = c.aerosol_optprop(wvls)
UVI=pd.DataFrame({'Date':[],'UVI':[]})
for t in range(3,5):
    if args.lat and args.lon:
        """crating DataFrame for the atmospheric_file and gases_files"""
        data=pd.DataFrame({'z(km)':[],'p(mb)':[],'T(K)':[],'air(# * cm-3)':[]})
        ozone=pd.DataFrame({'z(km)':[],'O3 (mass mixing ratio) [kg kg^-1]':[]})
        no2=pd.DataFrame({'z(km)':[],'NO2 (mass mixing ratio) [kg kg^-1]':[]})
        h2o=pd.DataFrame({'z(km)':[],'H2O (mass mixing ratio) [kg kg^-1]':[]})
    
        for z in range(0,60):
            data=data.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
                      'p(mb)':('{0:.5f}'.format(c.P_mlvl.reshape((8,14,14,60))[t,x,y,z]/100)),
                      'T(K)':('{0:.3f}'.format(c.cams_ml.t[t,z,x,y].values)),
                      'air(# * cm-3)':('{0:.7}'.format((c.P_mlvl.reshape((8,14,14,60))[t,x,y,z]/1.e+6/c.cams_ml.t[t,z,x,y].values/SI.k)))
                      },ignore_index=True)
    
            ozone=ozone.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
                                'O3 (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.go3[t,z,x,y].values)
                                },ignore_index=True)
        
            no2=no2.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
                            'NO2 (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.no2[t,z,x,y].values)
                            },ignore_index=True)
        
            h2o=h2o.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
                            'H2O (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.q[t,z,x,y].values/(1-c.cams_ml.q[t,z,x,y].values))
                            },ignore_index=True)

        """ create path for saving the variables files"""
        atmo_file=atmo_path+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'.dat'
        mol_file_o3=atmo_path+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'-ozone.dat'
        mol_file_no2=atmo_path+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'-no2.dat'
        mol_file_h2o=atmo_path+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'-h2o.dat'
        
        """save the atmospheric data file as .dat"""
        data.to_csv(atmo_file, columns=['z(km)','p(mb)','T(K)','air(# * cm-3)'],
                    sep=' ', encoding='utf-8', header=False,index=False)
        ozone.to_csv(mol_file_o3, columns=['z(km)','O3 (mass mixing ratio) [kg kg^-1]'],
                     sep=' ', encoding='utf-8', header=False,index=False)
        no2.to_csv(mol_file_no2, columns=['z(km)','NO2 (mass mixing ratio) [kg kg^-1]'],
                   sep=' ', encoding='utf-8', header=False,index=False)
        h2o.to_csv(mol_file_h2o, columns=['z(km)','H2O (mass mixing ratio) [kg kg^-1]'],
                   sep=' ', encoding='utf-8', header=False,index=False)

        """crating DataFrame for the aerosol_file"""
        aerosol_file=pd.DataFrame({'z(km)':[],'aer_layer':[]})
        for z in range(0,60,2):
            layer_data=pd.DataFrame({'wavelength':[],'extintion coeffient [km-1]':[],'single scattering albedo':[],
                                     '0':[],'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],})
            
            for n in range(0,len(wvls)):
                layer_data=layer_data.append({'wavelength':wvls[n],
                                              'extintion coeffient [km-1]':('{0:.7}'.format(np.array(AP_ml.ext).reshape((8,14,14,60,len(wvls)))[t,x,y,z,n])),
                                              'single scattering albedo':('{0:.7}'.format(np.array(AP_ml.ssa).reshape((8,14,14,60,len(wvls)))[t,x,y,z,n])),
                                              '0':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**0)),
                                              '1':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**1)),
                                              '2':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**2)),
                                              '3':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**3)),
                                              '4':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**4)),
                                              '5':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**5)),
                                              '6':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**6))
                                              },ignore_index=True)
                
            layer_file=atmo_path+'Aerosol/layers/'+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'-z'+str(z)+'.LAYER'    
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
    
            aerosol_file=aerosol_file.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
                                              'aer_layer': layer_file
                                              },ignore_index=True)
    
            aerosol_file.to_csv(atmo_path+'Aerosol/'+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y), 
                                columns=['z(km)','aer_layer'],
                                sep=' ', encoding='utf-8', header=False,index=False)
        """write the paths and input files for libradtran for each day and hs"""
        input_file="/vols/satellite/home/bayer/libradtran/Libradtran-files/input_files/"+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+".txt"
        output_file="/vols/satellite/home/bayer/libradtran/Libradtran-files/output_files/"+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+".txt"
        
        f=open(input_file, "w+")
        f.write
        f.write("\nday_of_year "+times[t].strftime('%j'))   # Correct for Earth-Sun distance
        f.write("\ndata_files_path /home/nbayer/libRadtran-2.0.3/data/")
        f.write("\natmosphere_file midlatitude_summer")
        #f.write("\natmosphere_file "+atmo_file)
        f.write("\nsource solar  ../solar_flux/kurudz_0.1nm.dat per_nm ")   #line identifies the location of the extraterrestrial solar flux file which defines the spectral resolution.
        # f.write("\nmol_file O3 "+mol_file_o3+ ' mmr')
        # f.write("\nmol_file NO2 "+mol_file_no2+ ' mmr')
        # f.write("\nmol_file H2O "+mol_file_h2o+ ' mmr')
        # f.write("\npressure "+str(c.cams_sfc.psfc[t,x,y].values/100))
        f.write("\nsza "+str(c.sza.reshape((8,14,14))[t,x,y]))
        # f.write("\nphi0 "+str(c.azi.reshape((8,14,14))[t,x,y]))
        f.write("\nmol_abs_param sbdart #spectralcalculation resolver, should be the best option for UV Indax calculations")
        f.write("\naerosol_default")    #switch the use of aerosol data on
        f.write("\naerosol_file explicit "+atmo_path+'Aerosol/'+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y))
        #f.write("\nck_lowtran_absorption O4 off")
        # f.write("\nrte_solver disort")
        # f.write("\nrte_solver fdisort2")
        # f.write("\nrte_solver twostr")
        # f.write("\nrte_solver rodents")
        f.write("\nrte_solver polradtran^")
        f.write("\ndisort_intcor moments")
        # f.write("\nno_absorption mol")
        # f.write("\nno_scattering mol")
        f.write("\nwavelength "+str(min(wvls))+" "+ str(max(wvls)))    #Wavelength range [nm]
        f.write("\noutput_process per_nm")
        f.write("\nverbose") 
        f.close()
        
        Librad_path='/home/nbayer/libRadtran-2.0.3/bin/uvspec'
        """Running LibRadTran with the input_file frim the step above and saving it in the output_file"""
        """run the script from /home/nbayer/libRadtran-2.0.3/bin/ """
        os.system(Librad_path+" < "+input_file+" > "+output_file)
    else:
        print('introduce lat and lon')
        break
        # for x in range(0,14):
        #     for y in range(0,14):
        #         """crating DataFrame for the atmospheric_file and gases_files"""
        #         data=pd.DataFrame({'z(km)':[],'p(mb)':[],'T(K)':[],'air(# * cm-3)':[]})
        #         ozone=pd.DataFrame({'z(km)':[],'O3 (mass mixing ratio) [kg kg^-1]':[]})
        #         no2=pd.DataFrame({'z(km)':[],'NO2 (mass mixing ratio) [kg kg^-1]':[]})
        #         h2o=pd.DataFrame({'z(km)':[],'H2O (mass mixing ratio) [kg kg^-1]':[]})
                
        #         for z in range(0,60):
        #             data=data.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
        #                               'p(mb)':('{0:.5f}'.format(c.P_mlvl.reshape((8,14,14,60))[t,x,y,z]/100)),
        #                               'T(K)':('{0:.3f}'.format(c.cams_ml.t[t,z,x,y].values)),
        #                               'air(# * cm-3)':('{0:.7}'.format((c.P_mlvl.reshape((8,14,14,60))[t,x,y,z]/1.e+6/c.cams_ml.t[t,z,x,y].values/SI.k)))
        #                               },ignore_index=True)
        
        #             ozone=ozone.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
        #                                 'O3 (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.go3[t,z,x,y].values)
        #                                 },ignore_index=True)
                    
        #             no2=no2.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
        #                             'NO2 (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.no2[t,z,x,y].values)
        #                             },ignore_index=True)
                    
        #             h2o=h2o.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
        #                             'H2O (mass mixing ratio) [kg kg^-1]':'{0:.7}'.format(c.cams_ml.q[t,z,x,y].values/(1-c.cams_ml.q[t,z,x,y].values))
        #                             },ignore_index=True)

        #         """ create path for saving the variables files"""
        #         atmo_file=atmo_path+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'.dat'
        #         mol_file_o3=atmo_path+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'-ozone.dat'
        #         mol_file_no2=atmo_path+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'-no2.dat'
        #         mol_file_h2o=atmo_path+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'-h2o.dat'
                
        #         """save the atmospheric data file as .dat"""
        #         data.to_csv(atmo_file, columns=['z(km)','p(mb)','T(K)','air(# * cm-3)'],
        #                     sep=' ', encoding='utf-8', header=False,index=False)
        #         ozone.to_csv(mol_file_o3, columns=['z(km)','O3 (mass mixing ratio) [kg kg^-1]'],
        #                      sep=' ', encoding='utf-8', header=False,index=False)
        #         no2.to_csv(mol_file_no2, columns=['z(km)','NO2 (mass mixing ratio) [kg kg^-1]'],
        #                    sep=' ', encoding='utf-8', header=False,index=False)
        #         h2o.to_csv(mol_file_h2o, columns=['z(km)','H2O (mass mixing ratio) [kg kg^-1]'],
        #                    sep=' ', encoding='utf-8', header=False,index=False)
                
        #         """crating DataFrame for the aerosol_file"""
        #         aerosol_file=pd.DataFrame({'z(km)':[],'aer_layer':[]})
        #         for z in range(0,60,2):
        #             layer_data=pd.DataFrame({'wavelength':[],'extintion coeffient [km-1]':[],'single scattering albedo':[],
        #                                      '0':[],'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],})
                    
        #             for n in range(0,len(wvls)):
        #                 layer_data=layer_data.append({'wavelength':wvls[n],
        #                                               'extintion coeffient [km-1]':('{0:.7}'.format(np.array(AP_ml.ext).reshape((8,14,14,60,len(wvls)))[t,x,y,z,n])),
        #                                               'single scattering albedo':('{0:.7}'.format(np.array(AP_ml.ssa).reshape((8,14,14,60,len(wvls)))[t,x,y,z,n])),
        #                                               '0':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**0)),
        #                                               '1':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**1)),
        #                                               '2':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**2)),
        #                                               '3':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**3)),
        #                                               '4':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**4)),
        #                                               '5':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**5)),
        #                                               '6':('{0:.4f}'.format(AP_ml.g.values.reshape((8,14,14,60,len(wvls)))[t,x,y,z,n]**6))
        #                                               },ignore_index=True)
                        
        #             layer_file=atmo_path+'Aerosol/layers/'+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+'-z'+str(z)+'.LAYER'    
        #             layer_data.to_csv(layer_file, columns=['wavelength',
        #                                                    'extintion coeffient [km-1]',
        #                                                    'single scattering albedo',
        #                                                    '0',
        #                                                    '1',
        #                                                    '2',
        #                                                    '3',
        #                                                    '4',
        #                                                    '5',
        #                                                    '6'],sep=' ', encoding='utf-8', header=False,index=False)
                    
        #             aerosol_file=aerosol_file.append({'z(km)':('{0:.3f}'.format(c.z_mlvl.reshape((8,14,14,60))[t,x,y,z]/1000)),
        #                                               'aer_layer': layer_file
        #                                               },ignore_index=True)
                    
        #             aerosol_file.to_csv(atmo_path+'Aerosol/'+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y), 
        #                                 columns=['z(km)','aer_layer'],
        #                                 sep=' ', encoding='utf-8', header=False,index=False)
                    
                # """write the paths and input files for libradtran for each day and hs"""
                # input_file="/vols/satellite/home/bayer/libradtran/Libradtran-files/input_files/"+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+".txt"
                # output_file="/vols/satellite/home/bayer/libradtran/Libradtran-files/output_files/"+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y)+".txt"
                
                # f=open(input_file, "w+")
                # f.write
                # f.write("\ndata_files_path /home/nbayer/libRadtran-2.0.3/data/")
                # f.write("\natmosphere_file "+atmo_file)
                # f.write("\nsource solar  ../solar_flux/kurudz_0.1nm.dat per_nm ")   #line identifies the location of the extraterrestrial solar flux file which defines the spectral resolution.
                # f.write("\nmol_file O3 "+mol_file_o3+ ' mmr')
                # f.write("\nmol_file NO2 "+mol_file_no2+ ' mmr')
                # f.write("\nmol_file H2O "+mol_file_h2o+ ' mmr')
                # f.write("\npressure "+str(c.cams_sfc.psfc[t,x,y].values/100))
                # f.write("\nsza "+str(c.sza.reshape((8,14,14))[t,x,y]))
                # # f.write("\nphi0 "+str(c.azi.reshape((8,14,14))[t,x,y]))
                # f.write("\nmol_abs_param sbdart #spectralcalculation resolver, should be the best option for UV Indax calculations")
                # f.write("\naerosol_default")    #switch the use of aerosol data on
                # f.write("\naerosol_file explicit "+atmo_path+'Aerosol/'+times[t].strftime('%Y%m%d:%H')+'lat_lon:'+str(x)+'_'+str(y))
                # #f.write("\nck_lowtran_absorption O4 off")
                # f.write("\nrte_solver disort")
                # f.write("\ndisort_intcor moments")
                # # f.write("\nno_absorption mol")
                # # f.write("\nno_scattering mol")
                # f.write("\nwavelength "+str(min(wvls))+" "+ str(max(wvls)))    #Wavelength range [nm]
                # f.write("\noutput_process per_nm")
                # f.write("\nverbose") 
                # f.close()
            
                # Librad_path='/home/nbayer/libRadtran-2.0.3/bin/uvspec'
                # """Running LibRadTran with the input_file frim the step above and saving it in the output_file"""
                # """run the script from /home/nbayer/libRadtran-2.0.3/bin/ """
                # os.system(Librad_path+" < "+input_file+" > "+output_file) 


    if os.stat(output_file).st_size>100:
        """Reading the output_file and dividing the columns in variables"""
        data_out=np.loadtxt(output_file,dtype=float)
        """calculating and adding in the array the total irradiance values"""
        data_out=np.c_[data_out,data_out[:,1]+data_out[:,2],(data_out[:,4]+data_out[:,5]+data_out[:,6])]
        data_out=np.c_[data_out,np.ones(len(data_out)),np.ones(len(data_out))]
        
            
        for f in range(0,len(data_out)):
            if data_out[f,0]<=298:
                data_out[f,9]=1
                data_out[f,10]=data_out[f,9]*(data_out[f,1]+data_out[f,2])
            elif int(data_out[f,0]) in range(299,328):
                data_out[f,9]=10**(0.094*(298-int(data_out[f,0])))
                data_out[f,10]=data_out[f,9]*(data_out[f,1]+data_out[f,2])
            elif data_out[f,0] in range(328,401):
                data_out[f,9]=10**(0.015*(139-int(data_out[f,0])))
                data_out[f,10]=data_out[f,9]*(data_out[f,1]+data_out[f,2])
            else:
                data_out[f,9]=0
                data_out[f,10]=data_out[f,9]*(data_out[f,1]+data_out[f,2])
        
        UVI=UVI.append({'Date':times[t].strftime('%Y%m%d:%H'),'UVI': (np.sum(data_out[:,10])/25)},
                            ignore_index=True)
        # plt.step(data_out[:,0],data_out[:,7], label='Total downward irradiance '+times[t].strftime('%H'))
        plt.plot(data_out[:,0],data_out[:,7], label='Total downward irradiance '+times[t].strftime('%H'))
        # plt.step(data_out[:,0],data_out[:,1], label='Direct irradiance '+times[t].strftime('%H'))
        plt.plot(data_out[:,0],data_out[:,1], label='Direct irradiance '+times[t].strftime('%H'))
        plt.legend()
        plt.title(times[t].strftime('%Y%m%d:%H'))
        # plt.show()
        plt.savefig(times[t].strftime('%Y%m%d:%H')+ '_plot2.png', dpi=300)
        plt.close()
        # data[:,1].plot(title=str(output_file[47:60]),xlim=[280,500],ylim=[0,2000])
    else:
        UVI=UVI.append({'Date':times[t].strftime('%Y%m%d:%H'),'UVI': 'NaN'},
                            ignore_index=True)
        
print(UVI)
      
"""read netCDF"""
nc = nc4.Dataset('/vols/satellite/home/bayer/uv/netCDF/20190723.nc','r')
time=nc.variables['time'][:]
time=num2date(time[:],units='seconds since 1970-01-01T00:00:00')
spect=np.array(nc.variables['spect'][:]*1000)
uvind=np.array(nc.variables['uvind'][:])
 
nc=xr.open_dataset('/vols/satellite/home/bayer/uv/netCDF/20190723.nc')  
plt.plot(np.arange(290, 400.01, 0.1),nc.spect[6,:])

for p in range(0,len(time)):
    if time[p].strftime('%H%M')=='1200':
        break
    
fig=plt.figure()
# datetime.fromtimetamp(nc.variables['time'][195]).strftime("%B %d, %Y %I:%M:%S")
plt.plot(spect[p,0:400],'r',label='Melpitz '+time[p].strftime('%m-%d %H:%M'))
plt.plot(data_out[:,0],data_out[:,1],'b', label='Direct irradiance [Libradtran]')
plt.plot(data_out[:,0],data_out[:,1]+data_out[:,2],'g', label='Global irradiance [Libradtran]')
plt.plot([], [], ' ',label='UV-Index [Libradtran]= '+str(round(UVI['UVI'][len(UVI)-1],2)))
plt.plot([], [], ' ',label='UV-Index [Spectrometer]= '+str(round(uvind[p],2)))
# plt.xlim(501)
plt.legend()
plt.grid()
plt.title(time[p].strftime("%B %d, %Y %I:%M:%S"))
# plt.xticks(280,500)
plt.show()









