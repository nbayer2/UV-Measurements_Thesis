#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 12:59:26 2019

@author: bayer

__author__ = "Nicolas Bayer"
__maintainer__ = "Nicolas Bayer"
__email__ = "bayer@tropos.de"
__status__ = "Production"
"""
import os
import read_bts2048rh as bts
import BTS2plot
#from Submodule import PlotHeatmap as plthtmp
import BTS2NetCDF 
import argparse
import json
import configparser
import platform
import pandas as pd
import xarray as xr
import logging



# get name of directory where main script is located
current_dirname = os.path.dirname(os.path.realpath(__file__))

# define log_path_file + create dir
log_path_file = current_dirname + "/log/uv_processing.log"
if not os.path.isdir(  os.path.dirname( log_path_file ) ):
    os.makedirs( os.path.dirname( log_path_file ) )
    print('Create directory     : ' + os.path.dirname(log_path_file) )


"""Insert de initial and final dates as strings as 20190107(year:2019/month:01/day:07)"""

"""for calling the function from the terminal"""
parser = argparse.ArgumentParser(description='Process UV radiometer measurements.') 
parser.add_argument('-s', type=str, dest='id', # la variable se guarda en args.id como string
                    help='Insert the initial date as 20190107(y:2019 m:01 d:07)')
parser.add_argument('-e', type=str, dest='fd',
                    help='Insert the final date as 20190107(y:2019 m:01 d:07)')
parser.add_argument('-i', '--image', action='store_true', 
                    help="create images files")
parser.add_argument('-n', '--netcdf', action='store_true', 
                    help="create netCDF files")
parser.add_argument('-st', '--statistics', action='store_true', 
                    help="create statistics of missing files")
parser.add_argument('-l', '--loglevel', default='INFO',
                    help="define loglevel of screen INFO (default) | WARNING | ERROR ")
args = parser.parse_args()


# create logger with 'UV'
logger = logging.getLogger('uv-processing')
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
fh = logging.FileHandler( log_path_file )
fh.setLevel(logging.DEBUG)

# create/check level
screen_level = logging.getLevelName(args.loglevel)

# create console handler with a higher log level
ch = logging.StreamHandler()
#ch.setLevel(logging.WARNING)
ch.setLevel(screen_level)

# create formatter and add it to the handlers
formatter = logging.Formatter(fmt='%(asctime)s | %(name)-24s | %(levelname)-8s | %(message)s | %(module)s (%(lineno)d)', datefmt='%Y-%m-%d %H:%M:%S',)

# add formatter to the handlers
fh.setFormatter(formatter)
ch.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)


# add first log mesage
logger.info('Start uv-processing')



"""Break in case the dates weren't correct"""
if len(args.id)!=8 or len(args.fd)!=8:
    logger.error('Wrong date introduced, 8 digits required!')
    exit()
if int(args.id)>int(args.fd):
    logger.error('Wrong date, start date < = end date required!')
    exit()
    
"""Check python version"""
python_version = platform.python_version().split(".")
if int(python_version[0]) < 3:
    logger.error( "Your python version is: " + platform.python_version() )
    logger.error( "Script will be terminated cause python version < 3 is required !" )
    exit()
  
  
"""Check python Submodule is already installed"""
if (args.statistics):
    try:
        from Submodule import PlotHeatmap as plthtmp
    except ImportError:
        logger.error('There was no such module installed: PlotHeatmap')
        exit()
    
    

def statistic(i8date,f8date):
    methodbts = "global" 
    
    config_file     = current_dirname + '/config/config.ini'
    
    
   
    """Read config file"""
    config = configparser.ConfigParser()
    
    
    
    """Read private config file"""
    if not os.path.isfile( config_file ):
        logger.error('File config/config.ini not exists, use DEFAULT config/templates/config.ini instead!' + config_file)
        exit()
    else:
        config.read( config_file )

    """Check if directories etc exists"""
    if not os.path.isdir( config.get('DEFAULT','main_path') ):
        logger.error('Path main_path not exists '+ config.get('DEFAULT','main_path'))
        quit()
        
    if not os.path.isdir( config.get('DEFAULT','image_path') ):
        logger.error('Path image_path not exists '+ config.get('DEFAULT','image_path'))
        quit()
        
    if not os.path.isdir( config.get('DEFAULT','netCDF_path') ):
        logger.error('Path netcdf_path not exists '+ config.get('DEFAULT','netCDF_path'))
        quit()
        
    json_file       = current_dirname + "/" + config.get('DEFAULT','json_file')
    if not os.path.isfile( json_file ):
        logger.error( 'File json not exists '+ json_file )
        quit()


        
    """Load content of json_file to python variable cfjson"""
    cfjson={}
    with open( json_file ) as f:
            cfjson= json.load(f)
    
    
    """ a lookup value dict for missing files"""
    dict_lookup_missing_value = {
        "file_not_exists": 1,
        "file_empty": 0.8,
        "file_less_than_1mb": 0.3,
        "file_size_ok" : 0
    }
    
    
    """add dataframe to plot missing files"""
    missing_files_key_name = "missing data"
    df = pd.DataFrame({'date' : [], missing_files_key_name : [] })
    
    """pandas time counter vector instead of loop"""    
    dates=pd.date_range(args.id, args.fd,freq='1D', name=str, normalize=False) 

    # pd.to_numeric(dates.str.replace('-',''))
    for date in dates:
        
        """Compose PathFileName of or0-File"""
        or0_file = config.get('STATION','station_prefix') + \
                        date.strftime('%y').zfill(2)+date.strftime('%m')+ date.strftime('%d')+".OR0"
        path_file = ""
        
        
        """path_file is dependent on the main_path_tree"""
        if (config.get('DEFAULT','main_path_tree') == 'flat') :
            path_file = config.get('DEFAULT','main_path') + or0_file
        elif (config.get('DEFAULT','main_path_tree') == 'yyyy/mm/dd/') :
            path_file = config.get('DEFAULT','main_path') + date.strftime('%Y/%m/%d/') + or0_file
            
        
        """Compose PathFileName of image-File and add to config to use in plot module"""
        image_path_file = config.get('DEFAULT','image_path') + eval( config.get('DEFAULT','image_subpath_file_regex') )
        config.set("DEFAULT", "image_path_file", image_path_file)
        
        
        """Compose PathFileName of netcdf-File and add to config"""        
        netcdf_path_file = config.get('DEFAULT','netcdf_path') + eval( config.get('DEFAULT','netcdf_subpath_file_regex') )
        config.set("DEFAULT", "netcdf_path_file", netcdf_path_file)


        """check file exists"""
        if os.path.isfile(path_file):  # see if the .OR0 file exist
            if os.stat(path_file).st_size<1:  # controls if the file is not empty
                logger.warn('file is empty '+path_file)
                if args.statistics:
                    df = df.append({'date': date.date(), missing_files_key_name : dict_lookup_missing_value["file_empty"]}, ignore_index=True)
            else:
                if args.statistics:
                    if os.stat(path_file).st_size<1048576:  # controls if the file is less than 1mb
                        df = df.append({'date': date.date(), missing_files_key_name : dict_lookup_missing_value["file_less_than_1mb"]}, ignore_index=True)
                    else:
                        df = df.append({'date': date.date(), missing_files_key_name : dict_lookup_missing_value["file_size_ok"]}, ignore_index=True)
                        
                """Check if netcdf (sub)directory exists"""
                if args.netcdf or args.image:
                    
                    """checking if the directory already exists, create subdir"""
                    if not os.path.isdir( os.path.dirname(netcdf_path_file) ):
                        os.makedirs( os.path.dirname(netcdf_path_file) )
                        logger.info('Create directory     : ' + os.path.dirname(netcdf_path_file) )
                
                """Obtanin the directory data from the OR0 files"""
                if args.netcdf:
                    d_bts1day=bts.read_oro_bts(path_file,methodbts,date.strftime('%Y%m%d'))
                    if args.netcdf:
                        
                        """checking if the file does already exist, and delete if so."""
                        if os.path.isfile( netcdf_path_file ):
                            os.remove( netcdf_path_file )
                            BTS2NetCDF.netCDF_file(d_bts1day,netcdf_path_file,cfjson) 
                        else:
                            """Save the data processed by the bts function in a netCDF file"""
                            BTS2NetCDF.netCDF_file(d_bts1day,netcdf_path_file,cfjson)
                
                if args.image:
                    """checking if the file does already exist"""
                    """ If it exist, opens data in xarray.
                        If not, creates netCDF file and load data in xarray"""
                    if os.path.isfile( netcdf_path_file ):
                        nc=xr.open_dataset( netcdf_path_file )
                    else:
                        d_bts1day=bts.read_oro_bts(path_file,methodbts,date.strftime('%Y%m%d')) 
                        BTS2NetCDF.netCDF_file(d_bts1day,netcdf_path_file,cfjson)
                        nc=xr.open_dataset(netcdf_path_file)
                    
                    """Plotting data"""
                    BTS2plot.plotme(nc,date,config)
                    

        else:
            logger.error("File not exist "+ path_file)
            if args.statistics:
                df = df.append({'date': date.date(), missing_files_key_name : dict_lookup_missing_value["file_not_exists"]}, ignore_index=True)
        
        """plot statistics"""
        if (args.statistics):
            """generate filename"""
            picture_filename = \
                config.get('DEFAULT','image_path') + config.get('STATION','station_prefix') + '_' + \
                str( df['date'][0].strftime('%Y%m%d') ) + \
                '_' + \
                str(df['date'][len(df.index)-1].strftime('%Y%m%d') ) + \
                '_missing_data'
                
            """plot first or second half of year"""
            if (date.strftime('%m%d')=='0630' or date.strftime('%-m%d')=='1231'):
                logger.info('Plot ' + picture_filename )
                
                """transform column date to datetime"""
                df['date'] =  pd.to_datetime(df['date'])
                
                plthtmp.main(
                    {
                    'data_import_type' : 'DataFrame',
                    'picture_filename' : picture_filename,
                    'DataFrame' : df
                    }
                )
                
                """init new dataframe"""
                df = pd.DataFrame({'date' : [], missing_files_key_name : [] })
                
            #plot period after the last half year
            elif (date.strftime('%Y%m%d') == f8date):
                logger.info('Plot short period: ' + picture_filename )
                    
                # transform column date to datetime
                df['date'] =  pd.to_datetime(df['date'])
                
                plthtmp.main(
                    {
                    'data_import_type' : 'DataFrame',
                    'picture_filename' : picture_filename,
                    'DataFrame' : df
                    }
                )
    logger.info('End uv-processing')

#####################################################################################                                                            
statistic(args.id,args.fd)






