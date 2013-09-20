#Copyright (C) 2013 Greg Miller <gmill002@gmail.com>
"""
Urban Phenology - Create x,y,t netcdf file, sample source data to create variables 

Usage: 
  python create_nc.py config.yaml

Inputs: 
  config.yaml: one or more yaml files specifying configuration parameters

Outputs:
  x,y,t netcdf data file with created variables 

Requires:
- numpy
- netCDF4
- arcpy 10.1, spatial analyst
- pyYAML
    
Procedure:

1. Initialize
  - create single config file listing all parameters
  - create x,y,t netcdf file using MODIS gdb spatial reference and time series dates 

2. MODIS EVI data
  - Extract MODIS gdb to netcdf file as 3-d variable, filter using QC data 

3. Census income data

4. Land cover    
  - Reclassify categories
  - aggregate into 250m grid using minimum majority percentage
  - add to netcdf as 2-d variable

5. Precipitation
  - aggregate (sum) daily weather into our time series dates
  - create derived lag or accumulate variables

    
Example Config File
===================

abbrev  : abq                                     # abbrev, file prefix
name    : Albuquerque                             # full name 
gdb     : I:\\urbphen\\ws_done\\datasets\\abq.gdb # source gdb containing EVI rasters ie abq_2004121 
outpath : I:\\urbphen\\netcdf\\abq                # created output folder
yaml    : I:\\urbphen\\netcdf\\abq\\info.yaml     # output of input config parameters
netcdf  : I:\\urbphen\\netcdf\\abq\\data.nc       # output netcdf data file

climate_type: semiarid   # not important
climate_rank: 3          # controls ordering on some plots
rainyearmonth : 12       # rain year cutoff month 
rainyearday : 31         # rain year cutoff day
maxt        : null       # max number of time steps
importmodis : True       # (True|False) import modis evi data 
airportwx   : True       # (True|False) import airport weather
prism_ppt   : False      # (True|False) import prism ppt

AnnualPrism:   # add prism annual ppt total on first period of month, 0 for others 
  var: APRI    # Output variable 
  path: 'I:\\urbphen\\AnnualPrism\\us_ppt' # input file

JoinCensus:                                             # extract data from census shapefiles 
  tractshp : I:\\urbphen\\JoinedCensusTracts\\abq.shp   # input census tract shapefile
  extractfields:                                        # sample these fields to netcdf variables 
    - FAVINC0
    - BLTBEF80
    - BLTAFT80
 
NLCD_Prep: # Generalize landcover from 30m to 250m using majority category if > THRESH otherwise NA
  file: I:\\urbphen\\NLCD2001\\save\\abq0  # input file 
  outfile : I:\\urbphen\\NLCD2001\\abq1    # created file 
  thresh : 0.70                            # aggregation threshold
  origvals: [11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95] # reclassify category mapping
  newvals : [10, 10, 21, 21, 23, 23, 30, 40, 40, 40, 50, 50, 81, 82, 90, 90]
  
landcover:                           # Create landcover variable 
  file: I:\\urbphen\\NLCD2001\\abq1  # input file
  name : NLCD                        # created variable name
  dtype: i1                          # created variable data type
  units: land cover category         # units description string 
  description: Land cover aggregated from NLCD 2001

scale       : 0.035                             # scale for plotting maps in R 
nlcdcol: I:\\urbphen\\NLCD2006\\nlcdinfo.csv    # csv file mapping land cover codes,labels,colors
excludevars:                                    # skip summarization for these variables 
  - xcoord
  - ycoord
  - tcoord
  - FAVINC0
  - BLTBEF80
  - BLTAFT80
  - NLCD

"""
import time
import datetime
import numpy as num
import os
import pdb
import re
import string
import yaml
import shutil
import csv
import subprocess
import sys

from netCDF4 import Dataset
import arcpy

arcpy.SetProduct("ArcInfo")
arcpy.CheckOutExtension("spatial")
arcpy.CheckOutExtension("3D")
arcpy.env.overwriteOutput = True

#netcdf fill value by data type string
NETCDF_FORMAT = 'NETCDF3_CLASSIC'
missvals = {'S1': '\x00', 
            'f4': 9.96920996839e+36,
            'f8': 9.96920996839e+36, 
            'i1': -127,
            'i2': -32767, 
            'i4': -2147483647,
            'i8': -9223372036854775806, 
            'u1': 255}

def RunYamlConfig(cfglist):
    """gather run parameters from a sequence of config dictionaries. """
    info = {}
    for cfgfile in cfglist:
        f = open(cfgfile, 'r')
        indict = yaml.load(f) 
        info.update(indict)
        f.close()
    
    print 'Run Parameters \n\n---\n%s---\n' % yaml.dump(info)
    RunConfig(info)
            
def RunConfig(info):
    """main processing function"""
    
    if 'ClipModisGDB' in info:
        ClipModisGDB(info, **info['ClipModisGDB'])
        return None
    
    if 'importmodis' in info:
        ImportModisGDB(info)
    
    if 'JoinCensus' in info:
        JoinCensus(info, **info['JoinCensus'])
        
    if 'NLCD_Prep' in info:
        NLCD_Prep(info, **info['NLCD_Prep'])
    
    if 'landcover' in info:
        ExtractCover(info, **info['landcover'])
        
    if info['airportwx']:
        AggregateDailyWx(info)
    
    if info['prism_ppt']:
        try:
            AddPrismWeather(info)
        except:
            print 'spatial precip failed'
        
    LagVar(info)
    AccumVar(info)
    
    if 'AnnualPrism' in info:
        AnnualPrism(info)
    
    if 'r' in info:
        for Rcmd in info['r']: 
            RunSubprocess(Rcmd,info)    
    
    WriteYaml(info)
    return None

def RunSubprocess(Rcmd, info):
    cmd = 'RScript "%s" "%s"' % (Rcmd, info['outpath'])
    print cmd
    proc = subprocess.Popen(cmd, shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while proc.poll() is None:
        time.sleep(10)
    
def ClipModisGDB(info, refraster, outpath, outgdb):
    """ clip input gdb to another extent and store as new gdb"""
    arcpy.env.workspace = info['gdb']
    rasters = arcpy.ListRasters('*')
    arcpy.env.snapraster = refraster
    print 'Clip and copy %d rasters' % len(rasters)
    newgdb = os.path.join(outpath, outgdb)
    if not arcpy.Exists(newgdb):
        print 'create Geodatabase %s' % newgdb
        arcpy.CreateFileGDB_management(outpath, outgdb)
    else:
        arcpy.env.workspace = newgdb
        newrasters = arcpy.ListRasters('*') 
        print '%s exists (%d rasters)' % (newgdb, len(newrasters))
        rasters = [r for r in rasters if r not in newrasters]
        print 'Clip and copy %d rasters' % len(rasters)
    
    arcpy.env.workspace = newgdb 
    for raster in rasters:
        arcpy.env.workspace = info['gdb']
        raster_clip = arcpy.sa.ExtractByMask(raster, refraster)
        arcpy.env.workspace = newgdb
        raster_clip.save(raster)
        print 'Clipped %s' % raster
    
def ImportModisGDB(info):
    """Import MODIS evi time series from gdb containing time stamped rasters """
    print 'Create Dataset "%s" at %s' % (info['abbrev'], info['outpath'])
    
    if os.path.exists(info['outpath']):
        shutil.rmtree(info['outpath'])
    
    os.mkdir(info['outpath'])
    
    #list the modis evi rasters in the gdb
    files, qafiles, times = ListEVIRasters(info['abbrev'], info['gdb'], info['maxt'])
    
    # get parameters required for creating the netcdf file 
    info = GetRefGrid(info, files[0], times)

    # create a netcdf file to hold the (x,y,t) data
    rootgrp = CreateNCDF(info['netcdf'], info['xcoords'], info['ycoords'], info['tcoords'])

    # add the modis evi data
    evi = rootgrp.createVariable('evi', 'i2', ('t','x','y'),fill_value=missvals['i2'])
    evi.description = 'MODIS evi'
    evi.units = '?'
    
    for (i, (raster,qa)) in enumerate(zip(files, qafiles)):
        qa_raw = arcpy.RasterToNumPyArray(qa)
        qa_int = num.array(UINT_To_Binary(qa_raw, start=0, end=2), dtype='i')
        qa_fail = qa_int >= 10
        data = arcpy.RasterToNumPyArray(raster)
        data[ qa_fail ] = missvals['i2']
        evi[i,:,:] = data
        print '\tgrid to netcdf %s ' % (raster)
    #
    rootgrp.close()
    WriteYaml(info)

def ListEVIRasters(abbrev, gdb, maxt):
    """list the modis evi rasters and dates in the gdb"""
    arcpy.env.workspace = gdb
    txt = string.join(arcpy.ListRasters('*'), sep=', ')
    print 'Import Modis from GDB %s' % arcpy.env.workspace

    files = re.findall('%s_evi_[0-9]{7}' % abbrev, txt)
    files.sort()
    
    qafiles = re.findall('%s_evi_qc_[0-9]{7}' % abbrev, txt)
    qafiles.sort()
    
    times = [int(re.findall('[0-9]{7}',e)[0]) for e in files]
    qatimes = [int(re.findall('[0-9]{7}',e)[0]) for e in qafiles]
    assert len(times) == len(qatimes), "Length of times and QA times do not match"
    
    if maxt != None:
        files = files[:maxt]
        qafiles = qafiles[:maxt]
        times = times[:maxt]
    
    return files, qafiles, times

def GetRefGrid(info, raster, times):
    """set the reference grid and times parameters"""
    arcpy.env.workspace = info['gdb'] 
    info['refraster'] = os.path.join(info['outpath'], 'ref')
    arcpy.CopyRaster_management(raster, info['refraster'])
    
    desc = arcpy.Describe(info['refraster'])
    info['spatialref'] = desc.spatialReference.name
    info['cellsize'] = int(desc.meanCellHeight)
    
    extent = desc.Extent
    info['xmin'] = int(extent.XMin)
    info['xmax'] = int(extent.XMax)
    info['ymin'] = int(extent.YMin)
    info['ymax'] = int(extent.YMax)
    
    info['xcoords'] = range(info['xmin'], info['xmax'], info['cellsize'])
    info['ycoords'] = range(info['ymin'], info['ymax'], info['cellsize'])
    info['tcoords'] = times
    info['refshape'] = [len(info['xcoords']),len(info['ycoords']),len(info['tcoords']),]

    info['dates'] = times
    #[datetime.datetime.strftime(datetime.datetime.strptime(str(t),'%Y%j'),'%Y-%m-%d') for t in times]
    
    return info

def CreateNCDF(fname, xcoords, ycoords, tcoords):
    """create a netcdf data file given x,y,t coordinates"""
    rootgrp = Dataset(fname, 'w', format=NETCDF_FORMAT)
    rootgrp.description = 'Urban Phenology Data'

    x = rootgrp.createDimension('x', len(xcoords))
    xcoord = rootgrp.createVariable('xcoord','i4',('x',))
    xcoord[:] = xcoords
    xcoord.units = 'meters east utm'

    y = rootgrp.createDimension('y', len(ycoords))
    ycoord = rootgrp.createVariable('ycoord','i4',('y',))
    ycoord[:] = ycoords
    ycoord.units = 'meters north utm'

    t = rootgrp.createDimension('t', None)
    times = rootgrp.createVariable('tcoord','i4',('t',))
    times[:] = tcoords
    times.units = '[YYYY][DayOfYear]'
    print rootgrp
    return rootgrp

def NLCD_Prep(xinfo, file='',outfile='', thresh=0, origvals=[], newvals=[]):
    """extract land cover raster"""
    info = ReadYaml(xinfo['yaml'])
    arcpy.env.workspace = info['outpath']
    #arcpy.Copy_management(file,os.path.join() )
    
    reclass = file +'_r'
    rmp = arcpy.sa.RemapValue([[o,n] for o, n in zip(origvals,newvals)])
    outReclass1 = arcpy.sa.Reclassify(file, "Value", rmp, 'NODATA')
    outReclass1.save(reclass)
    
    CalculateMajorityPercentage(reclass, outfile, info['refraster'], list(set(newvals)), thresh)

def CalculateMajorityPercentage(inraster, outraster, snapraster, values, thresh):
    """for each category, calculate the percentage that is majority type"""
    print "calculate majority percentage"
    arcpy.env.snapraster = snapraster
    outrc = ['%s%d' % (inraster, value) for value in values]
    for value, rcfile in zip(values,outrc):        
        isval = arcpy.sa.EqualTo(inraster, value)
        outBlockStat = arcpy.sa.BlockStatistics(isval, arcpy.sa.NbrRectangle(250, 250, "MAP"), 'MEAN', "DATA")
        outBlockStat.save(rcfile)
        print '\t' + rcfile
        
    arcpy.Mosaic_management(string.join(outrc, sep=';'), outrc[0], "MAXIMUM")
    keep = arcpy.sa.Reclassify(outrc[0], "Value", arcpy.sa.RemapRange([[0, thresh, 0],[thresh, 1, 1]]),'NODATA')
    refined = arcpy.sa.Times(keep, inraster)
    refined.save(outraster)
    
    [arcpy.Delete_management(rcfile) for rcfile in outrc]
        
def ExtractCover(xinfo, file='', name='', dtype='', units='', description=''):
    """extract land cover raster"""    
        # then aggregate as part of reprojection
    info = ReadYaml(xinfo['yaml'])
    arcpy.env.snapraster = info['refraster']
    prjrast = file+'p'
    clip = os.path.join(info['outpath'], name)
    
    arcpy.ProjectRaster_management(file, prjrast, info['refraster'], 'MAJORITY', cell_size=info['cellsize'])
    outExtractByMask = arcpy.sa.ExtractByMask(prjrast,  info['refraster'])
    outExtractByMask.save(clip)

    # create netcdf variable
    dev = arcpy.RasterToNumPyArray(clip)
    rootgrp = Dataset(info['netcdf'], 'a')
    lc = rootgrp.createVariable(name, dtype, ('t','x','y'),fill_value=missvals[dtype])
    lc.description = description
    lc.units = units
    for i,t in enumerate(info['tcoords']):
        lc[i,:,:] = dev
    
    rootgrp.close()
    WriteYaml(info)

def AggregateDailyWx(info):
    """aggregate daily weather station data
    [day, v1, v2, v3] to [startdate, enddate, v1]"""
        
    wxinfo = {
           'PRCP': {
                    'fun':'sum',
                    'units':'0.1 millimeters precipitation ',
                    'desc': '16 day precipitation total',
                    },
           'TMAX': {
                    'fun':'max',
                    'units': '0.1 degrees C (max)',
                    'desc': '16 day maximum temperature',
                    },
           'TMIN':{
                   'fun':'min',
                   'units':'0.1 degrees C (min)',
                   'desc':'16 day minimum temperature',
                   },
           }

    print 'aggregate wx'
    os.chdir("I:\\urbphen\\airport_wx")
    info = ReadYaml(info['yaml'])
    #coordinate the time bins
    bins = [datetime.datetime.strptime(str(t),'%Y%j') for t in info['tcoords'] ]
    bins.insert(0, bins[0] - (bins[1] - bins[0]) )
    bins = [int(datetime.datetime.strftime(t,'%Y%m%d')) for t in bins ]
    reftimes = bins[1:]
    
    #load daily wx data
    f = open('%s.csv' % info['abbrev'], 'r')
    reader = csv.reader(f)
    lines = [row for row in reader]
    f.close()
    wxvars = lines[0][1:] #names
    data = num.array(lines[1:])
    wxtimes = num.array(data[:,0],dtype='i') 
    wxdata = num.array(data[:,1:],dtype='i')

    # assign the bins
    idx = num.digitize(wxtimes,bins) - 1
    
    rootgrp = Dataset(info['netcdf'], 'a')
    
    for v, vname in enumerate(wxvars):

        #aggregate the data
        vinfo = wxinfo[vname]     
        if vinfo['fun'] == 'sum':
            data = [num.sum(wxdata[idx == i, v]) for i in num.arange(len(reftimes))]
        if vinfo['fun'] == 'min':
            data = [num.min(wxdata[idx == i, v]) for i in num.arange(len(reftimes))]
        if vinfo['fun'] == 'max':
            data = [num.max(wxdata[idx == i, v]) for i in num.arange(len(reftimes))]

        data = num.array(data,dtype='i')
        # write a summary file
        csvfile = os.path.abspath('%s_%s.csv' % (info['abbrev'], vname))
        bincount = num.bincount(idx+1)[1:-1]
        csvdata = num.transpose(num.vstack((bins[:-1],bins[1:], bincount, data)))
        num.savetxt(csvfile, csvdata, delimiter=',',fmt='%i')

        wx = rootgrp.createVariable(vname, 'i2', ('t','x','y'),fill_value=missvals['i2'])
        wx.description = vinfo['desc']
        wx.units = vinfo['units']
        for i , t in enumerate(info['tcoords']):
            
            wx[i,:,:] = data[i]
            
        del wx
        print '\t' + vname
    
    rootgrp.close()
    WriteYaml(info)
    
    return None

def AddPrismWeather(xinfo, varname='PRCPS', dtype='i4',normvar='PRCP'):
    print 'Add Prism Precip'
    info = ReadYaml(xinfo['yaml'])
    
    rootgrp = Dataset(info['netcdf'],'a')
    bv = rootgrp.variables[normvar]
    nv = rootgrp.createVariable(varname, bv.dtype, bv.dimensions, fill_value=missvals[dtype])
    nv.units = bv.units
    nv.description = bv.description
    
    abbrev = info['abbrev']
    arcpy.env.snapraster = info['refraster']
    arcpy.env.workspace = info['gdb']
    txt = string.join(arcpy.ListRasters('*'), sep=', ')
    files = re.findall('%s_ppt_[0-9]{6}' % abbrev, txt)
    rast = 0.0
    for raster in files:
        if arcpy.Exists(raster):
            rast = rast + arcpy.Raster(raster)
            print '\t' + raster
    
    rast.save('tempraster')    
    arcpy.Resample_management('tempraster', 'tempresample', info['cellsize'], "NEAREST")
    outExtractByMask = arcpy.sa.ExtractByMask('tempresample',  info['refraster'])
    outExtractByMask.save('tempclip')
    varsum = arcpy.RasterToNumPyArray('tempclip')
    varsum = varsum / num.mean(varsum)
    varsum = num.transpose(varsum)
    arcpy.Delete_management('tempraster')
    arcpy.Delete_management('tempresample')
    arcpy.Delete_management('tempclip')
    nv.description = nv.description + ' spatially normalized with PRISM average'

    for i in range(bv.shape[0]):
        nv[i,:,:] = num.array(bv[i,:,:] * varsum, dtype=dtype)
            
    rootgrp.close()
    WriteYaml(info)
        
def AnnualPrism(info):
    print 'Extract Prism PPT Annual Totals'
    assert arcpy.Exists(info['refraster']), 'reference raster not found'
    
    srcdir = info['AnnualPrism']['path']
    abbrev = info['abbrev']
    arcpy.env.workspace = srcdir
    arcpy.env.snapRaster = info['refraster']    
    arcpy.env.extent = info['refraster']
    arcpy.env.cellSize = info['refraster']
    rootgrp = Dataset(info['netcdf'],'a')
    nv = rootgrp.createVariable(info['AnnualPrism']['var'], 'i4', ('t','x','y',))
    nv[:] = 0
    years = [int(str(d)[:4]) for d in rootgrp.variables['tcoord'][:] ]
    uyears = num.sort(num.unique(years))
    for y in uyears:
        grid = os.path.join(srcdir, "us_ppt_{}14.asc".format(y))
        assert arcpy.Exists(grid), 'grid {} does not exist {}'.format(grid)
        
        egrid = os.path.join(srcdir, "u{}e".format(y))
        if arcpy.Exists(egrid): 
            arcpy.Delete_management(egrid)
        
        rgrid = os.path.join(srcdir, "u{}r".format(y))
        if arcpy.Exists(rgrid): 
            arcpy.Delete_management(rgrid)   
        
        pgrid = os.path.join(srcdir, "u{}p".format(y))
        if arcpy.Exists(pgrid): 
            arcpy.Delete_management(pgrid)
                
        arcpy.ProjectRaster_management(grid, pgrid, info['refraster'],'NEAREST')
        outExtract = arcpy.sa.ExtractByMask(pgrid, info['refraster'])
        outExtract.save(egrid)
        
        arcpy.Resample_management(egrid, rgrid, info['cellsize'], "NEAREST")
        
        print info['refshape']
        aa = arcpy.RasterToNumPyArray(rgrid)
        aa = num.transpose(aa)
        
        print aa.shape
        nv[years.index(y),:,:] = aa
        print y
            
    rootgrp.close()
    WriteYaml(info)
    
def UINT_To_Binary(inarray,start=0,end=16,width=16):
    """extract specific digits from the binary representation numbers in a n array
    
    inarray  : input array 
    width    : number of bits
    start;end: indices to extract 
    """
    return num.reshape([num.binary_repr(i,width=width)[start:end] for i in inarray.ravel()], inarray.shape)

def LagVar(info):
    """create a derived variable lagged by n steps in dim 3 
    
    mx[i,j,k] = mx0[i,j,k-lag]
    calculates only one timestep at a time which is slow but low memory"""
    info = ReadYaml(info['yaml'])
    
    vtype = 'L'
    
    if 'lag' not in info:
        info['lag'] = []

    rootgrp = Dataset(info['netcdf'], 'a')
    for (v,lagvar) in enumerate(info['lag']):
        bv = rootgrp.variables[lagvar['var']]
        info['lag'][v]['names'] = [lagvar['var']]
        info['lag'][v]['type'] = vtype
        
        for lag in lagvar['times']:
            varname = '%s%s%d' % (lagvar['var'], vtype, lag)
            nv = rootgrp.createVariable(varname, bv.dtype, bv.dimensions, fill_value=missvals['i2'])
            nv.units = bv.units
            nv.description = bv.description + ' with %d (x16) lag'          
            
            idx = num.array([i for i in range(bv.shape[0]) if (i-lag) >= 0])
            try: # as whole matrix lag
                nv[idx,:,:] = bv[(idx-lag),:,:]
            except: # do it one time slice at a time
                print 'one time step at a time'
                for i in idx:
                    nv[i,:,:] = bv[(i-lag),:,:]
            
            info['lag'][v]['names'].append(varname)
            
            print 'lag variable %s' % varname
    
    
        info['lag'][v]['times'].insert(0,0)
                        
    rootgrp.close()        
    WriteYaml(info)
    return None

def AccumVar(info):
    """create a derived variable with same dimensions, datatype 
    mx[i,j,k] = mx0[i,j,(k-accum):k]
    """
    vtype = 'A'
    
    info = ReadYaml(info['yaml'])
    if 'accumulate' not in info:
        info['accumulate'] = []
        
    rootgrp = Dataset(info['netcdf'], 'a')
    for (v,accumvar) in enumerate(info['accumulate']):
        bv = rootgrp.variables[accumvar['var']]
        info['accumulate'][v]['names'] = [accumvar['var']]
        info['accumulate'][v]['type'] = vtype
        
        for accum in accumvar['times']:
            varname ='%sA%d' % (accumvar['var'], accum)
            nv = rootgrp.createVariable(varname, bv.dtype, bv.dimensions,fill_value=missvals['i2'])
            nv.units = bv.units
            nv.description = bv.description + ' summed over '+str(accum)+'(x 16 day)'
            idx = [i for i in range(bv.shape[0]) if (i-accum) >= 0]
            for i in idx:
                nv[i,:,:,] = num.sum(bv[(i-accum):i,:,:], axis=0)
            
            info['accumulate'][v]['names'].append(varname)
            print 'Accumulate variable %s' % varname
    
        info['accumulate'][v]['times'].insert(0,0)
    
    rootgrp.close()        
    WriteYaml(info)

def JoinCensus(info, tractshp='', extractfields=[]):
    """given a polygon shapefile, extract a field as a raster then add to netcdf"""
    info = ReadYaml(info['yaml'])
    arcpy.env.snapraster = info['refraster']
    arcpy.env.workspace = os.path.dirname(tractshp)
    
    for field in extractfields:
        arcpy.PolygonToRaster_conversion(tractshp, field, field, "MAXIMUM_AREA","#", info['cellsize']) 
        
        outExtractByMask = arcpy.sa.ExtractByMask(field,  info['refraster'])
        outExtractByMask.save(field)
        data = num.array(arcpy.RasterToNumPyArray(field), dtype='i4')
    
        rootgrp = Dataset(info['netcdf'], 'a')
        lc = rootgrp.createVariable(field, 'i4', ('x','y'),fill_value=missvals['i4'])
        lc[:,:] = data
        print lc
        arcpy.Delete_management(field)
        
    rootgrp.close()
        
def ReadYaml(infile):
    """read a yaml config file"""
    f = open(infile,'r')
    info = yaml.load(f)
    f.close()
    return info

def WriteYaml(info):
    """write a yaml config file"""
    f = open(info['yaml'], 'w')
    yaml.dump(info, f, default_flow_style=False)
    f.close()
    

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        RunYamlConfig(sys.argv[1:])       
