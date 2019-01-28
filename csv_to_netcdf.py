#! usr/bin/env python

import netCDF4 as nc
import numpy as np

# This script converts PlotDigitizer outputfile to netcdf foramt, to be read by zpr_plotter
# Program info : http://plotdigitizer.sourceforge.net/

# PlotDigitizer .csv file name
fname = "AlP_gapT.csv"

#Axis units as converted
units = ['K','eV']

# NetCDF output filename
out_fname = "AlP_gapT.nc"

with open(fname,'r') as f:

    for i in range(4):
        f.readline()

    ndata = np.int(f.readline())

    line = f.readline().split(',')
    ax1 = line[0]
    ax2 = line[1].split('\n')[0]

    data1 = np.zeros((ndata))
    data2 = np.zeros((ndata))

    for i in range(ndata):
        line = f.readline().split(',')
        data1[i] = np.float(line[0])
        data2[i] = np.float(line[1])

f.close()

with nc.Dataset(out_fname,"w") as dts:

    dts.createDimension('number_of_points',ndata)
    
    data = dts.createVariable('ax1','d', ('number_of_points'))
    data[:] = data1
    data.long_name = ax1
    data.units = units[0]

    data = dts.createVariable('ax2','d', ('number_of_points'))
    data[:] = data2
    data.long_name = ax2
    data.units = units[1]




