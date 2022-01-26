#!/usr/bin/env python

from cdffile import CDFfile
import netCDF4 as nc

class EpcGkk2File(CDFfile):

    def read_nc(self, fname=None):

        # Open the OUT.nc file and read it
        fname = fname if fname else self.fname

        super(EpcGkk2File, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as root:

            self.kpoints = root.variables['reduced_coordinates_of_kpoints'][:,:]
            self.nmode = len(root.dimensions['number_of_modes'])
            self.nqpt = len(root.dimensions['number_of_qpoints'])
            self.nband = len(root.dimensions['max_number_of_states'])
            self.qpoints = root.variables['reduced_coordinates_of_qpoints'][:,:]

            self.gkk2_fan = root.variables['gkk2_fan_on_qgrid'][:,:,:,:,:] # (v,k,n,n',q)
            self.gkk2_ddw = root.variables['gkk2_ddw_on_qgrid'][:,:,:,:,:] # (v,k,n,n',q)


