#! usr/bin/env python

import numpy as np
import constants as cst
import os
import warnings
import netCDF4 as nc


# For qpoint coordinates, from ASCII file
class QPTfile(object):

    def __init__(self, fname, nqpt, *args, **kwargs):

        self.nqpt = nqpt
        self.fname = fname
        self.qpt_array = self.read_qpoints(self.fname)

    def read_qpoints(self, fname=None):

        fname = fname if fname else self.fname

        arr = np.zeros((self.nqpt,3))
        # Read qpoint list text file and convert to array
        with open(fname, 'r') as f:
        
            lines = f.readlines()
            
            for iline, line in enumerate(lines):

                line = line.split('\n')[0]
                point = np.array([np.float(line.split()[0]), np.float(line.split()[1]), np.float(line.split()[2])])
                arr[iline,:] = point

        return arr

    def get_qnorm(self):

        arr = np.zeros((self.nqpt))

        for i, qpt in enumerate(self.qpt_array):

            arr[i] = np.linalg.norm(qpt)

        return arr
        

# Base class for netCDF files 
class CDFfile(object):

   def __init__(self, fname=None, read=True):

       self.fname = fname

       if read and self.fname:
           self.read_nc()

   def read_nc(self, fname=None):

       if not os.path.exists(fname):
           raise OSError('File not found : {}'.format(fname))

# For ZPR values of each qpoint for direct or indirect gap, extracted using zpr_analyzer.py
# Contribution computed with ElectronPhononCoupling, sending only Gamma and 1 qpoint at a time, with weights (0, 1)
# Contribution of each qpoint is NOT weighted accurately yet.
class ZPRfile(CDFfile):
    
    def __init__(self, *args, **kwargs):

        super(ZPRfile, self).__init__(*args, **kwargs)
        self.eig0 = None
        self.kpoints = None
        self.temp = None
        self.zp_ren = None
        self.ren_units = None
        self.td_ren = None
        self.smearing = None
        self.omega_se = None
        self.contr_qpt = None
        self.self_energy = None
        self.broadening = None
        self.spectral_function = None
        self.fan_g2 = None

    # Open the EPC file and read it
    def read_nc(self, fname=None):                                                                                                                                                  

        fname = fname if fname else self.fname
        super(ZPRfile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as ncdata:

#            self.eig0 = ncdata.variables['eigenvalues'][:,:,:]
#            self.kpoints = ncdata.variables['reduced_coordinates_of_kpoints'][:,:]
            self.temp = ncdata.variables['temperatures'][:]
#            self.td_ren = ncdata.variables['indirect_gap_renormalization_band'][:,:]
            self.zp_ren = ncdata.variables['indirect_gap_renormalization_band'][:,:] # ntemp,(val/cond)
            self.ren_units = ncdata.variables['indirect_gap_renormalization_band'].getncattr('units')
            self.zp_ren_modes = ncdata.variables['gap_renormalization_by_modes'][:,:] # nmodes, VBM/CBM
            self.indirect_zp_ren_modes = ncdata.variables['indirect_gap_renormalization_by_modes'][:,:] # nmodes, VBM/CBM

            self.zp_eigcorr_modes = ncdata.variables['reduced_eigenvalue_corrections_modes'][:,:,:,:]# spin, kpt, band, mode

            self.nqpt = len(ncdata.dimensions['number_of_qpoints'])
            #self.contr_qpt = ncdata.variables['qpt_contribution'][:,:,:,:,:]
            #self.contr_qpt_modes = ncdata.variables['qpt_contribution_modes'][:,:,:,:,:]
#            self.self_energy = ncdata.variables['self_energy'][:,:,:,:,:]
#            self.broadening = ncdata.variables['zero_point_broadening'][:,:,:]
#            self.spectral_function = ncdata.variables['spectral_function'][:,:,:,:]
#            self.smearing = ncdata.variables['smearing'][:]
#            self.omega_se = ncdata.variables['omegase'][:]

            #self.qpt_weight = ncdata.variables['qpoint_weight'][:]
            #self.fan_g2 = ncdata.variables['fan_g2'][:,:,:,:,:,:,:]
            #self.ddw_g2 = ncdata.variables['ddw_g2'][:,:,:,:,:,:,:]
            #self.deltaE_ddw = ncdata.variables['deltaE_ddw'][:,:,:,:,:,:]
            #self.fan_occterm = ncdata.variables['fan_occterm'][:,:,:,:,:,:,:,:]


    #        print(self.eig0)
    #        print(self.zp_ren)

    @property
    def nsppol(self):
        return self.eig0.shape[0] if self.eig0 is not None else None

    @property
    def nkpt(self):
        return self.eig0.shape[1] if self.eig0 is not None else None


    @property
    def max_band(self):
        return self.eig0.shape[2] if self.eig0 is not None else None

    @property
    def ntemp(self):
        return np.int(self.temp.shape[0]) if self.temp is not None else None

    @property
    def nfreq(self):
        return self.self_energy.shape[3] if self.self_energy is not None else None

    @property
    def nmodes(self):
        if self.zp_ren_modes is not None:
            return self.zp_ren_modes.shape[0] 
        elif self.fan_g2 is not None:
            return self.fan_g2.shape[-1]
        else:
            return None

    def read_files(): 
        
        f = zprdata
        if f.fname:
            f.read_nc()


class QptContribution(object):

    def __init__(self,

            qpt_fname = None,     
            zpr_flist = None,     
            nqpt = 1,   
            wtq = [1.0],
            units = 'eV', 
            indirect = False,
            check_sum = False,
            kpts_gap = [0,0],
            valence = None,

            **kwargs):

        if not qpt_fname:
            raise Exception('Must provide a file for qpt_fname')
        if not zpr_flist:
            raise Exception('Must provide a file for zpr_flist')
        if not valence:
            raise Exception('Must provide index for last occupied band')

        # Retreive passed arg values:
        self.qpt_fname = qpt_fname
        self.zpr_flist = zpr_flist
        self.nqpt = nqpt
        self.units = units
        self.indirect = indirect
        self.check_sum = check_sum
        self.valence = valence
        self.kpts_gap = kpts_gap

        if len(zpr_flist) != self.nqpt:
            raise Exception('Must provide nqpt files in zpr_flist')

        if len(wtq) != self.nqpt:
            raise Exception('Must provide nqpt weights in the wtq list')

        if self.indirect and self.kpts_gap[0] == self.kpts_gap[1]:
            raise Exception('Indirect gap must be between 2 different kpoints')

        # set and normalize qpoint weights
        self.set_weights(wtq)

        # Read the qpoints coordinates and get their norm
        qptfile = QPTfile(self.qpt_fname,self.nqpt)
        self.qpoints = qptfile.read_qpoints()
        self.qnorm = qptfile.get_qnorm()

        # For each qpoint file, retrieve the mode-split contribution to the VBM and CBM ZPR
        for i, fname in enumerate(self.zpr_flist):

            zprfile = ZPRfile(fname)
            zprfile.read_nc()

            if i == 0:

                self.nmodes = zprfile.nmodes
                self.contribution = np.zeros((self.nqpt,self.nmodes,2))

            self.contribution[i,:,0] = zprfile.zp_eigcorr_modes[0,self.kpts_gap[0],self.valence-1,:]
            self.contribution[i,:,1] = zprfile.zp_eigcorr_modes[0,self.kpts_gap[1],self.valence,:]
        self.convert_units()
        
        self.split_contribution()

        if check_sum:
            self.check_zpr()

    def convert_units(self):

        if self.units == 'eV':
            self.contribution = self.contribution*cst.ha_to_ev
        elif self.units == 'meV':
            self.contribution = self.contribution*cst.ha_to_ev*1E3
        else:
            raise Exception('Output units must be eV or meV')


    def set_weights(self, wtq, normalize=True):

        if normalize:
            self.wtq = np.array(wtq)/sum(wtq)
        else:
            self.wtq = np.array(wtq)

    def split_contribution(self):

        natom = int(self.nmodes/3.)
        print(natom)
        self.acoustic = np.einsum('qvn->qn',self.contribution[:,:3,:])
        self.TO = np.einsum('qvn->qn',self.contribution[:,3:(2*natom-1)+3,:])
        self.LO = np.einsum('qvn->qn',self.contribution[:,-(natom-1):,:])

    def check_zpr(self):

        # Sum the different weighted contribution, to validate that nothing is missing from the total ZPR
        zpr = np.einsum('q,qvn->n',self.wtq,self.contribution)
        print('Total ZPR is {} {}'.format(zpr[1]-zpr[0],self.units))
        print('    VBM : {} {}'.format(zpr[0],self.units))
        print('    CBM : {} {}'.format(zpr[1],self.units))

        print('acoustic : {}'.format(np.einsum('q,qn->n', self.wtq,self.acoustic)))
        print('TO : {}'.format(np.einsum('q,qn->n', self.wtq,self.TO)))
        print('LO : {}'.format(np.einsum('q,qn->n', self.wtq,self.LO)))


#    def write_netcdf(self):

#        create_directory(self.out_fname)

#        with nc.Dataset(self.nc_output, "w") as dts:



##############################################################

def create_directory(fname):

    dirname = os.path.dirname(fname)

    if not dirname:
        return
    if not os.path.exists(dirname):
        os.system('mkdir -p ' + dirname)


def compute(

        qpt_fname = None,
        zpr_flist = None,
        nqpt = 1,
        wtq = [1.0],
        units = 'eV',
        indirect = False,
        check_sum = False,
        valence = None,
        kpts_gap = None,

        **kwargs):

    contribution = QptContribution(

            qpt_fname = qpt_fname,
            zpr_flist = zpr_flist,
            nqpt = nqpt,
            wtq = wtq,
            units = units,
            indirect = indirect,
            check_sum = check_sum,
            valence = valence,
            kpts_gap = kpts_gap,

            **kwargs)

    #contribution.write_netcdf()

    return contribution

