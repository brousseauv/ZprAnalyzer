#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import sys
import netCDF4 as nc
import sys
import os
import warnings
import itertools as itt
from scipy.interpolate import interp1d

############################
# This code retreives zpr and temperature-dependent corrections to the electronic eigenvalues computed with ABINIT and post-processed with G.Antonius's ElectronPhononCoupling module

# For a given subset of kpts, it averages the correction on all symmetry-equivalent kpts and reduces the grid subset to a 2D array foollowing a defined path.

# The output file gives the corrected electronic eigenvalues at different temperratures on the symmetry-reduced path.

# V.Brousseau 2017/07

############################

# Variables and classes
class VariableContainer:pass

cst = VariableContainer()
cst.ha_to_ev = np.float(27.211396132)
cst.ev_to_ha = np.float(1./cst.ha_to_ev)
cst.tolkpts = np.float(0.00001)

# Base class for netCDF files 
class CDFfile(object):
    
    def __init__(self, fname=None, read=True):

        self.fname = fname

        if read and self.fname:
            self.read_nc()

    def read_nc(self, fname=None):
            
        if not os.path.exists(fname):
            raise OSError('File not found : {}'.format(fname))            

class OUTfile(CDFfile):

    def __init__(self, *args, **kwargs):

        super(OUTfile, self).__init__(*args, **kwargs)
        self.symrel = None

    def read_nc(self, fname=None):

        fname = fname if fname else self.fname
        super(OUTfile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as ncdata:

            self.symrel = ncdata.variables['symrel'][:]

    @property
    def nsym(self):
        return int(self.symrel.shape[0]/9) if self.symrel is not None else None


class EIGfile(CDFfile):

    def __init__(self, *args, **kwargs):

        super(EIGfile, self).__init__(*args, **kwargs)
        self.eig0 = None
        self.kpoints = None

    # Open the EIG.nc file and read it
    def read_nc(self, fname=None):

        fname = fname if fname else self.fname
        super(EIGfile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as ncdata:
            
            self.eig0 = ncdata.variables['Eigenvalues'][:,:,:]
            self.kpoints = ncdata.variables['Kptns'][:,:]

    @property
    def nsppol(self):
        return self.eig0.shape[0] if self.eig0 is not None else None

    @property
    def nkpt(self):
        return self.eig0.shape[1] if self.eig0 is not None else None
    
    @property
    def max_band(self):
        return self.eig0.shape[2] if self.eig0 is not None else None

    
class EPCfile(CDFfile):
    
    def __init__(self, *args, **kwargs):

        super(EPCfile, self).__init__(*args, **kwargs)
        self.eig0 = None
        self.kpoints = None
        self.natom = None
        self.temp = None
        self.zp_ren = None
        self.td_ren = None
        self.smearing = None
        self.omega_se = None
        self.contr_qpt = None
        self.self_energy = None
        self.broadening = None
        self.spectral_function = None
        self.fan_g2 = None
        self.zpr_ren_modes = None
        self.qpt_weight = None
        self.fan_occ = None
        self.fan_unocc = None
        self.dw_occ = None
        self.dw_unocc = None

    # Open the EPC file and read it
    def read_nc(self, fname=None):

        fname = fname if fname else self.fname
        super(EPCfile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as ncdata:

            self.natom = len(ncdata.dimensions['number_of_atoms'])
            self.eig0 = ncdata.variables['eigenvalues'][:,:,:]
            self.kpoints = ncdata.variables['reduced_coordinates_of_kpoints'][:,:]
            self.temp = ncdata.variables['temperatures'][:]
            self.td_ren = ncdata.variables['temperature_dependent_renormalization'][:,:,:,:]
        # ONLY IF RAN WITH MY DEVELOP BRANCH
        #    self.fan_occ = ncdata.variables['temperature_dependent_renormalization_fan_occ'][:,:,:,:]
        #    self.fan_unocc = ncdata.variables['temperature_dependent_renormalization_fan_unocc'][:,:,:,:]
        #    self.ddw_occ = ncdata.variables['temperature_dependent_renormalization_ddw_occ'][:,:,:,:]
        #    self.ddw_unocc = ncdata.variables['temperature_dependent_renormalization_ddw_unocc'][:,:,:,:]

            self.zp_ren = ncdata.variables['zero_point_renormalization'][:,:,:]
            self.zp_ren_modes = ncdata.variables['zero_point_renormalization_by_modes'][:,:,:,:]
            self.nqpt = len(ncdata.dimensions['number_of_qpoints'])
            #self.contr_qpt = ncdata.variables['qpt_contribution'][:,:,:,:,:]
            #self.contr_qpt_modes = ncdata.variables['qpt_contribution_modes'][:,:,:,:,:]
            self.self_energy = ncdata.variables['self_energy'][:,:,:,:,:]
            self.broadening = ncdata.variables['zero_point_broadening'][:,:,:]
            self.spectral_function = ncdata.variables['spectral_function'][:,:,:,:]
            self.smearing = ncdata.variables['smearing'][:]
            self.omega_se = ncdata.variables['omegase'][:]

#            self.qpt_weight = ncdata.variables['qpoint_weight'][:]
#            self.fan_g2 = ncdata.variables['fan_g2'][:,:,:,:,:,:,:]
#            self.ddw_g2 = ncdata.variables['ddw_g2'][:,:,:,:,:,:,:]
#            self.deltaE_ddw = ncdata.variables['deltaE_ddw'][:,:,:,:,:,:]
#            self.fan_occterm = ncdata.variables['fan_occterm'][:,:,:,:,:,:,:,:]
#            self.ddw_tdep = ncdata.variables['ddw_tdep'][:,:,:]

               
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
            return self.fan_g2.shape[4]
        else:
            return None




class EigenvalueCorrections(object):

    #Input files
    epc_fname = None
    out_fname = None

    #Parameters
    nsppol = None
    nkpt = None
    valence = None
    max_band = None
    eig0 = None
    eigcorr = None
    td_ren = None
    temp = None
    ntemp = None
    grid = None
    nsym = None
    nfreq = None
    symmat = None
    smearing = None
    omega_se = None
    symrel = None
    units = 'eV'
    kpoints = None
    nqpt = None
    contr_qpt = None
    self_energy = None
    broadening = None
    spectral_function = None
    qpt_weight = None

    reduced_nkpt = None
    reduced_kpts_index = None
    reduced_kpath = None
    reduced_eig0 = None
    reduced_eigcorr = None
    reduced_td_ren = None
    reduced_contr_qpt = None
    gap_renorm = None
    gap_loc = None
    gap_energy = None
    split_kpt = None
    nband_contr = None

    #Options
    bands_to_print = None
    points_to_print = None
    reduce_path = False
    contribution = False
    full_contribution = False
    split_occupied_subspace = False
    gap_ren = False
    temperature = False
    senergy = False
    broad = False
    spectral = False
    indirect = False
    split = False,
    modes = False
    split_contribution = False
    spline = False
    spline_factor = 5

    def __init__(self,

                #Options
                bands_to_print = None,
                points_to_print = None,
                reduce_path = False,
                contribution = False,
                full_contribution = False,
                gap_ren = False,
                temperature = False,
                senergy = False,
                broad = False,
                spectral = False,
                indirect = False,
                split = False,
                modes = False,
                split_contribution = False,
                split_occupied_subspace = False,
                spline = False,
                spline_factor = 5,

                #Parameters
                nsppol = None,
                nkpt = None,
                valence = None,
                max_band = None,
                eig0 = None,
                eigcorr = None,
                td_ren = None,
                temp = None,
                ntemp = None,
                grid = None,
                nsym = None,
                nfreq = None,
                symrel = None,
                symmat = None,
                smearing = None,
                omega_se = None,
                qpt_weight = None,
                nband_contr = None,

                units = 'eV',
                kpoints = None,
                nqpt = None,
                contr_qpt = None,
                self_energy = None,
                broadening = None,
                spectral_function = None,

                reduced_nkpt = None,
                reduced_kpts_index = None,
                reduced_kpath = None,
                reduced_eig0 = None,
                reduced_eigcorr = None,
                reduced_td_ren = None,
                rootname = 'eigcorr.out',
                reduced_contr_qpt = None,
                gap_loc = None,
                gap_renorm = None,
                gap_energy = None,
                split_kpt = None,

                #Input files
                epc_fname = None,
                out_fname = None,

                **kwargs):

        # Check if filenames are present
        if not epc_fname:
            raise Exception('Must provide a file for epc_fname')
        if not out_fname:
            raise Exception('Must provide a file for out_fname')

        # Set input files 
        self.epc_fname = epc_fname
        self.out_fname = out_fname

        # Define files classes
        self.epc = EPCfile(epc_fname, read=False)
        self.out = OUTfile(out_fname, read=False)

        # Define options
        self.reduce_path = reduce_path
        self.bands_to_print = bands_to_print
        self.points_to_print = points_to_print
        self.contribution = contribution
        self.full_contribution = full_contribution
        self.split_contribution = split_contribution
        self.split_occupied_subspace = split_occupied_subspace
        self.gap_ren = gap_ren
        self.units = units
        self.temperature = temperature
        self.valence = valence
        self.grid = grid
        self.indirect = indirect
        self.split_kpt = split_kpt
        self.split = split
        self.modes = modes
        self.nband_contr = nband_contr
        self.spline = spline
        self.spline_factor = spline_factor

        # Units warning
        if self.units is not 'eV' and self.units is not 'meV':
            raise Exception('Units must be eV or meV')
        if not self.valence:
            raise Exception('Must provide the last valence band number as valence')

        # Read input files
        self.read_files()

        # Set parameters
        self.nsppol = self.epc.nsppol
        self.eig0 = self.epc.eig0
        self.kpoints = self.epc.kpoints
        self.natom = self.epc.natom
        self.nqpt = self.epc.nqpt
        self.temp = self.epc.temp
        self.nkpt = self.epc.nkpt
        self.ntemp = self.epc.ntemp
        self.max_band = self.epc.max_band
        self.nfreq = self.epc.nfreq
        self.smearing = self.epc.smearing
        self.omega_se = self.epc.omega_se

        self.set_rootname(rootname)

        self.nsym = self.out.nsym
        self.symrel = self.out.symrel

        if self.temperature:
            if self.modes:
                raise Exception('Mode decomposition only implemented at T=0 theory. Put Temperature to False.')
            else:
                self.td_ren = self.epc.td_ren
        else:
            if self.modes:
                self.td_ren = self.epc.zp_ren_modes
                self.td_ren = np.einsum('vskn->skn',self.td_ren)
            else:
                self.td_ren = self.epc.zp_ren
            self.td_ren = np.expand_dims(self.td_ren, axis=3) 

            # Change ntemp to 1
            self.ntemp = 1
            self.temp = np.array([0.])

        self.senergy = senergy
        self.spectral = spectral
        self.broad = broad
        self.self_energy = self.epc.self_energy
        self.broadening = self.epc.broadening
        self.spectral_function = self.epc.spectral_function

        # Reshape symmetry matrices
        self.symmat = self.reshape_symmat(self.symrel)

        # Compute eingenvalue corrections for each temperature
        self.eigcorr = self.get_corrections(self.nkpt, self.eig0, self.td_ren)

        # Get contribution from each qpt and/or mode
        if self.contribution :
            if self.modes:
                self.contr_qpt_modes = self.epc.contr_qpt_modes
                self.contr_qpt_modes = np.einsum('ijklm -> iklmj',self.contr_qpt_modes) #(nsppol,nkpt,nband,nqpt,nmodes)
                self.zp_ren_modes = self.epc.zp_ren_modes
                self.zp_ren_modes = np.einsum('ijkl->jkli',self.zp_ren_modes) #(nsppol,nkpt,nband,nmodes)
                self.nmodes = self.epc.nmodes

            else:
                self.contr_qpt = self.epc.contr_qpt
        else:
            if self.modes:
                self.zp_ren_modes = self.epc.zp_ren_modes
                self.zp_ren_modes = np.einsum('ijkl->jkli',self.zp_ren_modes) #(nsppol,nkpt,nband,nmodes)
                self.nmodes = self.epc.nmodes

        if self.full_contribution or self.split_contribution: 
            self.fan_g2 = self.epc.fan_g2
            self.ddw_g2 = self.epc.ddw_g2
            self.deltaE_ddw = self.epc.deltaE_ddw
            self.fan_occterm = self.epc.fan_occterm
            self.qpt_weight = self.epc.qpt_weight
            self.ddw_tdep = self.epc.ddw_tdep

        if self.split_occupied_subspace:
            self.fan_occ = self.epc.fan_occ
            self.fan_unocc = self.epc.fan_unocc
            self.ddw_occ = self.epc.ddw_occ
            self.ddw_unocc = self.epc.ddw_unocc

        # Average on symmetry-equivalent kpoints
        if self.reduce_path :
            self.reduced_nkpt, self.reduced_kpts_index, self.reduced_kpath = self.reduce_symmetry_equivalent(self.kpoints)
            self.reduced_eig0 = self.average_kpts(self.reduced_kpts_index, self.eig0, 1)

            self.reduced_td_ren = self.average_kpts(self.reduced_kpts_index, self.td_ren, self.ntemp)

            if self.temperature is False:
                # Adjust array dimensions
                self.reduced_td_ren = np.expand_dims(self.reduced_td_ren, axis=3)

            self.reduced_eigcorr = self.get_corrections(self.reduced_nkpt, self.reduced_eig0, self.reduced_td_ren)
            
            if self.contribution:
                if self.modes:
                    self.reduced_zp_ren_modes = np.zeros((self.nsppol, self.reduced_nkpt, self.max_band, self.nmodes))
                    self.reduced_zp_ren_modes = self.average_kpts(self.reduced_kpts_index, self.zp_ren_modes, self.nmodes)
                    self.reduced_contr_qpt_modes = np.zeros((self.nsppol, self.reduced_nkpt, self.max_band, self.nqpt, self.nmodes))
                    for v in range(self.nmodes):
                        self.reduced_contr_qpt_modes[:,:,:,:,v] = self.average_kpts(self.reduced_kpts_index, self.contr_qpt_modes[:,:,:,:,v], self.nqpt) 
                else:
                    self.reduced_contr_qpt = np.zeros((self.nsppol,self.reduced_nkpt,self.max_band,self.nqpt,self.ntemp))
                    for t in range(self.ntemp):
                        self.reduced_contr_qpt[:,:,:,:,t] = self.average_kpts(self.reduced_kpts_index, self.contr_qpt[:,:,:,:,t], self.nqpt)
            else:
                if modes:
                    self.reduced_zp_ren_modes = np.zeros((self.nsppol, self.reduced_nkpt, self.max_band, self.nmodes))
                    self.reduced_zp_ren_modes = self.average_kpts(self.reduced_kpts_index, self.zp_ren_modes, self.nmodes)


            if self.full_contribution or self.split_contribution:
                if not self.nband_contr:
                    raise Exception('Must provide nband_contr')

                self.reduced_fan_g2 = np.zeros((self.nsppol,self.reduced_nkpt, self.nband_contr, self.nband_contr, self.nmodes,self.nqpt,2))
                for n,v,ri in itt.product(range(self.nband_contr), range(self.nmodes),range(2)):
                    self.reduced_fan_g2[:,:,:,n,v,:,ri] = self.average_kpts_contr(self.reduced_kpts_index, self.fan_g2[:,:,:,n,v,:,ri], self.nqpt)

                self.reduced_ddw_g2 = np.zeros((self.nsppol,self.reduced_nkpt, self.nband_contr, self.nband_contr, self.nmodes,self.nqpt,2))
                for n,v,ri in itt.product(range(self.nband_contr), range(self.nmodes),range(2)):
                    self.reduced_ddw_g2[:,:,:,n,v,:,ri] = self.average_kpts_contr(self.reduced_kpts_index, self.ddw_g2[:,:,:,n,v,:,ri], self.nqpt)

                self.reduced_deltaE_ddw = np.zeros((self.nsppol, self.reduced_nkpt, self.nband_contr, self.nband_contr, self.nqpt,2))
                for n,ri in itt.product(range(self.nband_contr),range(2)):
                    self.reduced_deltaE_ddw[:,:,:,n,:,ri] = self.average_kpts_contr(self.reduced_kpts_index,self.deltaE_ddw[:,:,:,n,:,ri],self.nqpt)

                self.reduced_fan_occterm = np.zeros((self.nsppol, self.reduced_nkpt, self.nband_contr, self.nband_contr, self.nmodes, self.ntemp, self.nqpt,2))
                for n,v,t,ri in itt.product(range(self.nband_contr), range(self.nmodes), range(self.ntemp), range(2)):
                    self.reduced_fan_occterm[:,:,:,n,v,t,:,ri] = self.average_kpts_contr(self.reduced_kpts_index,self.fan_occterm[:,:,:,n,v,t,:,ri], self.nqpt)


#                print(self.kpoints[60,:])
#                print(self.reduced_kpath[10,:])
#                diff = np.abs(self.contr_qpt[0,60,27,:,:]-self.reduced_contr_qpt[0,10,27,:,:])
#                for iqpt,t in itt.product(range(self.nqpt),range(self.ntemp)):
#                    if diff[iqpt,t] > 1E-6:
#                        print(iqpt,t,diff)
                
                
#                for aa in range(2):
#                    if aa==0:
#                        print('VB')
#                    else:
#                        print('CB')
#                    for q in range(self.nqpt):
#
#                        print(q+1)
#                        #print('A:',self.contr_qpt[0,60,27+aa,q,-1])
#                        #print(self.reduced_kpath[9,:])
#                        #print('A-1 :',self.reduced_contr_qpt[0,9,27+aa,q,-1])
#                        #print('A-2 :',self.reduced_contr_qpt[0,8,27+aa,q,-1])
#                        #print('A-3 :',self.reduced_contr_qpt[0,7,27+aa,q,-1])
#                        print('around gap :',self.reduced_contr_qpt[0,5,27+aa,q,-1])
#
#                        #for dd in range(6):
#                            #print(q+1,self.kpoints[54+dd,:])
#                            #print(np.abs(self.contr_qpt[0,54+dd,27,q,0]-self.reduced_contr_qpt[0,9,27,q,0]))
#                            #print(self.contr_qpt[0,54+dd,27,q,0])
#                    print(np.sum(self.reduced_contr_qpt[0,5,27+aa,:,-1]))
#                    print(self.reduced_td_ren[0,5,27+aa,-1])

            if self.senergy:
                self.reduced_self_energy= np.zeros((self.nsppol, self.reduced_nkpt, self.max_band, self.nfreq,2))
                self.reduced_self_energy[:,:,:,:,0] = self.average_kpts(self.reduced_kpts_index, self.self_energy[:,:,:,:,0], self.nfreq)
                self.reduced_self_energy[:,:,:,:,1] = self.average_kpts(self.reduced_kpts_index, self.self_energy[:,:,:,:,1], self.nfreq)

            if self.split_occupied_subspace:
               self.reduced_fan_occ = self.average_kpts(self.reduced_kpts_index, self.fan_occ, self.ntemp) 
               self.reduced_fan_unocc = self.average_kpts(self.reduced_kpts_index, self.fan_unocc, self.ntemp) 
               self.reduced_ddw_occ = self.average_kpts(self.reduced_kpts_index, self.ddw_occ, self.ntemp) 
               self.reduced_ddw_unocc = self.average_kpts(self.reduced_kpts_index, self.ddw_unocc, self.ntemp) 



        
    ###########

    def read_files(self):
        
        for f in (self.epc, self.out):
            if f.fname:
                f.read_nc()

    def set_rootname(self, root):
        # Set root for output names
        self.rootname = root
        if self.reduce_path :
            self.rootname = self.rootname+'_reduced'
                
    def get_corrections(self, npoints, eig0, ren):

        # compute corrected eigenvalues for each temperature
        arr = np.zeros((self.nsppol, npoints, self.max_band, self.ntemp))
        for ispin in range(self.nsppol):
            for ikpt in range(npoints):
                for iband in range(self.max_band):
                    for itemp in range(self.ntemp):

                        arr[ispin, ikpt, iband, itemp] = (eig0[ispin, ikpt, iband] + ren[ispin, ikpt, iband, itemp])             

        return arr

    def reshape_symmat(self, arr):
        # Take the symrel array and reshape it in (len(symrel))/9 3x3 matrices
        matrice = np.zeros((self.nsym,3,3))
        arr = np.reshape(arr, (self.nsym, 9))
        for i in range(self.nsym):
            matrice[i,:,:] = np.reshape(arr[i],(3,3))
        
        return matrice

    def reduce_symmetry_equivalent(self, arr):
        # Check if a given kpoint is equivalent by symmetry to a previous one in the list.
        # Returns a list of the indices of the symmetry-inequivalent kpoints in self.kpoints
        indices = []
        keep = []
        matrice = self.symmat
        
        for ikpt, kpt in enumerate(self.kpoints):
            
            if ikpt == 0:
                indices.append([ikpt])
                keep.append(kpt)
            
            else:
                found = False

                for imat in matrice:
                    if found is False:
                        # Apply the inverse symmetry matrix to kpoint
                        imat_inv = np.linalg.inv(imat)
                        symkpt = np.einsum('ji,i->j',imat_inv,kpt)
    
                        for j in range(len(keep)):
                            arr = keep[j]
                            # Check if kpoint is already in keep
                            if np.array_equal(arr,symkpt) is True:
                                found = True
                                indices[j].append(ikpt)
                                
                if found is False:
                    indices.append([ikpt])
                    keep.append(kpt)
        
        # Construct the reduced_kpath like self.kpoints
        kpath = np.empty((len(indices), 3))
        
        for i in range(len(kpath)):
            kpath[i,:] = keep[i]
        
        return  len(indices), indices, kpath

    def average_kpts(self, index, corr, dim):
        # Create reduced array and average corrections on symmetry-equivalent kpoints
        if dim == 1:
            reduced_corr = np.zeros((self.nsppol, len(index), self.max_band))

            for i in range(len(index)):
                arr = index[i]
    
                for indice in arr:
                    for a in range(self.nsppol):
                        for b in range(self.max_band):
                            reduced_corr[a,i,b] += corr[a,indice,b]
                    
                reduced_corr[:,i,:] = reduced_corr[:,i,:]/len(arr)

        else:
            reduced_corr = np.zeros((self.nsppol, len(index), self.max_band, dim))
        
            for i in range(len(index)):
                arr = index[i]
    
                for indice in arr:
                    for a in range(self.nsppol):
                        for b in range(self.max_band):
                            for c in range(dim):
                                reduced_corr[a,i,b,c] += corr[a,indice,b,c]
    
                reduced_corr[:,i,:,:] = reduced_corr[:,i,:,:]/len(arr)

        return reduced_corr
        
    def average_kpts_contr(self, index, corr, dim):
        # Create reduced array and average corrections on symmetry-equivalent kpoints
        if dim == 1:
            reduced_corr = np.zeros((self.nsppol, len(index), self.nband_contr))

            for i in range(len(index)):
                arr = index[i]
    
                for indice in arr:
                    for a in range(self.nsppol):
                        for b in range(self.nband_contr):
                            reduced_corr[a,i,b] += corr[a,indice,b]
                    
                reduced_corr[:,i,:] = reduced_corr[:,i,:]/len(arr)

        else:
            reduced_corr = np.zeros((self.nsppol, len(index), self.nband_contr, dim))
        
            for i in range(len(index)):
                arr = index[i]
    
                for indice in arr:
                    for a in range(self.nsppol):
                        for b in range(self.nband_contr):
                            for c in range(dim):
                                reduced_corr[a,i,b,c] += corr[a,indice,b,c]
    
                reduced_corr[:,i,:,:] = reduced_corr[:,i,:,:]/len(arr)

        return reduced_corr
 
    def write_renormalization(self):

        # Write the renormalized array in a text file
        create_directory(self.eig_dat)

        with open(self.eig_dat, 'w') as f:

            if self.reduce_path:
                f.write("Corrected eigenvalues (eV) on symmetry-reduced path for {} kpoints ".format(self.reduced_nkpt))
                f.write("for temperatures :\n ")
#                if self.temperature:
                for T in self.temp:
                    f.write("{:>8.1f} ".format(T)+"K     ")
                f.write("\n")
                
                for ikpt, kpt in enumerate(self.reduced_kpath):
                    f.write("Kpt : {0[0]} {0[1]} {0[2]}\n".format(kpt))

                    if self.bands_to_print:
                         for band in self.bands_to_print:
                            for itemp in range(self.ntemp):

                                corr = self.reduced_eigcorr[0, ikpt, band-1, itemp]*cst.ha_to_ev
                                f.write("{:>12.8f}   ".format(corr))        
                            f.write("\n")
                       
                    else:
                        for iband in range(self.max_band):
                            for itemp in range(self.ntemp):

                                corr = self.reduced_eigcorr[0, ikpt, iband, itemp]*cst.ha_to_ev
                                f.write("{:>12.8f}   ".format(corr))        
                            f.write("\n")


            else:
                f.write("Corrected eigenvalues (eV) for {} kpoints ".format(self.nkpt))
                f.write("for temperatures :\n ")
#                if self.temperature:
                for T in self.temp:
                    f.write("{:>8.1f} ".format(T)+"K     ")
                f.write("\n")
                
                for ikpt, kpt in enumerate(self.kpoints):
                    f.write("Kpt : {0[0]} {0[1]} {0[2]}\n".format(kpt))

                    if self.bands_to_print:
                        for iband in range(self.max_band):
                            for itemp in range(self.ntemp):

                                corr = self.eigcorr[0, ikpt, iband, itemp]*cst.ha_to_ev
                                f.write("{:>12.8f}   ".format(corr))        
                            f.write("\n")

                    else:
                        for iband in range(self.max_band):
                            for itemp in range(self.ntemp):

                                corr = self.eigcorr[0, ikpt, iband, itemp]*cst.ha_to_ev
                                f.write("{:>12.8f}   ".format(corr))        
                            f.write("\n")


            if self.gap_ren :
                # Write direct gap renormalization info
                if self.reduce_path:
                    self.gap_loc, self.gap_energy, self.gap_renorm = self.get_gap_info(self.reduced_eig0, self.reduced_eigcorr,False)
                    if self.split:
                        self.gap_loc2, self.gap_energy2, self.gap_renorm2 = self.get_gap_info(self.reduced_eig0, self.reduced_eigcorr,True)
                    if self.spline:
                        kpts = self.kpoints_fine
                    else:
                        kpts = self.reduced_kpath

                    if self.modes:
                        self.gap_renorm_modes = self.get_mode_info(self.gap_loc,self.reduced_zp_ren_modes,False)

                        if self.split:
                            self.gap_renorm_modes_split = self.get_mode_info(self.gap_loc2,self.reduced_zp_ren_modes,True)

                    if self.indirect:
                        self.indirect_gap_loc, self.indirect_gap_energy, self.indirect_gap_renorm, self.indirect_gap_energy_band, self.indirect_gap_renorm_band = self.get_indirect_gap_info(self.reduced_eig0, self.reduced_eigcorr, False)

                        if self.modes:
                            self.indirect_gap_renorm_modes = self.get_indirect_mode_info(self.indirect_gap_loc, self.reduced_zp_ren_modes,False)

                        if self.split:
                            self.indirect_gap_loc2, self.indirect_gap_energy2, self.indirect_gap_renorm2, self.indirect_gap_energy_band2, self.indirect_gap_renorm_band2 = self.get_indirect_gap_info(self.reduced_eig0, self.reduced_eigcorr, True)

                            if self.modes:
                                self.indirect_gap_renorm_modes_split = self.get_indirect_mode_info(self.indirect_gap_loc2, self.reduced_zp_ren_modes, True)
                
                else:
                    # The non-reduced option does not have as many possibilities as the reduced one... anyway if I have only gamma, it will reduce to itself! So, i could simply
                    # remove this option ???
                    self.gap_loc, self.gap_energy, self.gap_renorm = self.get_gap_info(self.eig0, self.eigcorr,False)

                    if self.spline:
                        kpts = self.kpoints_fine
                    else:
                        kpts = self.kpoints

                create_directory(self.gap_dat)

                with open(self.gap_dat, 'w') as f:
            
                    f.write('Direct gap location and renormalization :\n\n')
                    f.write('Unperturbed gap at kpoint {0[0]} {0[1]} {0[2]}'.format(kpts[self.gap_loc[0]]))
                    f.write(' with gap energy {:>10.8f} {}\n\n'.format(self.gap_energy[0],self.units)) 

                    for itemp, T in enumerate(self.temp):
        
                        f.write('At T={:>6.1f} K, '.format(self.temp[itemp]))
                        f.write('gap at kpoint {0[0]} {0[1]} {0[2]}'.format(kpts[self.gap_loc[itemp+1]]))
                        f.write(' with gap energy {:>10.8f} {}\n'.format(self.gap_energy[itemp+1],self.units))
                        f.write('Gap renormalization : {:>10.8f} {}\n\n'.format(self.gap_renorm[itemp],self.units))

                    if self.indirect:
                        f.write('--------------------------------------------------------------\n\n') 
                        f.write('Indirect gap location and renormalization :\n\n')
                        f.write('Unperturbed gap at kpoints {0[0]} {0[1]} {0[2]} (VB)  '.format(kpts[self.indirect_gap_loc[0,0]]))
                        f.write('{0[0]} {0[1]} {0[2]} (CB)'.format(kpts[self.indirect_gap_loc[0,1]]))
                        f.write(' with gap energy {:>10.8f} {}\n\n'.format(self.indirect_gap_energy[0],self.units)) 

                        for itemp, T in enumerate(self.temp):
            
                            f.write('At T={:>6.1f} K, '.format(self.temp[itemp]))
                            f.write('gap at kpoints {0[0]} {0[1]} {0[2]} (VB)  '.format(kpts[self.indirect_gap_loc[itemp+1,0]]))
                            f.write('{0[0]} {0[1]} {0[2]} (CB)'.format(kpts[self.indirect_gap_loc[itemp+1,1]]))
                            f.write(' with gap energy {:>10.8f} {}\n'.format(self.indirect_gap_energy[itemp+1],self.units))
                            f.write('Gap renormalization : {:>10.8f} {}\n'.format(self.indirect_gap_renorm[itemp],self.units))
                            f.write('VB max renormalization : {:>10.8f} {}\n'.format(self.indirect_gap_renorm_band[itemp,0], self.units))
                            f.write('CB min renormalization : {:>10.8f} {}\n\n'.format(self.indirect_gap_renorm_band[itemp,1], self.units))

                if self.split:
                    create_directory(self.split_gap_dat)
    
                    with open(self.split_gap_dat, 'w') as f:
                
                        for d in range(2):
                            if d==0:
                                f.write('***************** H-A path *******************\n\n')
                            else:
                                f.write('***************** A-L path *******************\n\n')

                            f.write('Direct gap location and renormalization :\n\n')
                            f.write('Unperturbed gap at kpoint {0[0]} {0[1]} {0[2]}'.format(kpts[self.gap_loc2[0,d]]))
                            f.write(' with gap energy {:>10.8f} {}\n\n'.format(self.gap_energy2[0,d],self.units)) 
        
                            for itemp, T in enumerate(self.temp):
                
                                f.write('At T={:>6.1f} K, '.format(self.temp[itemp]))
                                f.write('gap at kpoint {0[0]} {0[1]} {0[2]}'.format(kpts[self.gap_loc2[itemp+1,d]]))
                                f.write(' with gap energy {:>10.8f} {}\n'.format(self.gap_energy2[itemp+1,d],self.units))
                                f.write('Gap renormalization : {:>10.8f} {}\n\n'.format(self.gap_renorm2[itemp,d],self.units))
        
                            if self.indirect:
                                f.write('--------------------------------------------------------------\n\n') 
                                f.write('Indirect gap location and renormalization :\n\n')
                                f.write('Unperturbed gap at kpoints {0[0]} {0[1]} {0[2]} (VB)  '.format(kpts[self.indirect_gap_loc2[0,d,0]]))
                                f.write('{0[0]} {0[1]} {0[2]} (CB)'.format(kpts[self.indirect_gap_loc2[0,d,1]]))
                                f.write(' with gap energy {:>10.8f} {}\n\n'.format(self.indirect_gap_energy2[0,d],self.units)) 
        
                                for itemp, T in enumerate(self.temp):
                    
                                    f.write('At T={:>6.1f} K, '.format(self.temp[itemp]))
                                    f.write('gap at kpoints {0[0]} {0[1]} {0[2]} (VB)  '.format(kpts[self.indirect_gap_loc2[itemp+1,d,0]]))
                                    f.write('{0[0]} {0[1]} {0[2]} (CB)'.format(kpts[self.indirect_gap_loc2[itemp+1,d,1]]))
                                    f.write(' with gap energy {:>10.8f} {}\n'.format(self.indirect_gap_energy2[itemp+1,d],self.units))
                                    f.write('Gap renormalization : {:>10.8f} {}\n'.format(self.indirect_gap_renorm2[itemp,d],self.units))
                                    f.write('VB max renormalization : {:>10.8f} {}\n'.format(self.indirect_gap_renorm_band2[itemp,d,0], self.units))
                                    f.write('CB min renormalization : {:>10.8f} {}\n\n'.format(self.indirect_gap_renorm_band2[itemp,d,1], self.units))


    def get_param_index(self, arr, param):

        # Retreives the index of a specific value in an array. Here intended for finding a kpoint index.
        for x in range(len(arr)):
            if np.array_equal(param, arr[x]):
                return  x
        print('Param not found in array')
        return None

    def get_gap_info(self,unpert,pert,split):
        # Gets direct gap information

        # Finds the location of the unperturbed and perturbed direct gap + computes the gap energy
        val = self.valence - 1
        cond = val + 1

        if self.spline:
            # Interpolate kpoints
            self.interpolate_kpoints()
            # Interpolate bare eigenvalues
            unpert = self.interpolate_eigenvalues(unpert)
            # Interpolated eig(T)
            tmp = np.zeros((self.nsppol,self.nkpt_fine,self.max_band,self.ntemp))
            self.spline_reduced_td_ren = np.zeros_like(tmp)
            self.spline_reduced_eigcorr = np.zeros_like(tmp)

            for t in range(self.ntemp):
                tmp[:,:,:,t] = self.interpolate_eigenvalues(pert[:,:,:,t])
                self.spline_reduced_td_ren[:,:,:,t] = self.interpolate_corrections(self.reduced_td_ren[:,:,:,t])
                self.spline_reduced_eigcorr[:,:,:,t] = self.interpolate_corrections(self.reduced_eigcorr[:,:,:,t])

#            import matplotlib.pyplot as plt
#            fig, arr = plt.subplots(2,1,squeeze=True)
#
#            kptarr = np.arange(self.nkpt_fine)
#            col = ['gray','purple','blue','green','yellow','red']
#            for t in range(self.ntemp):
#                arr[0].plot(kptarr, self.spline_reduced_td_ren[0,:,cond,t]*cst.ha_to_ev*1E3, linestyle='dashed',color = col[t])
#                arr[1].plot(kptarr, self.spline_reduced_td_ren[0,:,val,t]*cst.ha_to_ev*1E3, linestyle = 'dotted',color = col[t])
#                arr[0].plot(kptarr[::25], self.reduced_td_ren[0,:,cond,t]*cst.ha_to_ev*1E3, linestyle='None',marker='o',color = col[t])
#                arr[1].plot(kptarr[::25], self.reduced_td_ren[0,:,val,t]*cst.ha_to_ev*1E3, linestyle='None',marker='s',color = col[t])
#            #    arr[0].plot(kptarr[::25], self.reduced_eigcorr[0,:,cond,t]*cst.ha_to_ev*1E3, linestyle='None',marker='o',color = col[t])
#            #    arr[1].plot(kptarr[::25], self.reduced_eigcorr[0,:,val,t]*cst.ha_to_ev*1E3, linestyle='None',marker='s',color = col[t])
#
#            for i in range(2):
#                lim = arr[i].get_ylim()
#                arr[i].plot(262*np.ones((10)), np.linspace(lim[0],lim[1],10),'k')
#            plt.suptitle('5GPa')
#
#            #plt.savefig('quadratic_spline_for_eigcorr_3gpa.png')
#            #plt.savefig('cubic_spline_for_eigcorr_3gpa.png')
#            plt.savefig('cubic_spline_for_correction_5gpa.png')
#            #plt.savefig('quadratic_spline_for_correction_3gpa.png')
#
#
#
#            plt.show()


            pert = tmp

        if split == False:

#            print(pert[:,:,27:28,:])
#            print(unpert[:,:,28])
#            print(unpert[:,:,27])
#            print(self.td_ren[:,:,27:28,:])

            loc = np.zeros((self.ntemp+1),dtype=int)
            ener = np.zeros((self.ntemp+1), dtype=float)
            ren = np.zeros((self.ntemp), dtype=float)
        
            arr0 = unpert[0,:,cond] - unpert[0,:,val]
            loc[0] = np.argmin(arr0)
            ener[0] = np.amin(arr0)*cst.ha_to_ev
              
            for i in range(self.ntemp):
           
                arr = pert[0,:,cond,i]-pert[0,:,val,i]
                loc[i+1] = np.argmin(arr)
                ener[i+1] = np.amin(arr)*cst.ha_to_ev 
                ren[i] = ener[i+1]-ener[0]

#                test = (unpert[0,:,cond]-unpert[0,:,val]+self.reduced_td_ren[0,:,cond,i]-self.reduced_td_ren[0,:,val,i])
#                print(arr-test)
#                print(arr)
#                print(np.argmin(arr))
#                print(np.amin(arr))

        
            if self.units == 'meV':
                return loc, ener*1000, ren*1000
            elif self.units == 'eV':
                return loc, ener, ren

        elif split==True:
            loc = np.zeros((self.ntemp+1,2),dtype=int)
            ener = np.zeros((self.ntemp+1,2), dtype=float)
            ren = np.zeros((self.ntemp,2), dtype=float)
        
            # Define mid in kpt array
            if not self.split_kpt:
                raise Exception('split_kpt is not well-defined, as reduced_kpt_array cannot be splitted. Please check your input')
#            if self.spline:
            '''FIX ME does not work with all spline_factors. Sometimes it returns kmid=None... '''
            kmid = self.find_kpt_index(self.split_kpt)
 #           else:
 #               kmid = self.get_param_index(self.reduced_kpath, self.split_kpt)
#            kmid = int(self.reduced_nkpt/2.)+1

            #Treat left gap
            arr0 = unpert[0,:kmid,cond] - unpert[0,:kmid,val]
            loc[0,0] = np.argmin(arr0)
            ener[0,0] = np.amin(arr0)*cst.ha_to_ev
            diff = (unpert[:,:kmid,28]-unpert[:,:kmid,27])

             
            for i in range(self.ntemp):
                arr = pert[0,:kmid,cond,i]-pert[0,:kmid,val,i]
                loc[i+1,0] = np.argmin(arr)
                ener[i+1,0] = np.amin(arr)*cst.ha_to_ev 
                ren[i,0] = ener[i+1,0]-ener[0,0]
#                test = (unpert[0,:kmid,cond]-unpert[0,:kmid,val]+self.reduced_td_ren[0,:kmid,cond,i]-self.reduced_td_ren[0,:kmid,val,i])
#                print(np.argmin(arr))
#                print(ener[i+1,0])
       
            #Treat right gap
            arr0 = unpert[0,kmid:,cond] - unpert[0,kmid:,val]
            loc[0,1] = np.argmin(arr0)+kmid
            ener[0,1] = np.amin(arr0)*cst.ha_to_ev

            
            for i in range(self.ntemp):
                arr = pert[0,kmid:,cond,i]-pert[0,kmid:,val,i]
                loc[i+1,1] = np.argmin(arr)+kmid
                ener[i+1,1] = np.amin(arr)*cst.ha_to_ev 
                ren[i,1] = ener[i+1,1]-ener[0,1]
 #               test = (unpert[0,kmid:,cond]-unpert[0,kmid:,val]+self.reduced_td_ren[0,kmid:,cond,i]-self.reduced_td_ren[0,kmid:,val,i])
 #               print(np.argmin(arr)+kmid)
  #              print(ener[i+1,1])


            if self.units == 'meV':
                return loc, ener*1000, ren*1000
            elif self.units == 'eV':
                return loc, ener, ren

    def get_indirect_gap_info(self,unpert,pert,split):
        # Gets indirect gap information

        # Finds the location of the unperturbed and perturbed direct gap + computes the gap energy
        #mid = np.int(len(self.bands_to_print)/2)
        #val = self.bands_to_print[mid-2]
        #cond = self.bands_to_print[mid-1]
        #val = self.bands_to_print[0]-1
        #cond = self.bands_to_print[1]-1
        val = self.valence - 1
        cond = val + 1

        if self.spline:
            # Interpolate kpoints
            self.interpolate_kpoints()
            # Interpolate bare eigenvalues
            unpert = self.interpolate_eigenvalues(unpert)
            # Interpolated eig(T)
            tmp = np.zeros((self.nsppol,self.nkpt_fine,self.max_band,self.ntemp))
            for t in range(self.ntemp):
                tmp[:,:,:,t] = self.interpolate_eigenvalues(pert[:,:,:,t])

            pert = tmp


        if split == False:

#            print(pert[:,:,27:28,:])
#            print(unpert[:,:,28])
#            print(unpert[:,:,27])
#            print(self.td_ren[:,:,27:28,:])

            loc = np.zeros((self.ntemp+1,2),dtype=int)
            ener = np.zeros((self.ntemp+1), dtype=float)
            ener_band = np.zeros((self.ntemp+1,2), dtype=float)
            ren = np.zeros((self.ntemp), dtype=float)
            ren_band = np.zeros((self.ntemp,2), dtype=float)

            arrc = unpert[0,:,cond]
            arrv = unpert[0,:,val]
            loc[0,0] = np.argmax(arrv)
            loc[0,1] = np.argmin(arrc)
            ener_band[0,0] = arrv[loc[0,0]]*cst.ha_to_ev
            ener_band[0,1] = arrc[loc[0,1]]*cst.ha_to_ev 
            ener[0] = ener_band[0,1] - ener_band[0,0]
              
            for i in range(self.ntemp):
           
                arrc = pert[0,:,cond,i]
                arrv = pert[0,:,val,i]
                loc[i+1,0] = np.argmax(arrv)
                loc[i+1,1] = np.argmin(arrc)
                ener_band[i+1,0] = arrv[loc[i+1,0]]*cst.ha_to_ev 
                ener_band[i+1,1] = arrc[loc[i+1,1]]*cst.ha_to_ev 
                ener[i+1] = ener_band[i+1,1] - ener_band[i+1,0] 
                for j in range(2):
                    ren_band[i,j] = ener_band[i+1,j] - ener_band[0,j]
                ren[i] = ener[i+1]-ener[0]

#                test = (unpert[0,:,cond]-unpert[0,:,val]+self.reduced_td_ren[0,:,cond,i]-self.reduced_td_ren[0,:,val,i])
#                print(arr-test)
#                print(arr)
#                print(np.argmin(arr))
#                print(np.amin(arr))

        
            if self.units == 'meV':
                return loc, ener*1000, ren*1000, ener_band*1000, ren_band*1000
            elif self.units == 'eV':
                return loc, ener, ren, ener_band, ren_band

        elif split==True:
            loc = np.zeros((self.ntemp+1,2,2),dtype=int)
            ener = np.zeros((self.ntemp+1,2), dtype=float)
            ener_band = np.zeros((self.ntemp+1,2,2), dtype=float)
            ren = np.zeros((self.ntemp,2), dtype=float)
            ren_band = np.zeros((self.ntemp,2,2), dtype=float)

        
            # Define mid in kpt array
            if self.spline:
                kmid = self.get_param_index(self.kpoints_fine, self.split_kpt)
            else:
                kmid = self.get_param_index(self.reduced_kpath, self.split_kpt)

#            kmid = int(self.reduced_nkpt/2.)+1

            #Treat left gap
            arrc = unpert[0,:kmid,cond]
            arrv = unpert[0,:kmid,val]
            loc[0,0,0] = np.argmax(arrv)
            loc[0,0,1] = np.argmin(arrc)
            ener_band[0,0,0] = arrv[loc[0,0,0]]*cst.ha_to_ev
            ener_band[0,0,1] = arrc[loc[0,0,1]]*cst.ha_to_ev 
            ener[0,0] = ener_band[0,0,1] - ener_band[0,0,0]
             
            for i in range(self.ntemp):
                arrc = pert[0,:kmid,cond,i]
                arrv = pert[0,:kmid,val,i]
                loc[i+1,0,0] = np.argmax(arrv)
                loc[i+1,0,1] = np.argmin(arrc)
                ener_band[i+1,0,0] = arrv[loc[i+1,0,0]]*cst.ha_to_ev 
                ener_band[i+1,0,1] = arrc[loc[i+1,0,1]]*cst.ha_to_ev 
                ener[i+1,0] = ener_band[i+1,0,1] - ener_band[i+1,0,0] 
                for j in range(2):
                    ren_band[i,0,j] = ener_band[i+1,0,j] - ener_band[0,0,j]
                ren[i,0] = ener[i+1,0]-ener[0,0]


            #Treat right gap
            arrc = unpert[0,kmid:,cond]
            arrv = unpert[0,kmid:,val]
            loc[0,1,0] = np.argmax(arrv)
            loc[0,1,1] = np.argmin(arrc)
            ener_band[0,1,0] = arrv[loc[0,1,0]]*cst.ha_to_ev
            ener_band[0,1,1] = arrc[loc[0,1,1]]*cst.ha_to_ev 
            ener[0,1] = ener_band[0,1,1] - ener_band[0,1,0]
            loc[0,1,0] += kmid
            loc[0,1,1] += kmid
            
            for i in range(self.ntemp):
                arrc = pert[0,kmid:,cond,i]
                arrv = pert[0,kmid:,val,i]
                loc[i+1,1,0] = np.argmax(arrv)
                loc[i+1,1,1] = np.argmin(arrc)
                ener_band[i+1,1,0] = arrv[loc[i+1,1,0]]*cst.ha_to_ev 
                ener_band[i+1,1,1] = arrc[loc[i+1,1,1]]*cst.ha_to_ev 
                ener[i+1,1] = ener_band[i+1,1,1] - ener_band[i+1,1,0] 
                for j in range(2):
                    ren_band[i,1,j] = ener_band[i+1,1,j] - ener_band[0,1,j]
                ren[i,1] = ener[i+1,1]-ener[0,1]
                loc[i+1,1,0] += kmid
                loc[i+1,1,1] += kmid


            if self.units == 'meV':
                return loc, ener*1000, ren*1000, ener_band*1000, ren_band*1000
            elif self.units == 'eV':
                return loc, ener, ren, ener_band, ren_band

    def interpolate_kpoints(self):
        
        if self.reduce_path:
            nkpt = self.reduced_nkpt
            kpts_coarse = self.reduced_kpath
        else:
            nkpt = self.nkpt
            kpts_coarse = self.kpoints

        self.nkpt_fine = nkpt*self.spline_factor
        x = np.arange(nkpt)
        xfine = np.linspace(x[0],x[-1],self.nkpt_fine)
        kpts = np.zeros((self.nkpt_fine,3))

        for i in range(3):
            spl = interp1d(x, kpts_coarse[:,i], kind='linear')
            kpts[:,i] = spl(xfine)

        self.kpoints_fine = kpts

#        with open('kpts_fine_spline_pc2.txt','w') as f:
#
#            for k in range(self.nkpt_fine):
#                f.write('{:>9.7e}  {:>9.7e}  {:9.7e}\n'.format(kpts[k,0], kpts[k,1], kpts[k,2]))
#
#            f.close()
        return

    def interpolate_eigenvalues(self,eig0):

         if self.reduce_path:
             nkpt = self.reduced_nkpt
         else:
             nkpt = self.nkpt

         x = np.arange(nkpt)
         xfine = np.linspace(x[0],x[-1],self.nkpt_fine)

         if self.nsppol is None:
             self.nsppol = eig0.nsppol
         if self.max_band is None:
             self.max_band = eig0.max_band

         eig = np.zeros((self.nsppol,self.nkpt_fine,self.max_band))

         for j in range(self.max_band):
             spl = interp1d(x,eig0[0,:,j], kind='cubic')
             eig[0,:,j] = spl(xfine)

         return eig

    def interpolate_corrections(self,corr):

         if self.reduce_path:
             nkpt = self.reduced_nkpt
         else:
             nkpt = self.nkpt

         x = np.arange(nkpt)
         xfine = np.linspace(x[0],x[-1],self.nkpt_fine)

         if self.nsppol is None:
             self.nsppol = corr.nsppol
         if self.max_band is None:
             self.max_band = corr.max_band

         eig = np.zeros((self.nsppol,self.nkpt_fine,self.max_band))

         aindex = self.find_kpt_index(self.split_kpt)

         for j in range(self.max_band):
             #spl = interp1d(x[:11],corr[0,:11,j], kind='quadratic')
             spl = interp1d(x[:11],corr[0,:11,j], kind='cubic')
             eig[0,:aindex+1,j] = spl(xfine[:aindex+1])
             #spl = interp1d(x[10:],corr[0,10:,j], kind='quadratic')
             spl = interp1d(x[10:],corr[0,10:,j], kind='cubic')
             eig[0,aindex+1:,j] = spl(xfine[aindex+1:])

         return eig


    def find_kpt_index(self, loc):

        # This locates a given kpoint in a list, within a given tolerance
        loclst = list(loc)
        if self.spline:
            lst = [list(x) for x in self.kpoints_fine]
        elif self.reduce_path:
            lst = [list(x) for x in self.reduced_kpath]
        else:
            lst = [list(x) for x in self.kpoints]

        for k,kpoint in enumerate(lst):
            index = np.allclose(loclst,kpoint)
            if index==True:
                print(kpoint)
                return k


    def get_mode_info(self,loc, mode_ren, split):
        # Gets direct gap information, splitted into mode contribution

        # Finds the location of the unperturbed and perturbed direct gap + computes the gap energy
        val = self.valence - 1
        cond = val + 1

        if split == False:

            loc = loc[1]

            ren = np.zeros((self.nmodes,2), dtype=float) #mode, val/cond
            
            ren[:,0] = mode_ren[0,loc,val,:]*cst.ha_to_ev
            ren[:,1] = mode_ren[0,loc,cond,:]*cst.ha_to_ev
        
            if self.units == 'meV':
                return ren*1000
            elif self.units == 'eV':
                return ren

        elif split==True:

            loc = loc[1,:] # T=0, left/right

            ren = np.zeros((self.nmodes,2,2), dtype=float) #mode, val/cond, left/right 
        
            #Treat left and right gap
             
            for i in range(2):

                ren[:,0,i] = mode_ren[0,loc[i],val,:]*cst.ha_to_ev
                ren[:,1,i] = mode_ren[0,loc[i],cond,:]*cst.ha_to_ev
       
            if self.units == 'meV':
                return ren*1000
            elif self.units == 'eV':
                return ren

    def get_indirect_mode_info(self,loc,mode_ren,split):
        # Gets indirect gap information, splitted into mode contribution

        # Finds the location of the unperturbed and perturbed direct gap + computes the gap energy
        val = self.valence - 1
        cond = val + 1

        if split == False:

            loc = loc[1,:] #T=0, val/cond

            ren = np.zeros((self.nmodes,2), dtype=float) #mode, val/cond
            
            ren[:,0] = mode_ren[0,loc[0],val,:]*cst.ha_to_ev
            ren[:,1] = mode_ren[0,loc[1],cond,:]*cst.ha_to_ev
        
            if self.units == 'meV':
                return ren*1000
            elif self.units == 'eV':
                return ren

        elif split==True:

            loc = loc[1,:,:] # T=0, left/right, val/cond

            ren = np.zeros((self.nmodes,2,2), dtype=float) #mode, val/cond, left/right 
        
            #Treat left and right gap
             
            for i in range(2):

                ren[:,0,i] = mode_ren[0,loc[i,0],val,:]*cst.ha_to_ev
                ren[:,1,i] = mode_ren[0,loc[i,1],cond,:]*cst.ha_to_ev
       
            if self.units == 'meV':
                return ren*1000
            elif self.units == 'eV':
                return ren


              

    def write_corrections(self):
        # Write the renormalized array in a text file
        create_directory(self.zpr_dat)

        with open(self.zpr_dat, 'w') as f:

            if self.reduce_path:
                f.write("Eigenvalue corrections (eV) on symmetry-reduced path for {} kpoints ".format(self.reduced_nkpt))
                f.write("for temperatures :\n ")

                for T in self.temp:
                    f.write("{:>8.1f} ".format(T)+"K     ")
                f.write("\n")
                
                for ikpt, kpt in enumerate(self.reduced_kpath):
                    f.write("Kpt : {0[0]} {0[1]} {0[2]}\n".format(kpt))

                    if self.bands_to_print:
                         for band in self.bands_to_print:
                            for itemp in range(self.ntemp):

                                corr = self.reduced_td_ren[0, ikpt, band-1, itemp]*cst.ha_to_ev
                                f.write("{:>12.8f}   ".format(corr))        
                            f.write("\n")
                       
                    else:
                        for iband in range(self.max_band):
                            for itemp in range(self.ntemp):

                                corr = self.reduced_td_ren[0, ikpt, iband, itemp]*cst.ha_to_ev
                                f.write("{:>12.8f}   ".format(corr))        
                            f.write("\n")


            else:
                f.write("Eigenvalue corrections (eV) for {} kpoints ".format(self.nkpt))
                f.write("for temperatures :\n ")

                for T in self.temp:
                    f.write("{:>8.1f} ".format(T)+"K     ")
                f.write("\n")
                
                for ikpt, kpt in enumerate(self.kpoints):
                    f.write("Kpt : {0[0]} {0[1]} {0[2]}\n".format(kpt))

                    if self.bands_to_print:
                        for band in self.bands_to_print:
                            for itemp in range(self.ntemp):

                                corr = self.td_ren[0, ikpt, band-1, itemp]*cst.ha_to_ev
                                f.write("{:>12.8f}   ".format(corr))        
                            f.write("\n")

                    else:
                        for iband in range(self.max_band):
                            for itemp in range(self.ntemp):

                                corr = self.td_ren[0, ikpt, iband, itemp]*cst.ha_to_ev
                                f.write("{:>12.8f}   ".format(corr))        
                            f.write("\n")


        if self.reduce_path:
            if self.points_to_print is not None:
                create_directory(self.corr_dat)
    
                with open(self.corr_dat, 'w') as g:
    
                    g.write("Eigenvalue corrections (meV) on symmetry-reduced path, for temperatures :\n")
                    for T in self.temp:
                        g.write("{:>8.1f} ".format(T)+"K     ")
                    g.write("\n")
                    
                    for index in self.points_to_print:
                        g.write("Kpt : {0[0]} {0[1]} {0[2]}\n".format(self.reduced_kpath[index,:]))
    
                        if self.bands_to_print:
                            mid = np.int(0.5*len(self.bands_to_print))
                            band1 = self.bands_to_print[mid]-1
                            for itemp in range(self.ntemp):
                                corr = (self.reduced_td_ren[0, index, band1, itemp]-self.reduced_td_ren[0,index,band1-1,itemp])*cst.ha_to_ev*1E3
                                g.write("{:>12.8f}   ".format(corr))        
                            g.write("\n")
                           
                        else:
                            raise Exception('Please specify which 2 bands must be considered') 

    def write_contribution(self):
        
        # Write the renormalized array in a text file
        create_directory(self.contr_dat)

     #   with open(self.contr_dat, 'w') as f:

            # for bands_to_print, write the following : eigcorr(k,n,t), sum(contr_qpt on qpt))

    def write_netcdf(self):
        # Write corrected array to netCDF file
        create_directory(self.nc_output)

        # Write on a NC file with etsf-io name convention
        with nc.Dataset(self.nc_output, "w") as dts:

            # Add gridsize as an attribute
            #dts.gridsize = self.grid

            # Create dimensions
            dts.createDimension('number_of_spins', self.nsppol)
            dts.createDimension('number_of_atoms', self.natom) 
            dts.createDimension('number_of_kpoints', self.nkpt)
            if self.spline:
                dts.createDimension('number_of_fine_kpoints',self.nkpt_fine)
            dts.createDimension('number_of_bands', self.max_band)
            dts.createDimension('number_of_bands_contr', self.nband_contr)
            dts.createDimension('number_of_frequencies', self.nfreq)
            dts.createDimension('cartesian', 3)
            dts.createDimension('cplex', 2)
            dts.createDimension('one', 1)

            dts.createDimension('number_of_temperatures', self.ntemp)

            dts.createDimension('number_of_qpoints', self.nqpt)
            dts.createDimension('number_of_modes',3*self.natom)

            # Create and write variables
            data = dts.createVariable('spline_interpolation', 'i1')
            data[:] = self.spline

            ## Bare variables
            data = dts.createVariable('temperatures', 'd', ('number_of_temperatures'))
            data[:] = self.temp[:]

            data = dts.createVariable('smearing', 'd')
            data[:] = self.smearing

            data = dts.createVariable('reduced_coordinates_of_kpoints', 'd', ('number_of_kpoints', 'cartesian'))
            data[:,:] = self.kpoints[:,:]
        
            data = dts.createVariable('bare_eigenvalues', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands'))
            data.units = "Hartree"
            data[:,:,:] = self.eig0[:,:,:]

            ## Frequencies for self-energy and spectral function
            data = dts.createVariable('omega_se', 'd', ('number_of_frequencies'))
            data[:] = self.omega_se[:]


            ## ZP- and TD- renormalization
            data = dts.createVariable('eigenvalue_corrections', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands', 'number_of_temperatures'))
            data.units = "Hartree"
            data[:,:,:,:] = self.td_ren[:,:,:,:]

            if self.spline:
                data = dts.createVariable('eigenvalue_corrections_spline_interpolation', 'd', ('number_of_spins','number_of_fine_kpoints','number_of_bands','number_of_temperatures'))
                data.units = 'Hartree'
                data[:,:,:,:] = self.spline_reduced_td_ren[:,:,:,:]

                data = dts.createVariable('corrected_eigenvalue_spline_interpolation', 'd', ('number_of_spins','number_of_fine_kpoints','number_of_bands','number_of_temperatures'))
                data.units = 'Hartree'
                data[:,:,:,:] = self.spline_reduced_eigcorr[:,:,:,:]


            data = dts.createVariable('eigenvalue_corrections_modes', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands', 'number_of_modes'))
            data.units = "Hartree"
            if self.modes:
                data[:,:,:,:] = self.zp_ren_modes[:,:,:,:]


            data = dts.createVariable('corrected_eigenvalues', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands', 'number_of_temperatures'))
            data.units = "Hartree"
            data[:,:,:,:] = self.eigcorr[:,:,:,:]


            ## Qpoint contribution
            data = dts.createVariable('qpt_contribution', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands', 'number_of_qpoints','number_of_temperatures'))
            if self.contribution and not self.modes:
                data[:,:,:,:,:] = self.contr_qpt[:,:,:,:,:]
            data = dts.createVariable('qpt_contribution_modes', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands', 'number_of_qpoints','number_of_modes'))
            if self.contribution and self.modes :
                data[:,:,:,:,:] = self.contr_qpt_modes[:,:,:,:,:]


            ## ZP Self-energy
            data = dts.createVariable('self_energy', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands', 'number_of_frequencies', 'cplex'))
            if self.senergy:
                data[:,:,:,:,:] = self.self_energy[:,:,:,:,:]            

            ## ZP Spectral function
            data = dts.createVariable('spectral_function', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands', 'number_of_frequencies'))
            if self.spectral:
                data[:,:,:,:] = self.spectral_function[:,:,:,:]

            ## ZP Broadening
            data = dts.createVariable('broadening', 'd', ('number_of_spins', 'number_of_kpoints', 'number_of_bands'))
            if self.broad:
                data[:,:,:] = self.broadening[:,:,:]

            ## Direct gap location, energy and renormalization 
            ## For now, only for nspins =1
            data = dts.createVariable('unperturbed_gap_location', 'd', ('cartesian'))   
            if self.gap_ren:
                if self.reduce_path:
                    if self.spline:
                        data[:] = self.kpoints_fine[self.gap_loc[0]]
                    else:
                        data[:] = self.reduced_kpath[self.gap_loc[0]]
                else:
                    if self.spline:
                        data[:] = self.kpoints_fine[self.gap_loc[0]]
                    else:
                        data[:] = self.kpoints[self.gap_loc[0]]

            data = dts.createVariable('gap_location', 'd',('number_of_temperatures', 'cartesian'))
            if self.gap_ren:
                if self.reduce_path:
                    if self.spline:
                        data[:,:] =  self.kpoints_fine[self.gap_loc[1:]]
                    else:
                        data[:,:] =  self.reduced_kpath[self.gap_loc[1:]]
                else:
                    if self.spline:
                        data[:,:] = self.kpoints_fine[self.gap_loc[1:]]
                    else:
                        data[:,:] = self.kpoints[self.gap_loc[1:]]

            data = dts.createVariable('unperturbed_gap_energy', 'd', ('one'))
            data.units = self.units
            if self.gap_ren:
                data[:] = self.gap_energy[0]

            data = dts.createVariable('gap_energy', 'd', ('number_of_temperatures'))
            data.units = self.units
            if self.gap_ren:
                data[:] = self.gap_energy[1:]

            data = dts.createVariable('gap_renormalization', 'd', ('number_of_temperatures'))
            data.units = self.units
            if self.gap_ren:
                data[:] = self.gap_renorm[:]

            # Indirect gap data
            data = dts.createVariable('unperturbed_indirect_gap_location', 'd', ('cplex','cartesian')) #(vb-cb, cart)  
            if self.gap_ren and self.indirect:
                if self.reduce_path:
                    if self.spline:
                        data[:,:] = self.kpoints_fine[self.indirect_gap_loc[0,:]]
                    else:
                        data[:,:] = self.reduced_kpath[self.indirect_gap_loc[0,:]]
#                else:
#                    data[:] = self.kpoints[self.gap_loc[0]]

            data = dts.createVariable('indirect_gap_location', 'd',('number_of_temperatures', 'cplex','cartesian'))
            if self.gap_ren:
                if self.indirect and self.reduce_path:
                    if self.spline:
                        data[:,:,:] = self.kpoints_fine[self.indirect_gap_loc[1:,:]]
                    else:
                        data[:,:,:] =  self.reduced_kpath[self.indirect_gap_loc[1:,:]]
#            else:
#                data[:,:] = self.kpoints[self.gap_loc[1:]]

            data = dts.createVariable('unperturbed_indirect_gap_energy', 'd', ('one'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.indirect:
                    data[:] = self.indirect_gap_energy[0]

            data = dts.createVariable('indirect_gap_energy', 'd', ('number_of_temperatures'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.indirect:
                    data[:] = self.indirect_gap_energy[1:]

            data = dts.createVariable('indirect_gap_renormalization', 'd', ('number_of_temperatures'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.indirect:
                    data[:] = self.indirect_gap_renorm[:]

            data = dts.createVariable('indirect_gap_energy_band', 'd', ('number_of_temperatures','cplex'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.indirect:
                    data[:,:] = self.indirect_gap_energy_band[1:,:]

            data = dts.createVariable('indirect_gap_renormalization_band', 'd', ('number_of_temperatures','cplex'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.indirect:
                    data[:,:] = self.indirect_gap_renorm_band[:,:]

            data = dts.createVariable('unperturbed_indirect_gap_location_split', 'd', ('cplex','cplex','cartesian'))   # (left-right, val-cond)
            if self.gap_ren:
                if self.reduce_path and self.split:
                    if self.spline:
                        data[:,:,:] = self.kpoints_fine[self.indirect_gap_loc2[0,:,:]]
                    else:
                        data[:,:,:] = self.reduced_kpath[self.indirect_gap_loc2[0,:,:]]
#                else:
#                    data[:] = self.kpoints[self.gap_loc[0]]

            data = dts.createVariable('indirect_gap_location_split', 'd',('number_of_temperatures', 'cplex', 'cplex','cartesian'))
            if self.gap_ren:
                if self.reduce_path and self.split:
                    if self.spline:
                        data[:,:,:,:] =  self.kpoints_fine[self.indirect_gap_loc2[1:,:,:]] #(temp, cart, lfet-right, val-cond)  
                    else:
                        data[:,:,:,:] =  self.reduced_kpath[self.indirect_gap_loc2[1:,:,:]] #(temp, cart, lfet-right, val-cond)
#            else:
#                data[:,:] = self.kpoints[self.gap_loc[1:]]

            data = dts.createVariable('unperturbed_indirect_gap_energy_split', 'd', ('one','cplex'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.split:
                    data[:,:] = self.indirect_gap_energy2[0,:]

            data = dts.createVariable('indirect_gap_energy_split', 'd', ('number_of_temperatures','cplex'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.split:
                    data[:,:] = self.indirect_gap_energy2[1:,:]

            data = dts.createVariable('indirect_gap_renormalization_split', 'd', ('number_of_temperatures','cplex'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.split:
                    data[:,:] = self.indirect_gap_renorm2[:,:]

            data = dts.createVariable('indirect_gap_energy_band_split', 'd', ('number_of_temperatures','cplex','cplex'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.split:
                    data[:,:,:] = self.indirect_gap_energy_band2[1:,:,:]

            data = dts.createVariable('indirect_gap_renormalization_band_split', 'd', ('number_of_temperatures','cplex','cplex'))
            data.units = self.units
            if self.gap_ren:
                if self.reduce_path and self.split:
                    data[:,:,:] = self.indirect_gap_renorm_band2[:,:,:]

            # Mode decomposition data
            data = dts.createVariable('gap_renormalization_by_modes','d', ('number_of_modes','cplex'))
            data.units = self.units
            if self.modes:
                data[:,:] = self.gap_renorm_modes[:,:]

            data = dts.createVariable('gap_renormalization_by_modes_split','d', ('number_of_modes','cplex','cplex'))
            data.units = self.units
            if self.modes and self.split:
                data[:,:,:] = self.gap_renorm_modes_split[:,:,:]

            data = dts.createVariable('indirect_gap_renormalization_by_modes','d', ('number_of_modes','cplex'))
            data.units = self.units
            if self.modes:
                data[:,:] = self.indirect_gap_renorm_modes[:,:]

            data = dts.createVariable('indirect_gap_renormalization_by_modes_split','d', ('number_of_modes','cplex','cplex'))
            data.units = self.units
            if self.modes and self.split:
                data[:,:,:] = self.indirect_gap_renorm_modes_split[:,:,:]


            ## Reduced path renormalization
#            if self.reduce_path:
            
            dts.createDimension('number_of_reduced_kpoints', self.reduced_nkpt)

            data = dts.createVariable('reduced_coordinates_of_reduced_kpath', 'd', ('number_of_reduced_kpoints', 'cartesian'))
            if self.reduce_path:
                data[:,:] = self.reduced_kpath[:,:]

            data = dts.createVariable('reduced_bare_eigenvalues', 'd', ('number_of_spins', 'number_of_reduced_kpoints', 'number_of_bands'))
            if self.reduce_path:
                data[:,:,:] = self.reduced_eig0[:,:,:]

            data = dts.createVariable('reduced_eigenvalue_corrections', 'd', ('number_of_spins', 'number_of_reduced_kpoints', 'number_of_bands', 'number_of_temperatures'))
            if self.reduce_path:
                data[:,:,:,:] = self.reduced_td_ren[:,:,:,:]

            data = dts.createVariable('reduced_eigenvalue_corrections_modes', 'd', ('number_of_spins', 'number_of_reduced_kpoints', 'number_of_bands', 'number_of_modes'))
            if self.reduce_path and self.modes:
                data[:,:,:,:] = self.reduced_zp_ren_modes[:,:,:,:]

            data = dts.createVariable('reduced_corrected_eigenvalues', 'd', ('number_of_spins', 'number_of_reduced_kpoints', 'number_of_bands', 'number_of_temperatures'))
            if self.reduce_path:
                data[:,:,:,:] = self.reduced_eigcorr[:,:,:,:]


            data = dts.createVariable('reduced_qpt_contribution', 'd', ('number_of_spins', 'number_of_reduced_kpoints', 'number_of_bands', 'number_of_qpoints','number_of_temperatures'))
            if self.contribution and not self.modes:
                data[:,:,:,:,:] = self.reduced_contr_qpt[:,:,:,:,:]

            data = dts.createVariable('reduced_qpt_contribution_modes', 'd', ('number_of_spins', 'number_of_reduced_kpoints', 'number_of_bands', 'number_of_qpoints','number_of_modes'))
            if self.contribution and self.modes:
                data[:,:,:,:,:] = self.reduced_contr_qpt_modes[:,:,:,:,:]

            data = dts.createVariable('reduced_self_energy', 'd', ('number_of_spins', 'number_of_reduced_kpoints', 'number_of_bands', 'number_of_frequencies', 'cplex'))
            if self.senergy:
                data[:,:,:,:,:] = self.reduced_self_energy[:,:,:,:,:]

            ## Split quantities for left/right gaps
            data = dts.createVariable('unperturbed_gap_location_split', 'd', ('cplex','cartesian'))   
            if self.gap_ren and self.split:
                if self.spline:
                    data[:,:] = self.kpoints_fine[self.gap_loc2[0,:]]
                else:
                    data[:,:] = self.reduced_kpath[self.gap_loc2[0,:]]
    
            data = dts.createVariable('gap_location_split', 'd',('number_of_temperatures', 'cplex', 'cartesian'))
            if self.split and self.gap_ren:
                if self.spline:
                    data[:,:,:] = self.kpoints_fine[self.gap_loc2[1:,:]]
                else:
                    data[:,:,:] =  self.reduced_kpath[self.gap_loc2[1:,:]]
    
            data = dts.createVariable('gap_energy_split', 'd', ('number_of_temperatures', 'cplex'))
            data.units = self.units
            if self.split and self.gap_ren:
                data[:,:] = self.gap_energy2[1:,:]
    
            data = dts.createVariable('gap_renormalization_split', 'd', ('number_of_temperatures','cplex'))
            data.units = self.units
            if self.split and self.gap_ren:
                data[:,:] = self.gap_renorm2[:,:]

            # Full contribution reduced data
            data = dts.createVariable('reduced_fan_g2','d',('number_of_spins','number_of_reduced_kpoints','number_of_bands_contr','number_of_bands_contr','number_of_modes','number_of_qpoints','cplex'))
            if self.reduce_path and self.full_contribution:
                data[:,:,:,:,:,:,:] = self.reduced_fan_g2[:,:,:,:,:,:,:]

            data = dts.createVariable('reduced_ddw_g2','d',('number_of_spins','number_of_reduced_kpoints','number_of_bands_contr','number_of_bands_contr','number_of_modes','number_of_qpoints','cplex'))
            if self.reduce_path and self.full_contribution:
                data[:,:,:,:,:,:,:] = self.reduced_ddw_g2[:,:,:,:,:,:,:]

            data = dts.createVariable('reduced_deltaE_ddw','d',('number_of_spins','number_of_reduced_kpoints','number_of_bands_contr','number_of_bands_contr','number_of_qpoints','cplex'))
            if self.reduce_path and self.full_contribution:
                data[:,:,:,:,:,:] = self.reduced_deltaE_ddw[:,:,:,:,:,:]

            data = dts.createVariable('reduced_fan_occterm', 'd',('number_of_spins','number_of_reduced_kpoints','number_of_bands_contr','number_of_bands_contr','number_of_modes','number_of_temperatures',
                'number_of_qpoints','cplex'))
            if self.reduce_path and self.full_contribution:
                data[:,:,:,:,:,:,:,:] = self.reduced_fan_occterm[:,:,:,:,:,:,:,:]

            data = dts.createVariable('qpt_weight','d',('number_of_qpoints'))
            if self.qpt_weight is not None:
                data[:] = self.qpt_weight[:]

            data = dts.createVariable('ddw_tdep','d',('number_of_modes','number_of_temperatures','number_of_qpoints'))
            if self.full_contribution:
                data[:,:,:] = self.ddw_tdep[:,:,:]

        
            # Contribution from Fan/DW, occ/unocc
            data = dts.createVariable('reduced_fan_occ','d', ('number_of_spins','number_of_reduced_kpoints','number_of_bands','number_of_temperatures'))
            data.units = "Hartree"
            if self.reduce_path and self.split_occupied_subspace:
                data[:,:,:,:] = self.reduced_fan_occ[:,:,:,:]

            data = dts.createVariable('reduced_fan_unocc','d', ('number_of_spins','number_of_reduced_kpoints','number_of_bands','number_of_temperatures'))
            data.units = "Hartree"
            if self.reduce_path and self.split_occupied_subspace:
                data[:,:,:,:] = self.reduced_fan_unocc[:,:,:,:]

            data = dts.createVariable('reduced_ddw_occ','d', ('number_of_spins','number_of_reduced_kpoints','number_of_bands','number_of_temperatures'))
            data.units = "Hartree"
            if self.reduce_path and self.split_occupied_subspace:
                data[:,:,:,:] = self.reduced_ddw_occ[:,:,:,:]

            data = dts.createVariable('reduced_ddw_unocc','d', ('number_of_spins','number_of_reduced_kpoints','number_of_bands','number_of_temperatures'))
            data.units = "Hartree"
            if self.reduce_path and self.split_occupied_subspace:
                data[:,:,:,:] = self.reduced_ddw_unocc[:,:,:,:]
            
        return

    def write_gap_info(self):
        print('tada')

    @property
    def eig_dat(self):
        return str('OUT/dat/'+self.rootname+'_EIG.dat')
    
    @property
    def zpr_dat(self):
        return str('OUT/dat/'+self.rootname+'_ZPR.dat')

    @property
    def gap_dat(self):
        return str('OUT/dat/'+self.rootname+'_GAP.dat')

    @property
    def split_gap_dat(self):
        return str('OUT/dat/'+self.rootname+'_SPLIT.dat')

    @property
    def corr_dat(self):
        return str('OUT/dat/'+self.rootname+'_CORR.dat')

    @property
    def nc_output(self):
        return str('OUT/nc/'+self.rootname+'_ZPR.nc')

    @property
    def contr_dat(self):
        return str('OUT/dat/'+self.rootname+'_CONTR.dat')
    

# Create a directory if it does not exist
def create_directory(fname):

    dirname = os.path.dirname(fname)
    if not dirname:
        return
    if not os.path.exists(dirname):
        os.system('mkdir -p ' + dirname)

# Main function
def compute(
        #Input files
        epc_fname = None,
        out_fname = None,
        rootname = 'eigcorr.out',

        #Parameters
        nsppol = None,
        nkpt = None,
        max_band = None,
        eig0 = None,
        units = 'eV',
        kpoints = None,
        zp_ren = None,
        td_ren = None,
        temp = None,
        ntemp = None,
        grid = None,
        nband_contr = None,
        valence = None,

        #Options
        bands_to_print = None,
        points_to_print = None,
        split_kpt = None,
        reduce_path = False,
        temperature = False,
        contribution = False,
        full_contribution = False,
        gap_ren = False,
        senergy = False,
        broad = False,
        spectral = False,
        indirect = False,
        modes=False,
        split_contribution = False,
        split_occupied_subspace = False,
        spline = False,
        spline_factor = 5,

        **kwargs):

    # Compute corrected bandstructures for ZPM and/or temperature-dependent corrections
    corrections = EigenvalueCorrections(#passed arguments)
            out_fname = out_fname, 
            epc_fname = epc_fname,

            rootname = rootname,
            grid = grid,
            bands_to_print = bands_to_print,
            reduce_path = reduce_path,
            contribution = contribution,
            gap_ren = gap_ren,
            temperature = temperature,
            senergy = senergy,
            broad = broad,
            spectral = spectral,
            indirect = indirect,
            modes = modes,
            full_contribution = full_contribution,
            split_contribution = split_contribution,
            split_occupied_subspace = split_occupied_subspace,

            valence = valence,
            nsppol = nsppol,
            nkpt = nkpt,
            max_band = max_band,
            eig0 = eig0,
            temp = temp,
            td_ren = td_ren,
            ntemp = ntemp,
            units = units,
            kpoints = kpoints,
            points_to_print = points_to_print,
            split_kpt = split_kpt,
            nband_contr = nband_contr,
            spline = spline,
            spline_factor = spline_factor,

            **kwargs)
    
    # Write output file
    corrections.write_renormalization()
    corrections.write_corrections()
    corrections.write_contribution()
    corrections.write_netcdf()

    return corrections


###########################


