#! /usr/bin/env python

from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import netCDF4 as nc
import numpy as np
import os
import warnings
import itertools as itt
from copy import copy

############################
# this code plots the corrected bandstructure on a symmetrized kpath obtained with zpr_analyser.py

# Inputfile : <rootname>_ZPR.nc
# Output : png image

# V. Brousseau 2017/07
###########################

# Tex settings
#rc('text', usetex = True)
#rc('font', family = 'serif', weight = 'bold')


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

class GAPfile(CDFfile):

    def __init__(self, *args, **kwargs):

        super(GAPfile, self).__init__(*args, **kwargs)

        self.explicit_gap_energies = None
        self.explicit_pressures = None

    def read_nc(self, fname=None):

        fname = fname if fname else self.fname
        super(GAPfile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as ncdata:
        
            self.explicit_gap_energies = ncdata.variables['gap_energies'][:]
            self.explicit_pressures = ncdata.variables['pressures'][:]
            self.explicit_gap_energies_units = ncdata.variables['gap_energies'].getncattr('units')

class EXPfile(CDFfile):

    def __init__(self, *args, **kwargs):

        super(EXPfile, self).__init__(*args,**kwargs)

        self.xaxis = None
        self.yaxis = None
        self.ndata = None

    def read_nc(self, fname=None):

        fname = fname if fname else self.fname
        super(EXPfile, self).read_nc(fname)

        with nc.Dataset(fname,'r') as ncdata:

            self.xaxis = ncdata.variables['ax1'][:]
            self.yaxis = ncdata.variables['ax2'][:]
            self.ndata = len(self.xaxis)
            self.xaxis_units = ncdata.variables['ax1'].getncattr('units')
            self.yaxis_units = ncdata.variables['ax2'].getncattr('units')


class ZPRfile(CDFfile):
        
    def __init__(self, *args, **kwargs):

        super(ZPRfile, self).__init__(*args, **kwargs)

        self.eig0 = None
        self.kpoints = None
        self.temp = None
        self.eigcorr = None
        self.omega_se = None
        self.self_energy = None
        self.spectral_function = None
        self.broadening = None
        self.smearing = None

    # Open the ZPR.nc file and read it
    def read_nc(self, fname = None):
        
        fname = fname if fname else self.fname
        super(ZPRfile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as ncdata:

            self.eig0 = ncdata.variables['reduced_bare_eigenvalues'][:,:,:]
            self.kpoints = ncdata.variables['reduced_coordinates_of_reduced_kpath'][:,:]
            self.temp = ncdata.variables['temperatures'][:]
            self.eigcorr = ncdata.variables['reduced_corrected_eigenvalues'][:,:,:,:]
            self.correction = ncdata.variables['reduced_eigenvalue_corrections'][:,:,:,:]
            self.omega_se = ncdata.variables['omega_se'][:]
            self.self_energy = ncdata.variables['reduced_self_energy'][:,:,:,:,:]
            self.spectral_function = ncdata.variables['spectral_function'][:,:,:,:]
            self.broadening = ncdata.variables['broadening'][:,:,:]
            self.smearing = ncdata.variables['smearing'][:]

#            self.gridsize = ncdata.gridsize
            self.unperturbed_gap_location = ncdata.variables['unperturbed_gap_location'][:]
            self.unperturbed_gap_energy = ncdata.variables['unperturbed_gap_energy'][:]
            self.gap_location = ncdata.variables['gap_location'][:,:]
            self.gap_energy = ncdata.variables['gap_energy'][:]
            self.gap_energy_units = ncdata.variables['gap_energy'].getncattr('units')
            self.gap_renormalization = ncdata.variables['gap_renormalization'][:]
#            self.band_energy = ncdata.variables['band_energy'][:,:]
#            self.band_renormalization = ncdata.variables['band_renormalization'][:,:]

            self.eigcorr_modes = ncdata.variables['eigenvalue_corrections_modes'][:,:,:,:]

            self.unperturbed_gap_location_split = ncdata.variables['unperturbed_gap_location_split'][:,:]
            #self.unperturbed_gap_energy_split = ncdata.variables['unperturbed_gap_energy_split'][:,:]
            self.gap_location_split = ncdata.variables['gap_location_split'][:,:,:]
            self.gap_energy_split = ncdata.variables['gap_energy_split'][:,:]
            self.gap_energy_units_split = ncdata.variables['gap_energy_split'].getncattr('units')
            self.gap_renormalization_split = ncdata.variables['gap_renormalization_split'][:,:]
#            self.band_energy_split = ncdata.variables['band_energy_split'][:,:,:]
#            self.band_renormalization_split = ncdata.variables['band_renormalization_split'][:,:,:]

            self.unperturbed_indirect_gap_location = ncdata.variables['unperturbed_indirect_gap_location'][:,:]
            self.indirect_gap_location = ncdata.variables['indirect_gap_location'][:,:,:]
            self.unperturbed_indirect_gap_energy = ncdata.variables['unperturbed_indirect_gap_energy'][:]
            self.indirect_gap_energy = ncdata.variables['indirect_gap_energy'][:]
            self.indirect_gap_ren = ncdata.variables['indirect_gap_renormalization'][:]
            self.indirect_gap_energy_band = ncdata.variables['indirect_gap_energy_band'][:,:]
            self.indirect_gap_ren_band = ncdata.variables['indirect_gap_renormalization_band'][:,:]

            self.unperturbed_indirect_gap_location_split = ncdata.variables['unperturbed_indirect_gap_location_split'][:,:,:]
            self.indirect_gap_location_split = ncdata.variables['indirect_gap_location_split'][:,:,:,:]
            self.unperturbed_indirect_gap_energy_split = ncdata.variables['unperturbed_indirect_gap_energy_split'][:,:]
            self.indirect_gap_energy_split = ncdata.variables['indirect_gap_energy_split'][:,:]
            self.indirect_gap_ren_split = ncdata.variables['indirect_gap_renormalization_split'][:,:]
            self.indirect_gap_energy_band_split = ncdata.variables['indirect_gap_energy_band_split'][:,:,:]
            self.indirect_gap_ren_band_split = ncdata.variables['indirect_gap_renormalization_band_split'][:,:,:]

            ### Only for split contribution, VB and CB / modes separate ###

#            self.fan_g2 = ncdata.variables['reduced_fan_g2'][:,:,:,:,:,:,:] # spin kpt 2 2 mode qpt cplex
#            self.ddw_g2 = ncdata.variables['reduced_ddw_g2'][:,:,:,:,:,:,:] # same
#            self.deltaE_ddw = ncdata.variables['reduced_deltaE_ddw'][:,:,:,:,:,:] # spin kpt 2 2 qpt cplex
#            self.fan_occterm = ncdata.variables['reduced_fan_occterm'][:,:,:,:,:,:,:,:] # spin kpt 2 2 mode temp qpt cplex
#            self.qpt_weight = ncdata.variables['qpt_weight'][:]
#            self.ddw_tdep = ncdata.variables['ddw_tdep'][:,:,:] # mode temp qpt


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
        return np.int(self.omega_se.shape[0]) if self.omega_se is not None else None

    @property
    def nmodes(self):
        if self.fan_g2 is not None:
            return self.fan_g2.shape[4] if self.fan_g2 is not None else None
        elif self.eigcorr_modes is not None:
            return self.eigcorr_modes.shape[0]
        else:
            return None

    @property
    def nqpt(self):
        return self.fan_g2.shape[-2] if self.fan_g2 is not None else None

class ZPR_plotter(object):

    # Input file
    zpr_fnames = list()
    gap_fname = None
        
    # Parameters
    nsppol = None
    nkpt = None
    max_band = None
    ntemp = None
    temp = None
    eig0 = None
    eigcorr = None
    kpoints = None
    fermi = None
    fermi_td = None
    tmax = None
    scissor = None
    units = 'eV'
    gap_units = None
    pressure = None
    crit_pressure = None
    explicit_pressures = None
    
    gap_loc0 = None
    gap_energy0 = None
    gap_loc = None
    gap_energy = None
    gap_ren = None

    omega_se = None
    spectral_function = None
    self_energy = None
    broadening = None
    smearing = None
    
    # Options
    figsize = (20,9)
    color = None
    labels = None
    linestyle = None
    ylims = None
    yticks = None
    yminorticks = None
    xlims = None
    xticks = None
    xticks_labels = None
    xticks_alignment = None
    bands_to_print = None
    band_numbers = None
    point_for_se = None
    temp_to_print = 0.
    subplots = False
    main_title = None
    title = None
    savefile = None
    degen = None

    renormalization = True
    spectral = False
    senergy = False
    broad = False
    gap = False
    pgap = False
    vbcb = False
    separate_bands = False
    split = False
    split2 = False
    follow = False
    indirect = False
    unpert = False
    split_contribution = False
    modes = False
    verbose = False
    zero_gap_value = None
    experimental_data = None
    expdata = None
    zero_gap_units = 'eV'

    def __init__(self,

            # Options
            figsize = (20,9),
            color = None,
            linestyle = None,
            labels = None,
            ylims = None,
            ylimits = None,
            yticks = None,
            yminorticks = None,
            xlims = None,
            xticks = None,
            xticks_labels = None,
            xticks_alignment = None,
            bands_to_print = None,
            band_numbers = None,
            point_for_se = None,
            temp_to_print = 0.,
            subplots = False,
            main_title = None,
            title = None,
            savefile = None,

            gap_loc0 = None,
            gap_energy0 = None,
            gap_loc = None,
            gap_energy = None,
            gap_ren = None,


            omega_se = None,
            spectral_function = None,
            self_energy = None,
            broadening = None,
            smearing = None,

            renormalization = True,
            spectral = False,
            senergy = False,
            broad = False,
            gap = False,
            vbcb = False,
            pgap = False,
            separate_bands = False,
            split = False,
            split2 = False,
            follow=False,
            indirect = False,
            unpert = False,
            split_contribution = False,
            modes = False,
            degen = None,
            verbose = False,
            zero_gap_value = None,
            zero_gap_units = 'eV',

            # Parameters
            nsppol = None,
            nkpt = None,
            max_band = None,
            ntemp = None,
            temp = None,
            eig0 = None,
            eigcorr = None,
            kpoints = None,
            fermi = None,
            fermi_td = None,
            tmax = None,
            scissor = None,
            pressure = None,
            crit_pressure = None,
            explicit_pressures = None,
            experimental_data = None,
            expdata = None,

            units = 'eV',
            gap_units = None,
            rootname = 'zpr.png',

            # Input file
            zpr_fnames = list(),
            gap_fname = None,
            
            **kwargs):

        # Set input file
        self.zpr_fnames = zpr_fnames
        self.gap_fname = gap_fname


        # Define options
        self.figsize = figsize
        self.color = color
        self.labels = labels
        self.linestyle = linestyle
        self.ylims = ylims
        self.yticks = yticks
        self.yminorticks = yminorticks
        self.xlims = xlims
        self.xticks = xticks
        self.xticks_labels = xticks_labels
        self.xticks_alignment = xticks_alignment
        self.bands_to_print = bands_to_print
        self.band_numbers = band_numbers
        self.point_for_se = point_for_se
        self.temp_to_print = temp_to_print
        self.units = units
        self.subplots = subplots
        self.main_title = main_title
        self.title = title
        self.savefile = savefile
        self.scissor = scissor
        self.fermi = fermi
        self.fermi_td = fermi_td
        self.tmax = tmax
        self.pressure = pressure
        self.crit_pressure = crit_pressure
        self.verbose = verbose
        self.zero_gap_value = zero_gap_value 
        self.zero_gap_units = zero_gap_units
        self.experimental_data = experimental_data

        self.renormalization = renormalization
        self.spectral = spectral
        self.senergy = senergy
        self.broad = broad
        self.gap = gap
        self.vbcb = vbcb
        self.pgap = pgap
        self.separate_bands = separate_bands
        self.split = split
        self.split2 = split2
        self.follow = follow
        self.indirect = indirect
        self.unpert = unpert
        self.split_contribution = split_contribution
        self.modes = modes

        # Check if input is correct 
        self.check_input()

        self.set_rootname(rootname)

        if self.experimental_data is not None:
            self.expdata = EXPfile(self.experimental_data)
            self.expdata.read_nc()
    
    def check_input(self):
        
        # Checks that all required input is provided
        if not self.zpr_fnames:
            raise Exception('Must provide files for zpr_fnames')


        if sum(map(bool,[self.renormalization, self.senergy, self.spectral])) > 1:
            raise Exception('Please specify ONLY ONE quantity to plot from renormalization, self-energy and spectral function')

        if self.senergy:
            if not self.point_for_se:
                raise Exception('Must provide a kpoint index to treat in point_for_se')  
    
            if not self.bands_to_print:
                raise Exception('Must provide band numbers for self-energy in bands_to_print') 

        if self.units is not 'meV' and self.units is not 'eV':
            raise Exception('Units must be eV or meV')


    def read_file(self):
    
        f = self.zpr

        if f.fname:
            f.read_nc()

    def read_other_file(self):
        f = self.exgap

        if f.fname:
            f.read_nc()

    def read_expdata_file(self):

        with open(self.experimental_data,'r') as f:

            for i in range(4):
                f.readline()

            ndata = f.readline()
            f.readline()

            data = np.empty((ndata,2))

            for i in range(ndata):
                line = f.readline().split(',')
                data[i,0] = line[0]
                data[i,1] = line[1]

        return data
            
            

       

    def set_rootname(self, root):
        self.rootname = root

    
    def plot_zpr(self):

        file_qty = len(self.zpr_fnames)

        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
    
            self.omega_se = self.zpr.omega_se
            self.self_energy = self.zpr.self_energy
            self.spectral_function = self.zpr.spectral_function
            self.broadening = self.zpr.broadening
            self.smearing = self.zpr.smearing

            kpt_array = np.arange(self.nkpt)
            scissor_array = np.zeros(self.max_band)

            if self.tmax is None:
                self.tmax = self.temp[-1]    
    
            # Define the scissor shift for each band
            if self.scissor is not None:
                for iband in range(self.max_band):
                    if iband >= self.scissor[0]:
                        scissor_array[iband] = self.scissor[1]

            if self.units is 'eV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev
            elif self.units is 'meV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev*1000
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev*1000
                for i in range(2):
                    self.ylims[i] = self.ylims[i]*1000
                if self.scissor is not None: 
                    self.scissor[1] = self.scissor[1]*1000


            # this is just a quick check for relative error when using the double grid. For a given kpt index and band index, compute omega_max/(Emax-Ekn)
#            omax = 0.02
#            err = np.zeros((5,2))
#            for ik,ib in itt.product(range(5),range(2)):
#                err[ik,ib] = omax/np.abs(self.eig0[0,ik+3,ib+27]-self.eig0[0,ik+3,-1])*100
#             
#            print(err)
            

            # define index for max temperature
            tmax_index = self.find_temp_index()
#            lines = [

            # Set parameters from 1st file            
            if ifile==0:
                # Plot thermal corrections on the same plot or on different subplots
                if self.subplots :
                    plot_qty = tmax_index+1
                else:
                    plot_qty = 1

                # Ajust Fermi to eV or meV
                if self.units =='eV':
#                    self.fermi = self.fermi*cst.ha_to_ev
                    for m in range(len(self.fermi_td)):
                        self.fermi_td[m] = self.fermi_td[m]*cst.ha_to_ev
                elif self.units == 'meV':
 #                   self.fermi = self.fermi*cst.ha_to_ev*1000
                    for m in range(len(self.fermi_td)):
                        self.fermi_td[m] = self.fermi_td[m]*cst.ha_to_ev*1000

                # Define subplots grid
                fig, arr = plt.subplots(plot_qty, file_qty, squeeze = False, sharey = True, figsize = self.figsize)

                if self.ntemp > len(self.color):
                    raise Exception('Color vector not long enough! {} colors required.'.format(self.ntemp))
                if self.fermi_td is None:
                    raise Exception('Must provide fermi_td values')
                if file_qty > len(self.fermi_td):
                    raise Exception('fermi_td vector not long enough! {} fermi energy values required.'.format(file_qty))

            
            for iplot in range(plot_qty):

                if self.bands_to_print:
                    for iband in self.bands_to_print:
                        arr[iplot][ifile].plot(kpt_array, self.eig0[0,:,iband-1]-self.fermi_td[ifile]+scissor_array[iband-1], color='k',linewidth=1.5)
                        
                        if self.subplots:
                            arr[iplot][ifile].plot(kpt_array, self.eigcorr[0,:,iband-1,iplot]-self.fermi_td[ifile]+scissor_array[iband-1], color=self.color[iplot],linewidth=1.5)
                        else:
                            for itemp in range(tmax_index+1):
                                arr[iplot][ifile].plot(kpt_array, self.eigcorr[0,:,iband-1,itemp]-self.fermi_td[ifile]+scissor_array[iband-1], color = self.color[itemp],linewidth=1.5)
        
                else:
                    for iband in range(self.max_band):
                        arr[iplot][ifile].plot(kpt_array, self.eig0[0,:,iband]-self.fermi_td[ifile]+scissor_array[iband], color='k', label='unperturbed')
                        
                        if self.subplots:
                            arr[iplot][ifile].plot(kpt_array, self.eigcorr[0,:,iband,iplot]-self.fermi_td[ifile]+scissor_array[iband], color = self.color[iplot], label=str(self.temp[iplot])+'K')
                        else:
                            for itemp in range(tmax_index+1):
                                arr[iplot][ifile].plot(kpt_array, self.eigcorr[0,:,iband,itemp]-self.fermi_td[ifile]+scissor_array[iband], color = self.color[itemp], label=str(self.temp[itemp])+'K')
    
                # Set xlims, xticks, xtickslabels
                self.set_xaxis(arr[iplot][ifile],kpt_array)
        
                # Set visual references
                self.set_vrefs(arr[iplot][ifile], kpt_array, 0.) 
    
                # set legend
                if ifile == file_qty-1:
                    self.set_legend_ren(arr[iplot][ifile],iplot)

                # Set subfigure title
                self.set_title(arr[iplot][ifile], self.title[ifile])
    
            # Set y-axis
            #self.set_yaxis(arr[0][ifile], 'E-E(A-CBM) ({})'.format(self.units))
        self.set_yaxis(arr[0][0], r'E-E$_F$ ({})'.format(self.units))

        plt.subplots_adjust(left=0.08,bottom=0.12,right=0.94,top=0.90,wspace=0.08,hspace=0.20)
#        self.set_legend_ren2(fig,lines)

#        self.set_main_title(fig)
        self.save_figure(fig)

        plt.show()


    def plot_vbcb(self):

        # Plot valence band and conduction band corrections on the given reduced kpath
        # Later, add bstr on top for visual reference?
        # Later, add total gap correction on the bottom?
        file_qty = len(self.zpr_fnames)

        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
            self.gap_location = self.zpr.gap_location
    
#            self.omega_se = self.zpr.omega_se
#            self.self_energy = self.zpr.self_energy
#            self.spectral_function = self.zpr.spectral_function
            self.broadening = self.zpr.broadening
            self.smearing = self.zpr.smearing

            kpt_array = np.arange(self.nkpt)
            scissor_array = np.zeros(self.max_band)
    
            # Define the scissor shift for each band
            if self.scissor is not None:
                for iband in range(self.max_band):
                    if iband >= self.scissor[0]:
                        scissor_array[iband] = self.scissor[1]

            if self.units is 'eV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev
                self.correction = self.zpr.correction*cst.ha_to_ev
            elif self.units is 'meV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev*1000
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev*1000
                self.correction = self.zpr.correction*cst.ha_to_ev*1000
                if self.ylims is not None:
                    self.ylims = self.ylims*1000
                if self.scissor is not None: 
                    self.scissor[1] = self.scissor[1]*1000


            # Define which temperature to print
            temp_index = self.find_temp_index()

            
            # Set parameters from 1st file            
            if ifile==0:

                # Ajust Fermi to eV or meV
                if self.units =='eV':
                    self.fermi = self.fermi*cst.ha_to_ev
                elif self.units == 'meV':
                    self.fermi = self.fermi*cst.ha_to_ev*1000

                # Define subplots grid
                fig, arr = plt.subplots(4, 1, squeeze = False, sharex =True, sharey = False, figsize = self.figsize)
                plt.subplots_adjust(hspace=0.05, top=0.95, bottom=0.05)

                # Create list of ylimits
                self.ylimits = []

                if file_qty > len(self.color):
                    raise Exception('Color vector not long enough! {} colors required.'.format(file_qty))
                if file_qty > len(self.labels):
                    raise Exception('Labels vector not long enough! {} labels required.'.format(file_qty))
                if self.linestyle:
                    if file_qty > len(self.linestyle):
                        raise Exception('Linestyle vector not long enough! {} linestyle required.'.format(file_qty))


                if not self.band_numbers:
                    raise Exception('Must precise valence and conduction band number!')

        
                self.main_title = self.main_title+' {}K'.format(self.temp[temp_index])

                gap_index = self.find_gap_index(self.gap_location[temp_index])
                A_index = self.find_gap_index(np.array([0.,0.,0.5]))


                
                # Plot unperturbed eigenvalues
                for iband in self.bands_to_print:
                    arr[0,0].plot(kpt_array, self.eig0[0,:,iband-1]-self.fermi+scissor_array[iband-1], color='k')


            if temp_index is not None:
                if self.linestyle:
                    arr[1,0].plot(kpt_array, self.correction[0,:,self.band_numbers[1]-1,temp_index], color = self.color[ifile],label=self.labels[ifile],linestyle=self.linestyle[ifile])
                    arr[2,0].plot(kpt_array, self.correction[0,:,self.band_numbers[0]-1,temp_index], color = self.color[ifile],label=self.labels[ifile],linestyle=self.linestyle[ifile])
                    arr[3,0].plot(kpt_array, self.correction[0,:,self.band_numbers[1]-1,temp_index]-self.correction[0,:,self.band_numbers[0]-1,temp_index], color = self.color[ifile],label=self.labels[ifile],linestyle=self.linestyle[ifile])
                else:
                    arr[1,0].plot(kpt_array, self.correction[0,:,self.band_numbers[1]-1,temp_index], color = self.color[ifile],label=self.labels[ifile])
                    arr[2,0].plot(kpt_array, self.correction[0,:,self.band_numbers[0]-1,temp_index], color = self.color[ifile],label=self.labels[ifile])
                    arr[3,0].plot(kpt_array, self.correction[0,:,self.band_numbers[1]-1,temp_index]-self.correction[0,:,self.band_numbers[0]-1,temp_index], color = self.color[ifile],label=self.labels[ifile])



                if ifile==0:

                    self.ylimits.append(arr[0,0].get_ylim())
                    for i in range(3):
                        self.ylimits.append(list(arr[i+1,0].get_ylim()))
                    
                else:
                    for i in range(3):
                        self.adjust_ylimits(arr[i+1,0], i+1, temp_index)                 
               

        # Set xlims, xticks, xtickslabels
        self.set_xaxis(arr[2][0],kpt_array)
    
        # Set visual references
        for i in range(4):   
            self.set_vrefs(arr[i][0], kpt_array, 0.) 


        # set legend
        self.set_legend_vbcb(arr[3,0])


        # Set y-axis
#            self.set_yaxis(arr[0,0], 'E-E(A-CBM) ({})'.format(self.units))
        self.set_yaxis(arr[0,0], 'E-Ef (unpert., {})'.format(self.units))
        self.set_yaxis(arr[1,0], 'CB corr ({})'.format(self.units))
        self.set_yaxis(arr[2,0], 'VB corr ({})'.format(self.units))
        self.set_yaxis(arr[3,0], 'Total corr ({})'.format(self.units))

        for i in range(4):
            self.set_ylimits_vbcb(arr[i,0],self.ylimits[i])


        for i in range(4):
            self.set_hrefs(self.ylimits[i], arr[i,0], gap_index,'black')
            self.set_hrefs(self.ylimits[i], arr[i,0], A_index,'black')

        self.set_main_title(fig)
        
        plt.show() 
        self.save_figure(fig)


    def plot_gap(self):

        file_qty = len(self.zpr_fnames)

        # Define figure
        fig, _arr = plt.subplots(1,1, figsize=self.figsize, squeeze=False)

        if self.split:
            fig2, _arr2 = plt.subplots(1,1, figsize=self.figsize, squeeze=False)

        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
  
            # Add indirect/direct gap information
            if self.indirect:
                if self.split:
                    self.gap_ren = self.zpr.indirect_gap_ren_split
                    self.gap_units = self.zpr.gap_energy_units

                else:
                    self.gap_ren = self.zpr.indirect_gap_ren
                    self.gap_units = self.zpr.gap_energy_units

            else:
                if self.split:
                    self.gap_ren = self.zpr.gap_renormalization_split
                    self.gap_units = self.zpr.gap_energy_units

                else:
                    self.gap_ren = self.zpr.gap_renormalization
                    self.gap_units = self.zpr.gap_energy_units


            if self.split:
                if self.expdata:
                    raise Exception("Experimental data option not implemented for split yet")

                else:
                    _arr[0][0].plot(self.temp, self.gap_ren[:,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                    _arr2[0][0].plot(self.temp, self.gap_ren[:,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])

            else:
                if self.expdata:
                    if self.gap_units == self.expdata.yaxis_units:
                        # they are both the same, no need to convert
                        if self.zero_gap_value:
                            if self.zero_gap_units == self.gap_units:
                                shift = self.zero_gap_value-self.gap_ren[0]
                            else:
                                if self.zero_gap_units == 'eV' and self.gap_units == 'meV':
                                    shift = self.zero_gap_value*1E3-self.gap_ren[0]
                                elif self.zero_gap_units == 'meV' and self.gap_units == 'eV':
                                    shift = self.zero_gap_value*1E-3-self.gap_ren[0]
                                else:
                                    raise Exception('Please provide zero gap value in eV or meV')

                        else:
                            shift = 0.

                        _arr[0][0].plot(self.temp, self.gap_ren+shift, marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                        _arr[0][0].plot(self.expdata.xaxis, self.expdata.yaxis, marker='s', linestyle='None',color='black')

                    else:
                        # the only case I write for now is meV to eV. Add Hartree to eV later??
                        if self.gap_units == 'meV' and self.expdata.yaxis_units == 'eV':

                            if self.expdata.xaxis_units != 'K':
                                raise Exception('Temperature axis must be in Kelvin')

                            if self.zero_gap_value:
                                if self.zero_gap_units == 'eV':
                                    shift = self.zero_gap_value-self.gap_ren[0]*1E-3
                                elif self.zero_gap_units == 'meV':
                                    shift = (self.zero_gap_value-self.gap_ren[0])*1E-3
                                else:
                                    raise Exception('Please provide zero gap value in eV or meV')
                            else:
                                shift = 0.

                            _arr[0][0].plot(self.temp, self.gap_ren*1E-3+shift, marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                            _arr[0][0].plot(self.expdata.xaxis, self.expdata.yaxis, marker='s', linestyle='None',color='black')

       
                else:
                    _arr[0][0].plot(self.temp, self.gap_ren, marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])

        self.set_xaxis(_arr[0][0], self.temp)
#        self.set_vrefs(_arr[0][0], self.temp, 0.)
        if self.expdata:
            self.set_yaxis(_arr[0][0], 'Gap energy ({})'.format(self.expdata.yaxis_units))
        else:
            self.set_yaxis(_arr[0][0], 'Gap renormalization ({})'.format(self.gap_units))
        self.set_legend_gap(_arr[0][0])
        self.set_main_title(fig) 

        if self.split:
            self.set_xaxis(_arr2[0][0], self.temp)
            self.set_vrefs(_arr2[0][0], self.temp, 0.)
            self.set_yaxis(_arr2[0][0], 'Gap renormalization ({})'.format(self.gap_units))
            self.set_legend_gap(_arr2[0][0])

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])

            self.set_main_title(fig2)
               

        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure(fig)

        plt.show()
            



    def plot_gap_separate(self):

        file_qty = len(self.zpr_fnames)

        # Define figure
        fig, _arr = plt.subplots(3,1, figsize=self.figsize, sharex=True, squeeze=False)

        if self.split:
            fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, sharex=True, squeeze=False)

        if self.verbose:
            h = open('output/compared_renormalization_{}_lowt.dat'.format(self.main_title),'w')


        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints

            if self.verbose:
                if ifile == 0:
                    h.write('Conduction band, Valence band and Gap Renormalization ({}) for {}\n'.format(self.zpr.gap_energy_units,self.main_title))
                    h.write('{:12s}'.format(''))
                    for t in range(self.ntemp):
                        h.write('   {:>5.0f}K  '.format(self.temp[t]))
                h.write('\n{:15s}'.format(self.labels[ifile]))
                
            
            # Add indirect/direct gap information
            if self.indirect:
                if self.split:
                    self.band_ren = self.zpr.indirect_gap_ren_band_split
                    self.gap_ren = self.zpr.indirect_gap_ren_split
                    self.gap_units = self.zpr.gap_energy_units

                else:
                    self.band_ren = self.zpr.indirect_gap_ren_band
                    self.gap_ren = self.zpr.indirect_gap_ren
                    self.gap_units = self.zpr.gap_energy_units

            else:
                if self.split:
                    self.band_ren = self.zpr.band_renormalization_split
                    self.gap_ren = self.zpr.gap_renormalization_split
                    self.gap_units = self.zpr.gap_energy_units

                else:
                    self.band_ren = self.zpr.band_renormalization
                    self.gap_ren = self.gap_renormalization
                    self.gap_units = self.zpr.gap_energy_units


            if self.split:
                #conduction band
                _arr[0][0].plot(self.temp, self.band_ren[:,0,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                _arr2[0][0].plot(self.temp, self.band_ren[:,1,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                #Valence band
                _arr[1][0].plot(self.temp, self.band_ren[:,0,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                _arr2[1][0].plot(self.temp, self.band_ren[:,1,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                #Gap correction
                _arr[2][0].plot(self.temp, self.gap_ren[:,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                _arr2[2][0].plot(self.temp, self.gap_ren[:,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])


            else:
                _arr[0][0].plot(self.temp, self.band_ren[:,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                _arr[1][0].plot(self.temp, self.band_ren[:,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                _arr[2][0].plot(self.temp, self.gap_ren, marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])

            if self.verbose:
                side = 0 # 0 = HA, 1 = AL
                if ifile == 1: #Reference file for error calculation 
                    err = np.zeros((self.ntemp,3))
                    refzpr = np.zeros((self.ntemp,3))
                    refzpr[:,0] = self.band_ren[:,side,1] # left/right, val/cond
                    refzpr[:,1] = self.band_ren[:,side,0]
                    refzpr[:,2] = self.gap_ren[:,side]

                if self.split:
                    if ifile > 1:
                        err[:,0] = (self.band_ren[:,side,1] - refzpr[:,0])/refzpr[:,0]*100
                        err[:,1] = (self.band_ren[:,side,0] - refzpr[:,1])/refzpr[:,1]*100
                        err[:,2] = (self.gap_ren[:,side] - refzpr[:,2])/refzpr[:,2]*100
                        
                    h.write('\n{:15s}'.format('  CB'))
                    for t in range(self.ntemp):
                        h.write('{:>8.4f}   '.format(self.band_ren[t,side,1]))
                    if ifile>1:
                        h.write('\n{:15s}'.format(''))
                        for t in range(self.ntemp):
                            h.write(' {:>6.2f}%   '.format(err[t,0]))
                    h.write('\n{:15s}'.format('  VB'))
                    for t in range(self.ntemp):
                        h.write('{:>8.4f}   '.format(self.band_ren[t,side,0]))
                    if ifile>1:
                        h.write('\n{:15s}'.format(''))
                        for t in range(self.ntemp):
                            h.write(' {:>6.2f}%   '.format(err[t,1]))
                    h.write('\n{:15s}'.format('  Gap'))
                    for t in range(self.ntemp):
                        h.write('{:>8.4f}   '.format(self.gap_ren[t,side]))
                    if ifile>1:
                        h.write('\n{:15s}'.format(''))
                        for t in range(self.ntemp):
                            h.write(' {:>6.2f}%   '.format(err[t,2]))

                else:
                    print('error calculation not implemented')
                    h.write('\n{:15s}'.format('  CB'))
                    for t in range(self.ntemp):
                        h.write('{:>8.4f}   '.format(self.band_ren[t,1]))
                    h.write('\n{:15s}'.format('  VB'))
                    for t in range(self.ntemp):
                        h.write('{:>8.4f}   '.format(self.band_ren[t,0]))
                    h.write('\n{:15s}'.format('  Gap'))
                    for t in range(self.ntemp):
                        h.write('{:>8.4f}   '.format(self.gap_ren[t]))



        for i in range(3):
            self.set_vrefs(_arr[i][0], self.temp, 0.)

            if self.split:
                self.set_vrefs(_arr2[i][0], self.temp, 0.)

        
        self.set_xaxis(_arr[2][0], self.temp)
        self.set_yaxis(_arr[0][0], 'CB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))


        if self.split:
            self.set_xaxis(_arr2[2][0], self.temp)
            self.set_yaxis(_arr2[0][0], 'CB ren ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[1][0], 'VB ren ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units))

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])

        else:
            self.set_title(_arr[0][0], self.main_title)

###

        self.set_legend_gap(_arr[2][0])
        self.set_main_title(fig) 

        if self.split:
            self.set_legend_gap(_arr2[2][0])

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])

            self.set_main_title(fig2)

        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure(fig)

        plt.show()
 
    def plot_pgap(self):

        file_qty = len(self.zpr_fnames)

        # Define figure
        fig, _arr = plt.subplots(2,1, figsize=self.figsize, squeeze=False, sharex=True)
        #plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig2, _arr2 = plt.subplots(2,1, figsize=self.figsize, squeeze=False, sharex=True)
          #  plt.subplots_adjust(hspace=0.05, top=0.95) 



        # Read explicit data, if required
        if self.gap_fname is not None:

            self.exgap = GAPfile(self.gap_fname, read=False)
            self.read_other_file()


#            if self.units == 'meV':
            if self.exgap.explicit_gap_energies_units == 'eV':
                self.explicit_gap_energies = self.exgap.explicit_gap_energies*1000
            else:
                self.explicit_gap_energies = self.exgap.explicit_gap_energies
#            else:
#                self.explicit_gap_energies = self.exgap.explicit_gap_energies

            self.explicit_pressures = self.exgap.explicit_pressures

#            if len(self.explicit_gap_energies) != file_qty:
#                raise Exception('gap_fname must have the same number of pressures as zpr data')

#            _arr[1][0].plot(self.explicit_pressures, self.explicit_gap_energies, 'k', marker = '*', label='unperturbed (explicit)') 

            

        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)
            print(zpr_file)
            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
            
            self.energy0 = self.zpr.unperturbed_gap_energy
            self.gridsize = self.zpr.gridsize
            self.gap_units = self.zpr.gap_energy_units

            if self.split :
                self.loc0 = self.zpr.unperturbed_gap_location_split
                self.gap_location = self.zpr.gap_location_split
                self.gap_energy = self.zpr.gap_energy_split
                self.gap_ren = self.zpr.gap_renormalization_split

                gap_index0 = np.zeros((2),dtype = int)
                gap_index0[0] = self.find_gap_index(self.loc0[0])
                gap_index0[1] = self.find_gap_index(self.loc0[1])

            else:
                self.loc0 = self.zpr.unperturbed_gap_location
                self.gap_location = self.zpr.gap_location
                self.gap_energy = self.zpr.gap_energy
                self.gap_ren = self.zpr.gap_renormalization
                gap_index0 = self.find_gap_index(self.loc0)

            if ifile==0:
                if self.split:
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp,2))
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((file_qty,2))

                else:
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp))
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((file_qty))
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

            print(self.gap_energy)
            if self.split:
                self.full_gap_energy[ifile,:,:] = self.gap_energy
                self.full_gap_ren[ifile,:,:] = self.gap_ren
                self.full_energy0[ifile,:] = self.energy0

            else:
                self.full_gap_energy[ifile,:] = self.gap_energy
                self.full_gap_ren[ifile,:] = self.gap_ren
                self.full_energy0[ifile] = self.energy0
            
            if self.gap_fname is not None:
                self.full_energy0[ifile] = self.explicit_gap_energies[ifile]
                self.full_gap_energy[ifile,:] =  self.full_gap_ren[ifile,:] + np.ones(self.ntemp)*self.full_energy0[ifile]


        if self.crit_pressure is not None:
            crit_index = self.find_temp_index()
        else:
            crit_index = None
        print(crit_index)


        if self.follow:
            self.full_energy0[-1,:] = -self.full_energy0[-1,:]
            self.full_gap_energy[-1,:,:] = -self.full_gap_energy[-1,:,:]

        for T in range(self.ntemp):
            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    # figure 1
                    if T==self.temp[0]:
                       
                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black', label='unperturbed')
                    else:
                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black')
    
                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #                self.set_legend_gap(_arr[0][0])
                    _arr[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T,0], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
            #        _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='interpolated')
    
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp+1, fontsize=16)
            #        self.set_legend_pgap(_arr2[0][0])
    #                _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=1.5, color=self.color[T], label=self.ref_temp[T])
                    _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='o', linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,T,0], marker='o', linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_energy0[s:,0],marker='d', color='black')
             #       _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='black')

                    # figure 2
                    if T==self.temp[0]:
                       
                        _arr2[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='unperturbed')
                    else:
                        _arr2[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black')
    
                    _arr2[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #                self.set_legend_gap(_arr2[0][0])
                    _arr2[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T,1], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
            #        _arr22[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='interpolated')

                    _arr2[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp+1, fontsize=16)
   
            #        self.set_legend_pgap(_arr22[0][0])
    #                _arr2[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=1.5, color=self.color[T], label=self.ref_temp[T])
                    _arr2[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='o', linewidth=1.5, color=self.color[T])
                    _arr2[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,T,1], marker='o', linewidth=1.5, color=self.color[T])
                    _arr2[1][0].plot(self.pressure[s:], self.full_energy0[s:,1],marker='d', color='black')
             #       _arr22[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='black')



                else:
                    if T==self.temp[0]:
                       
                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='unperturbed')
                    else:
                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black')
    
                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #                self.set_legend_gap(_arr[0][0])
                    _arr[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
            #        _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='interpolated')
    
                    self.set_legend_pgap(_arr[1][0])
            #        self.set_legend_pgap(_arr2[0][0])
    #                _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=1.5, color=self.color[T], label=self.ref_temp[T])
                    _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o', linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='black')
             #       _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='black')


            else:

                if self.split:
                    # figure 1
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr[1][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

                    #figure 2
                    _arr2[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr2[1][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

                else:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr[1][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

        for i in range(2):
            self.set_vrefs(_arr[i][0], self.pressure, 0.)
            self.set_hrefs(self.ylims, _arr[i][0], self.crit_pressure,'black')
            self.set_vrefs(_arr2[i][0], self.pressure, 0.)
            self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'black')



        self.set_xaxis(_arr[1][0], self.pressure)
        self.set_yaxis(_arr[0][0], 'Gap renormalization ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'Gap energy ({})'.format(self.gap_units))

        if self.split:
            self.set_xaxis(_arr2[1][0], self.pressure)
            self.set_yaxis(_arr2[0][0], 'Gap renormalization ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[1][0], 'Gap energy ({})'.format(self.gap_units))

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])
        else:
            self.set_title(_arr[0][0], self.title[0])


        plt.show()
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure(fig)


#
#        if self.gap_fname is not None:
#            fig2, _arr2 = plt.subplots(1,1, figsize=self.figsize, squeeze=False)
#
#            _arr2[0][0].plot(self.explicit_pressures, self.explicit_gap_energies, 'k', marker = 'd', label='explicit') 
#
#            if crit_index is not None:
#
#                s = crit_index+1
#                _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='blue', label='interpolated')
#                _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='blue')
#
#            else:
#                _arr2[0][0].plot(self.pressure, self.full_energy0,marker='d', color='blue')
#
#            self.set_legend_pgap(_arr2[0][0])
#            self.set_xaxis(_arr2[0][0], self.pressure)
#            self.set_yaxis(_arr2[0][0], 'Gap energy ({})'.format(self.gap_units))
#            self.set_main_title(fig2)
#
#            plt.show()
#            self.save_figure2(fig2)

    def plot_pgap_indirect(self):


        file_qty = len(self.zpr_fnames)

        # Define figure
        only=True

        if only:
            fig, _arr = plt.subplots(1,1, figsize=self.figsize, squeeze=False, sharex=True)
        else:
            fig, _arr = plt.subplots(2,1, figsize=self.figsize, squeeze=False, sharex=True)
        #plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            if only:
                fig2, _arr2 = plt.subplots(1,1, figsize=self.figsize, squeeze=False, sharex=True)
            else:
                fig2, _arr2 = plt.subplots(2,1, figsize=self.figsize, squeeze=False, sharex=True)
          #  plt.subplots_adjust(hspace=0.05, top=0.95) 



        # Read explicit data, if required
        if self.gap_fname is not None:

            self.exgap = GAPfile(self.gap_fname, read=False)
            self.read_other_file()


            if self.exgap.explicit_gap_energies_units == 'eV':
                self.explicit_gap_energies = self.exgap.explicit_gap_energies*1000
            else:
                self.explicit_gap_energies = self.exgap.explicit_gap_energies

#            else:
#                self.explicit_gap_energies = self.exgap.explicit_gap_energies

            self.explicit_pressures = self.exgap.explicit_pressures

#            if len(self.explicit_gap_energies) != file_qty:
#                raise Exception('gap_fname must have the same number of pressures as zpr data')

#            _arr[1][0].plot(self.explicit_pressures, self.explicit_gap_energies, 'k', marker = '*', label='unperturbed (explicit)') 

            

        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)
            print(zpr_file)
            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
            
            #self.gridsize = self.zpr.gridsize
            self.gap_units = self.zpr.gap_energy_units

            if self.split :
                self.energy0 = self.zpr.unperturbed_indirect_gap_energy_split
                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
                self.gap_location = self.zpr.indirect_gap_location_split
                self.gap_energy = self.zpr.indirect_gap_energy_split
                self.gap_ren = self.zpr.indirect_gap_ren_split
                self.band_energy = self.zpr.indirect_gap_energy_band_split

                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                for a in range(2):
                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

            elif self.split2:
                if self.pressure[ifile] > self.crit_pressure:
                    side = 1    
                else:
                    side = 0     

                self.energy0 = self.zpr.unperturbed_indirect_gap_energy_split
                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
                self.gap_location = self.zpr.indirect_gap_location_split
                self.gap_energy = self.zpr.indirect_gap_energy_split
                self.gap_ren = self.zpr.indirect_gap_ren_split
                self.band_energy = self.zpr.indirect_gap_energy_band_split

                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                for a in range(2):
                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

 

            else:
                self.energy0 = self.zpr.unperturbed_indirect_gap_energy
                self.loc0 = self.zpr.unperturbed_indirect_gap_location
                self.gap_location = self.zpr.indirect_gap_location
                self.gap_energy = self.zpr.indirect_gap_energy
                self.gap_ren = self.zpr.indirect_gap_ren

                gap_index0 = np.zeros((2),dtype = int)
                gap_index0[0] = self.find_gap_index(self.loc0[0])
                gap_index0[1] = self.find_gap_index(self.loc0[1])

            # Initialize full arrays
            if ifile==0:
                if self.split:
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp,2))
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((len(self.explicit_pressures),2))

                if self.split2: # single array, left gap for trovoal phase and right gap for topol phase
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp))
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((file_qty))

                else:
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp))
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((file_qty))
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

#            print(self.gap_energy)
            if self.split:
                self.full_gap_energy[ifile,:,:] = self.gap_energy
                self.full_gap_ren[ifile,:,:] = self.gap_ren
                if not self.gap_fname:
                    self.full_energy0[ifile,:] = self.energy0
                else:
                    # find P in explicit pressures array
                    p_index = self.find_pressure_index(self.pressure[ifile])
                    self.full_energy0[ifile,:] = self.explicit_gap_energies[p_index]

            if self.split2:
                self.full_gap_ren[ifile,:] = self.gap_ren[:,side]
                self.full_gap_energy[ifile,:] = self.gap_energy[:,side]
                print(self.gap_fname)
                if not self.gap_fname:
                    self.full_energy0[ifile] = self.energy0[ifile,side]
                else:
                    # find P in explicit pressures array
                    p_index = self.find_pressure_index(self.pressure[ifile])
                    self.full_energy0[ifile] = self.explicit_gap_energies[p_index]
            else:
                self.full_gap_energy[ifile,:] = self.gap_energy
                self.full_gap_ren[ifile,:] = self.gap_ren
                if not self.gap_fname:
                    self.full_energy0[ifile] = self.energy0
                else:
                    p_index = self.find_pressure_index(self.pressure[ifile])
                    self.full_energy0[ifile] = self.explicit_energies[p_index]

            if self.gap_fname is not None:
                if self.split:
                    for a in range(2):
                #        self.full_energy0[ifile,a] = self.explicit_gap_energies[ifile]
                        self.full_gap_energy[ifile,:,a] =  self.full_gap_ren[ifile,:,a] + np.ones(self.ntemp)*self.full_energy0[ifile,a]
                elif self.split2:
                    self.full_gap_energy[ifile,:] = self.full_gap_ren[ifile,:] + np.ones(self.ntemp)*self.full_energy0[ifile]
                else:
                 #   self.full_energy0[ifile] = self.explicit_gap_energies[ifile]
                    self.full_gap_energy[ifile,:] =  self.full_gap_ren[ifile,:] + np.ones(self.ntemp)*self.full_energy0[ifile]


        if self.crit_pressure is not None:
            crit_index = self.find_temp_index()
            if self.gap_fname is not None:
                crit_index2 = self.find_temp_index2(self.explicit_pressures)
                print('crit_index2',crit_index2)
        else:
            crit_index = None


        if self.follow:
            self.full_energy0[-1,:] = -self.full_energy0[-1,:]
            self.full_gap_energy[-1,:,:] = -self.full_gap_energy[-1,:,:]

    # Extrapolate behavior linearly towards TPT... it IS a bit sketchy.
        pressure1 = np.arange(1.5,2.0,0.1)
        tmparr = np.arange(2.1,2.6,0.1) # change 2.2 to 2.6 foir full BZ
        self.extr_pressure2 = np.arange(2.3, 3.0, 0.1)

        # linearize unperturbed gap energy to extrapolate inside WSM phase
        y0,y1 = np.polyfit(self.explicit_pressures[15:20], self.explicit_gap_energies[15:20],1)
        tmp2 = y1 + tmparr*y0
        self.extr_pressure1 = np.concatenate((pressure1, tmparr))


        self.extr_full_gap_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp))
        self.extr_full_gap_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp))
        self.extr_full_gap_energy1 = np.zeros((len(self.extr_pressure1), self.ntemp))
        self.extr_full_gap_energy2 = np.zeros((len(self.extr_pressure2), self.ntemp))
        extr_energy01 = np.zeros((len(pressure1)))
        self.extr_full_energy02 = np.zeros((len(self.extr_pressure2)))

        # Get full_gap_energy0 from explicit data
        for p, pres in enumerate(pressure1):
            ind = self.find_pressure_index(pres)
            extr_energy01[p] = self.explicit_gap_energies[ind]

        self.extr_full_energy01 = np.concatenate((extr_energy01,tmp2))

        for p, pres in enumerate(self.extr_pressure2):
            ind = self.find_pressure_index(pres)
            self.extr_full_energy02[p] = self.explicit_gap_energies[ind]

        for T in range(self.ntemp):

            if self.split:
                #trivial side
                x0,x1 = np.polyfit(self.pressure[0:crit_index+1], self.full_gap_ren[0:crit_index+1,T,0], 1)
                self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]

                # topol side
                x0,x1 =  np.polyfit(self.pressure[crit_index+1:], self.full_gap_ren[crit_index+1:,T,0], 1)
                self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]

            elif self.split2:
                #trivial side
                x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_gap_ren[crit_index-1:crit_index+1,T], 1)
                self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]

                # topol side
                x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_gap_ren[crit_index+1:crit_index+3,T], 1)
                self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]



#        print(self.extr_pressure1)
#        print(self.extr_full_gap_energy1)
#        print(self.extr_pressure2)
#        print(self.extr_full_gap_energy2)

        for T in range(self.ntemp):
            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    # figure 1
                    if T==self.temp[0]:
                        ###### FIX ME : keep explicit_pressures and explicit_gap_energies, just get a new crit index for this array! 
                        if self.gap_fname is not None:
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='unperturbed')
                            if self.crit_index2 is not None:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='unperturbed')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='unperturbed')

                            else:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='unperturbed')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='unperturbed')

                        else:
                            if only:
                                _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='unperturbed')
                            else:
                                _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='unperturbed')
                    else:
                        if only:
                            _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black')
                        else:
                            _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black')
    
                    if only:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T,0], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    else:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
        #                self.set_legend_gap(_arr[0][0])
                        _arr[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T,0], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                #        _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='o', color='black', label='interpolated')
    
                    # extrapolated data
                    if only:
                        _arr[0,0].plot(self.extr_pressure1, self.extr_full_gap_energy1[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                    else:
                        _arr[0,0].plot(self.extr_pressure1[:-2], self.extr_full_gap_ren1[:-2,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure1, self.extr_full_gap_energy1[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


    ##-------------------------------------- passed TPT
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.30), ncol = (self.ntemp+1), fontsize=12)
                    if only:
                        handles, labels = _arr[0,0].get_legend_handles_labels()
                        _arr[0][0].legend(self.flip(handles, 4), self.flip(labels, 4), loc='lower center', bbox_to_anchor=(0.5,-0.50), ncol=4, fontsize=20, numpoints=1)
                    else:
                        handles, labels = _arr[1,0].get_legend_handles_labels()
                        _arr[1][0].legend(self.flip(handles, 4), self.flip(labels, 4), loc='lower center', bbox_to_anchor=(0.5,-0.50), ncol=4, fontsize=20, numpoints=1)
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.30), ncol = (self.ntemp+1), fontsize=22)

#                    _arr[1][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.27,1.62), ncol = 1, fontsize=22)

            #        self.set_legend_pgap(_arr2[0][0])
    #                _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=2.0, color=self.color[T], label=self.ref_temp[T])
                    if T==self.temp[0]:
                        if self.gap_fname is not None:
#
                            if crit_index2 is not None:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[crit_index2+2,:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label='unperturbed')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[crit_index2+2:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label='unperturbed')

                    if only:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T,0], marker='o', linewidth=2.0, color=self.color[T])
                    else:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='o', linewidth=2.0, color=self.color[T])
                        _arr[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,T,0], marker='o', linewidth=2.0, color=self.color[T])

                    if only:
                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                    else:
                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


                    if self.gap_fname is not None:
                        x=1
                        #_arr[1][0].plot(self.explicit_pressures[s:], self.explicit_gap_energies[s:],marker='d', color='black')
                    else:
                        if only:
                            _arr[0][0].plot(self.pressure[s:], self.full_energy0[s:,0],marker='o', color='black')
                        else:
                            _arr[1][0].plot(self.pressure[s:], self.full_energy0[s:,0],marker='o', color='black')
             #       _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='o', color='black')

                    # figure 2
                    if T==self.temp[0]:
                        if self.gap_fname is not None:
                            if only:
                                _arr2[0][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='unperturbed')
                            else:
                                arr2[1][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='unperturbed')
                        else:
                            if only:
                                _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='unperturbed')
                            else:
                                _arr2[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='unperturbed')
                    else:
                        if only:
                            _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black')
                        else:
                            _arr2[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black')
    
                    if only:
                        _arr2[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T,1], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    else:
                        _arr2[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
        #                self.set_legend_gap(_arr2[0][0])
                        _arr2[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T,1], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                #        _arr22[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='interpolated')

                    if only:
                        _arr2[0][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp+1, fontsize=14)
                    else:
                        _arr2[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp+1, fontsize=14)
   
            #        self.set_legend_pgap(_arr22[0][0])
    #                _arr2[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=2.0, color=self.color[T], label=self.ref_temp[T])
                    if only:
                        _arr2[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T,1], marker='o', linewidth=2.0, color=self.color[T])
                    else:
                        _arr2[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='o', linewidth=2.0, color=self.color[T])
                        _arr2[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,T,1], marker='o', linewidth=2.0, color=self.color[T])

                    if self.gap_fname is not None:
                        if only:
                            _arr2[0][0].plot(self.explicit_pressures[s:], self.full_energy0[s:,1],marker='d', color='black')
                        else:
                            _arr2[1][0].plot(self.explicit_pressures[s:], self.full_energy0[s:,1],marker='d', color='black')
                    else:
                        if only:
                            _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:,1],marker='d', color='black')
                        else:
                            _arr2[1][0].plot(self.pressure[s:], self.full_energy0[s:,1],marker='d', color='black')
             #       _arr22[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='black')


       #-------- crit_index but no split or split2
                else:
                    if T==self.temp[0]:
                        ###### FIX ME : keep explicit_pressures and explicit_gap_energies, just get a new crit index for this array! 
                        if self.gap_fname is not None:
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='unperturbed')
                            if crit_index2 is not None:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='unperturbed')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='unperturbed')
                            else:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='unperturbed')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='unperturbed')

                        else:
                            if only:
                                _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='unperturbed')
                            else:
                                _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='unperturbed')
#                    else:
#                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black')
    
                    if only:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    else:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
        #                self.set_legend_gap(_arr[0][0])
                        _arr[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                #        _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='o', color='black', label='interpolated')

    #                if only:
    #                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #                else:
    #                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #    #                self.set_legend_gap(_arr[0][0])
    #                    _arr[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #            #        _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='o', color='black', label='interpolated')
    
                    # extrapolated data
                    if only:
                        _arr[0,0].plot(self.extr_pressure1, self.extr_full_gap_energy1[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                    else:
                        _arr[0,0].plot(self.extr_pressure1[:-2], self.extr_full_gap_ren1[:-2,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure1, self.extr_full_gap_energy1[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


    ##-------------------------------------- passed TPT
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.30), ncol = (self.ntemp+1), fontsize=12)
                    if only:
                        handles, labels = _arr[0,0].get_legend_handles_labels()
                        _arr[0][0].legend(self.flip(handles, 4), self.flip(labels, 4), loc='lower center', bbox_to_anchor=(0.5,-0.50), ncol=4, fontsize=20, numpoints=1)
                    else:
                        handles, labels = _arr[1,0].get_legend_handles_labels()
                        _arr[1][0].legend(self.flip(handles, 4), self.flip(labels, 4), loc='lower center', bbox_to_anchor=(0.5,-0.50), ncol=4, fontsize=20, numpoints=1)
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.30), ncol = (self.ntemp+1), fontsize=22)

#                    _arr[1][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.27,1.62), ncol = 1, fontsize=22)

            #        self.set_legend_pgap(_arr2[0][0])
    #                _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=2.0, color=self.color[T], label=self.ref_temp[T])
                    if T==self.temp[0]:
                        if self.gap_fname is not None:
#
                            if crit_index2 is not None:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[crit_index2+2:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label=None)
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[crit_index2+2:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label=None)

                    if only:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o', linewidth=2.0, color=self.color[T])
                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])

                    else:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=2.0, color=self.color[T])
                        _arr[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,t], marker='o', linewidth=2.0, color=self.color[t])

                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


                    if self.gap_fname is not None:
                        x=1
                        #_arr[1][0].plot(self.explicit_pressures[s:], self.explicit_gap_energies[s:],marker='d', color='black')
                    else:
                        if only:
                            _arr[0][0].plot(self.pressure[s:], self.full_energy0[s:,0],marker='o', color='black')
                        else:
                            _arr[1][0].plot(self.pressure[s:], self.full_energy0[s:,0],marker='o', color='black')
             #       _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='o', color='black')


            #---------------------- No crit_index
            else:

                if only:
                    raise Exception('Only keyword not implemented yet')
                if self.split:
                    # figure 1
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr[1][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

                    #figure 2
                    _arr2[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr2[1][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

                else:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr[1][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

        if only:
            self.set_vrefs(_arr[0][0], self.pressure, 0.)
            self.set_hrefs(self.ylims, _arr[0][0], self.crit_pressure,'black')
            if self.split:
                self.set_vrefs(_arr2[0][0], self.pressure, 0.)
                self.set_hrefs(self.ylims, _arr2[0][0], self.crit_pressure,'black')

        else:
            for i in range(2):
        #            if i==0 :
                    self.set_vrefs(_arr[i][0], self.pressure, 0.)
                    self.set_hrefs(self.ylims, _arr[i][0], self.crit_pressure,'black')
                    if self.split:
                        self.set_vrefs(_arr2[i][0], self.pressure, 0.)
                        self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'black')

        limms = [[-80.,60.],[0.,350]]
#        limms = [[-60.,60.],[0.,80.],[-80.,40.]]

        if only:
            ylims = limms[1]    
            self.set_vrefs(_arr[0][0], self.pressure, 0.)
#            if i==0:
            self.set_hrefs(ylims, _arr[0][0], self.crit_pressure,'gray')
            self.set_hrefs(ylims, _arr[0][0], self.crit_pressure+0.2,'gray')
            _arr[0][0].fill([2.08,2.28,2.28,2.08],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
#            else:
            self.set_hrefs(ylims, _arr[0,0], 2.18,'orange')
            self.set_hrefs(ylims, _arr[0,0], 2.50,'orange')
            _arr[0,0].fill([2.18,2.50,2.50,2.18],[ylims[0],ylims[0],ylims[1],ylims[1]],'orange',alpha=0.2)
            self.set_hrefs(ylims, _arr[0,0], 2.68,'magenta')


        else:
            for i in range(2):
    #            lims = _arr[i,0].get_ylim() 
    #            ylims = lims
                ylims = limms[i]    
                self.set_vrefs(_arr[i][0], self.pressure, 0.)
    #            if i==0:
                self.set_hrefs(ylims, _arr[i][0], self.crit_pressure,'black')
                self.set_hrefs(ylims, _arr[i][0], self.crit_pressure+0.2,'black')
                _arr[i][0].fill([2.08,2.28,2.28,2.08],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
    #            else:
                self.set_hrefs(ylims, _arr[i,0], 2.28,'black')
                self.set_hrefs(ylims, _arr[i,0], 2.65,'black')
                _arr[i,0].fill([2.28,2.65,2.65,2.28],[ylims[0],ylims[0],ylims[1],ylims[1]],'orange',alpha=0.2)
 

        if only:
            self.set_xaxis(_arr[0][0], self.pressure)
            self.set_yaxis(_arr[0][0], 'Gap energy ({})'.format(self.gap_units))
            _arr[0,0].set_ylim(limms[1])

        else:
            self.set_xaxis(_arr[1][0], self.pressure)
            self.set_yaxis(_arr[0][0], 'Gap renormalization ({})'.format(self.gap_units))
            self.set_yaxis(_arr[1][0], 'Gap energy ({})'.format(self.gap_units))

            _arr[1,0].set_ylim(limms[1])


        if only:
            _arr[0,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=24) 
            _arr[0,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=24) 

#            _arr[0,0].text(self.explicit_pressures[5],290, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=26)
#            _arr[0,0].text(self.explicit_pressures[38],290, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=26)


            _arr[0,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=26, weight='bold', color='gray')
            _arr[0,0].text(self.explicit_pressures[25],250, r'WSM 300K', fontsize=26, weight='bold', color='orange')

        else:
            _arr[1,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=24) 
            _arr[1,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=24) 

#            _arr[1,0].text(self.explicit_pressures[5],290, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=26)
#            _arr[1,0].text(self.explicit_pressures[38],290, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=26)


            _arr[1,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=26, weight='bold', color='gray')
            _arr[1,0].text(self.explicit_pressures[25],250, r'WSM 300K', fontsize=26, weight='bold', color='orange')



        if self.split:
            if only:
                self.set_xaxis(_arr2[0][0], self.pressure)
                self.set_yaxis(_arr2[0][0], 'Gap energy ({})'.format(self.gap_units))
                self.set_title(_arr[0][0], self.title[0])

            else:
                self.set_xaxis(_arr2[1][0], self.pressure)
                self.set_yaxis(_arr2[0][0], 'Gap renormalization ({})'.format(self.gap_units))
                self.set_yaxis(_arr2[1][0], 'Gap energy ({})'.format(self.gap_units))

                self.set_title(_arr[0][0], self.title[0])
                self.set_title(_arr2[0][0], self.title[1])
        else:
            self.set_title(_arr[0][0], self.main_title)

#        fig.subplots_adjust(hspace=0.09, bottom=0.18)

#        plt.subplots_adjust(left = 0.1, bottom = 0.11, right=0.9, top = 0.93, wspace=0.2,hspace=0.07)
        if only:
            fig.subplots_adjust(bottom= 0.31, right = 0.94)
        else:
            fig.subplots_adjust(left = 0.1, bottom = 0.21, right=0.95, top = 0.95, wspace=0.2,hspace=0.07)
        plt.show()
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure(fig)


    def plot_pgap_indirect_only(self):


        file_qty = len(self.zpr_fnames)

        # Define figure
        fig, _arr = plt.subplots(1,1, figsize=self.figsize, squeeze=False, sharex=True)
        #plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig2, _arr2 = plt.subplots(1,1, figsize=self.figsize, squeeze=False, sharex=True)
          #  plt.subplots_adjust(hspace=0.05, top=0.95) 



        # Read explicit data, if required
        if self.gap_fname is not None:

            self.exgap = GAPfile(self.gap_fname, read=False)
            self.read_other_file()


            if self.exgap.explicit_gap_energies_units == 'eV':
                self.explicit_gap_energies = self.exgap.explicit_gap_energies*1000
            else:
                self.explicit_gap_energies = self.exgap.explicit_gap_energies

#            else:
#                self.explicit_gap_energies = self.exgap.explicit_gap_energies

            self.explicit_pressures = self.exgap.explicit_pressures

#            if len(self.explicit_gap_energies) != file_qty:
#                raise Exception('gap_fname must have the same number of pressures as zpr data')

#            _arr[1][0].plot(self.explicit_pressures, self.explicit_gap_energies, 'k', marker = '*', label='unperturbed (explicit)') 

            

        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)
            print(zpr_file)
            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
            
            #self.gridsize = self.zpr.gridsize
            self.gap_units = self.zpr.gap_energy_units

            if self.split :
                self.energy0 = self.zpr.unperturbed_indirect_gap_energy_split
                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
                self.gap_location = self.zpr.indirect_gap_location_split
                self.gap_energy = self.zpr.indirect_gap_energy_split
                self.gap_ren = self.zpr.indirect_gap_ren_split
                self.band_energy = self.zpr.indirect_gap_energy_band_split

                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                for a in range(2):
                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])
 

            else:
                self.energy0 = self.zpr.unperturbed_indirect_gap_energy
                self.loc0 = self.zpr.unperturbed_indirect_gap_location
                self.gap_location = self.zpr.indirect_gap_location
                self.gap_energy = self.zpr.indirect_gap_energy
                self.gap_ren = self.zpr.indirect_gap_ren

                gap_index0 = np.zeros((2),dtype = int)
                gap_index0[0] = self.find_gap_index(self.loc0[0])
                gap_index0[1] = self.find_gap_index(self.loc0[1])

            if ifile==0:
                if self.split:
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp,2))
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((len(self.explicit_pressures),2))

                else:
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp))
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((file_qty))
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

#            print(self.gap_energy)
            if self.split:
                self.full_gap_energy[ifile,:,:] = self.gap_energy
                self.full_gap_ren[ifile,:,:] = self.gap_ren
                if not self.gap_fname:
                    self.full_energy0[ifile,:] = self.energy0
                else:
                    # find P in explicit pressures array
                    p_index = self.find_pressure_index(self.pressure[ifile])
                    self.full_energy0[ifile,:] = self.explicit_gap_energies[p_index]
            else:
                self.full_gap_energy[ifile,:] = self.gap_energy
                self.full_gap_ren[ifile,:] = self.gap_ren
                if not self.gap_fname:
                    self.full_energy0[ifile] = self.energy0
                else:
                    p_index = self.find_pressure_index(self.pressure[ifile])
                    self.full_energy0[ifile] = self.explicit_energies[p_index]
            
            if self.gap_fname is not None:
                if self.split:
                    for a in range(2):
                #        self.full_energy0[ifile,a] = self.explicit_gap_energies[ifile]
                        self.full_gap_energy[ifile,:,a] =  self.full_gap_ren[ifile,:,a] + np.ones(self.ntemp)*self.full_energy0[ifile,a]
                else:
                 #   self.full_energy0[ifile] = self.explicit_gap_energies[ifile]
                    self.full_gap_energy[ifile,:] =  self.full_gap_ren[ifile,:] + np.ones(self.ntemp)*self.full_energy0[ifile]


        if self.crit_pressure is not None:
            crit_index = self.find_temp_index()
            if self.gap_fname is not None:
                crit_index2 = self.find_temp_index2(self.explicit_pressures)
        else:
            crit_index = None


        if self.follow:
            self.full_energy0[-1,:] = -self.full_energy0[-1,:]
            self.full_gap_energy[-1,:,:] = -self.full_gap_energy[-1,:,:]

    # Extrapolate behavior linearly towards TPT... it IS a bit sketchy.
        pressure1 = np.arange(1.5,2.0,0.1)
        tmparr = np.arange(2.1,2.6,0.1) # change 2.2 to 2.6 foir full BZ
        self.extr_pressure2 = np.arange(2.3, 3.6, 0.1)

        # linearize unperturbed gap energy to extrapolate inside WSM phase
        y0,y1 = np.polyfit(self.explicit_pressures[15:20], self.explicit_gap_energies[15:20],1)
        tmp2 = y1 + tmparr*y0
        self.extr_pressure1 = np.concatenate((pressure1, tmparr))


        self.extr_full_gap_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp))
        self.extr_full_gap_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp))
        self.extr_full_gap_energy1 = np.zeros((len(self.extr_pressure1), self.ntemp))
        self.extr_full_gap_energy2 = np.zeros((len(self.extr_pressure2), self.ntemp))
        extr_energy01 = np.zeros((len(pressure1)))
        self.extr_full_energy02 = np.zeros((len(self.extr_pressure2)))

        # Get full_gap_energy0 from explicit data
        for p, pres in enumerate(pressure1):
            ind = self.find_pressure_index(pres)
            extr_energy01[p] = self.explicit_gap_energies[ind]

        self.extr_full_energy01 = np.concatenate((extr_energy01,tmp2))

        for p, pres in enumerate(self.extr_pressure2):
            ind = self.find_pressure_index(pres)
            self.extr_full_energy02[p] = self.explicit_gap_energies[ind]

        for T in range(self.ntemp):

            #trivial side
            x0,x1 = np.polyfit(self.pressure[0:crit_index+1], self.full_gap_ren[0:crit_index+1,T,0], 1)
            self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
            self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]

            # topol side
            x0,x1 =  np.polyfit(self.pressure[crit_index+1:], self.full_gap_ren[crit_index+1:,T,0], 1)
            self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
            self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]

#        print(self.extr_pressure1)
#        print(self.extr_full_gap_energy1)
#        print(self.extr_pressure2)
#        print(self.extr_full_gap_energy2)

        for T in range(self.ntemp):
            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    # figure 1
                    if T==self.temp[0]:
                        ###### FIX ME : keep explicit_pressures and explicit_gap_energies, just get a new crit index for this array! 
                        if self.gap_fname is not None:
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='unperturbed')
                            _arr[0][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='unperturbed')


                        else:
                            _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='unperturbed')
                    else:
                        _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black')
    
#                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #                self.set_legend_gap(_arr[0][0])
                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T,0], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
            #        _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='o', color='black', label='interpolated')
    
                    # extrapolated data
 #                   _arr[0,0].plot(self.extr_pressure1, self.extr_full_gap_ren1[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                    _arr[0,0].plot(self.extr_pressure1, self.extr_full_gap_energy1[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


    ##-------------------------------------- passed TPT
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.30), ncol = (self.ntemp+1), fontsize=12)
                    handles, labels = _arr[0,0].get_legend_handles_labels()
                    _arr[0][0].legend(self.flip(handles, 4), self.flip(labels, 4), loc='lower center', bbox_to_anchor=(0.5,-0.3), ncol=4, fontsize=20, numpoints=1)
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.30), ncol = (self.ntemp+1), fontsize=22)

#                    _arr[1][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.27,1.62), ncol = 1, fontsize=22)

            #        self.set_legend_pgap(_arr2[0][0])
    #                _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=2.0, color=self.color[T], label=self.ref_temp[T])
#                    _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='o', linewidth=2.0, color=self.color[T])
                    if self.gap_fname is not None:
                        _arr[0][0].plot(self.explicit_pressures[crit_index2+2:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black')

                    _arr[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T,0], marker='o', linewidth=2.0, color=self.color[T])

 #                   _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                    _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


                    if self.gap_fname is not None:
                        x=1
                        #_arr[1][0].plot(self.explicit_pressures[s:], self.explicit_gap_energies[s:],marker='d', color='black')
                    else:
                        _arr[0][0].plot(self.pressure[s:], self.full_energy0[s:,0],marker='o', color='black')
             #       _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='o', color='black')

                    # figure 2
                    if T==self.temp[0]:
                        if self.gap_fname is not None:
                            _arr2[0][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='unperturbed')
                        else:
                            _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='unperturbed')
                    else:
                        _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black')
    
#                    _arr2[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #                self.set_legend_gap(_arr2[0][0])
                    _arr2[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T,1], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
            #        _arr22[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='interpolated')

                    _arr2[0][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp+1, fontsize=16)
   
            #        self.set_legend_pgap(_arr22[0][0])
    #                _arr2[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=2.0, color=self.color[T], label=self.ref_temp[T])
#                    _arr2[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='o', linewidth=2.0, color=self.color[T])
                    _arr2[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T,1], marker='o', linewidth=2.0, color=self.color[T])

                    if self.gap_fname is not None:
                        _arr2[0][0].plot(self.explicit_pressures[s:], self.full_energy0[s:,1],marker='d', color='black')
                    else:
                        _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:,1],marker='d', color='black')
             #       _arr22[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='black')



                else:
                    if T==self.temp[0]:
                        if self.gap_fname is not None:
                            _arr[0][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s],marker='d', color='black', label='unperturbed')
                        else:
                            _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='unperturbed')
                    else:
                        _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black')
    
#                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
    #                self.set_legend_gap(_arr[0][0])
                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
            #        _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='interpolated')
    
#                    self.set_legend_pgap(_arr[1][0])
            #        self.set_legend_pgap(_arr2[0][0])
    #                _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=2.0, color=self.color[T], label=self.ref_temp[T])
 #                   _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', linewidth=2.0, color=self.color[T])
                    _arr[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o', linewidth=2.0, color=self.color[T])
                    if self.gap_fname is not None:
                        _arr[0][0].plot(self.explicit_pressures[s:], self.full_energy0[s:],marker='d', color='black')
                    else:
                        _arr[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='black')
             #       _arr2[0][0].plot(self.pressure[s:], self.full_energy0[s:],marker='d', color='black')


            #---------------------- No crit_index
            else:
                if self.split:
                    # figure 1
#                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[0][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr[0][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

                    #figure 2
 #                   _arr2[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[0][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr2[0][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

                else:
  #                  _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='o', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[0][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    #_arr[1][0].plot(self.pressure, self.full_energy0,marker='d', color='black')

        for i in range(1):
#            if i==0 :
            self.set_vrefs(_arr[i][0], self.pressure, 0.)
            self.set_hrefs(self.ylims, _arr[i][0], self.crit_pressure,'black')
            self.set_vrefs(_arr2[i][0], self.pressure, 0.)
            self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'black')

        limms = [[0.,350]]
#        limms = [[-60.,60.],[0.,80.],[-80.,40.]]
        for i in range(1):
#            lims = _arr[i,0].get_ylim() 
#            ylims = lims
            ylims = limms[i]    
            self.set_vrefs(_arr[i][0], self.pressure, 0.)
#            if i==0:
        # For 0K WSM
        #    self.set_hrefs(ylims, _arr[i][0], self.crit_pressure,'black')
        #    self.set_hrefs(ylims, _arr[i][0], self.crit_pressure+0.2,'black')
        #    _arr[i][0].fill([2.08,2.28,2.28,2.08],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
#       # For 300K WSM
            self.set_hrefs(ylims, _arr[i,0], 2.28,'black')
            self.set_hrefs(ylims, _arr[i,0], 2.65,'black')
            _arr[i,0].fill([2.08,2.28,2.28,2.08],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
            _arr[i,0].fill([2.28,2.65,2.65,2.28],[ylims[0],ylims[0],ylims[1],ylims[1]],'orange',alpha=0.2)


        self.set_xaxis(_arr[0][0], self.pressure)
#        self.set_yaxis(_arr[0][0], 'Gap renormalization ({})'.format(self.gap_units))
        self.set_yaxis(_arr[0][0], 'Gap energy ({})'.format(self.gap_units))

        _arr[0,0].set_ylim(limms[0])
#        plt.subplots_adjust(left = 0.1, bottom = 0.11, right=0.9, top = 0.93, wspace=0.2,hspace=0.07)
        fig.subplots_adjust(left = 0.1, bottom = 0.21, right=0.95, top = 0.93, wspace=0.2,hspace=0.07)


        _arr[0,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=28) 
        _arr[0,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=28) 

        _arr[0,0].text(self.explicit_pressures[8],300, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=30)
        _arr[0,0].text(self.explicit_pressures[35],300, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=30)

        _arr[0,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=28, weight='bold', color='gray')
        _arr[0,0].text(self.explicit_pressures[26],250, r'WSM 300K', fontsize=28, weight='bold', color='orange')



        if self.split:
            self.set_xaxis(_arr2[0][0], self.pressure)
#            self.set_yaxis(_arr2[0][0], 'Gap renormalization ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[0][0], 'Gap energy ({})'.format(self.gap_units))

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])
        else:
            self.set_title(_arr[0][0], self.title[0])


        plt.show()
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure(fig)


    def flip(self, items, ncol):
        return itt.chain(*[items[i::ncol] for i in range(ncol)])
 
 
 
    def plot_mode_decomposition(self):
 
       file_qty = len(self.zpr_fnames) 
 
       # Define figure
       fig, _arr = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)

       if self.split:
           fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
           plt.subplots_adjust(hspace=0.05, top=0.95) 

       if not self.band_numbers:
           raise Exception('Must provide valence and conduction band numbers vis band_numbers')


       # Read and treat all input files
       for ifile, zpr_file in enumerate(self.zpr_fnames):
                
           # Define file class
           self.zpr = ZPRfile(zpr_file, read=False)

           # Read input file
           self.read_file()

           # Set parameters for this file
           self.nsppol = self.zpr.nsppol
           self.nkpt = self.zpr.nkpt
           self.max_band = self.zpr.max_band
           self.ntemp = self.zpr.ntemp
           self.temp = self.zpr.temp
           self.kpoints = self.zpr.kpoints
           self.nmodes =self.zpr.nmodes
           
           self.gap_units = self.zpr.gap_energy_units

           if self.units == 'eV':
               self.eig0 = self.zpr.eig0*cst.ha_to_ev
               self.eigcorr_modes = self.zpr.eigcorr_modes*cst.ha_to_ev
           if self.units == 'meV':
               self.eig0 = self.zpr.eig0*cst.ha_to_ev*1000
               self.eigcorr_modes = self.zpr.eigcorr_modes*cst.ha_to_ev*1000
             

           if self.split:
                if self.indirect:
                    self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
#                    self.gap_location = self.zpr.indirect_gap_location_split
                    self.gap_ren = self.zpr.indirect_gap_ren_split

                    gap_index0 = np.zeros((2,2),dtype = int) # left/right, val/cond

                    for i,j in itt.product(range(2),range(2)):
                        gap_index0[i,j] = self.find_gap_index(self.loc0[i,j])
 
                else:
                    self.loc0 = self.zpr.unperturbed_gap_location_split
                    self.gap_location = self.zpr.gap_location_split
                    self.gap_ren = self.zpr.gap_renormalization_split

                    gap_index0 = np.zeros((2),dtype = int)
                    gap_index0[0] = self.find_gap_index(self.loc0[0])
                    gap_index0[1] = self.find_gap_index(self.loc0[1])
               
               
           else:
                if self.indirect:
                    self.loc0 = self.zpr.unperturbed_indirect_gap_location
                    self.gap_location = self.zpr.indirect_gap_location
                    self.gap_ren = self.zpr.indirect_gap_ren

                    gap_index0 = np.zeros((2),dtype = int)
                    gap_index0[0] = self.find_gap_index(self.loc0[0]) #VB
                    gap_index0[1] = self.find_gap_index(self.loc0[1]) #CB

                else:
                    self.loc0 = self.zpr.unperturbed_gap_location
                    self.gap_location = self.zpr.gap_location
                    self.gap_ren = self.zpr.gap_renormalization
                    gap_index0 = self.find_gap_index(self.loc0)

           if ifile==0:
               if self.split:
                   self.full_mode_ren = np.zeros((file_qty,2,3,self.nmodes)) #nfiles, left/right, val/cond/gap, nmodes
               else:
                   self.full_mode_ren = np.zeros((file_qty,3,self.nmodes))

           if self.split:
               for j,m in itt.product(range(2), range(self.nmodes)):
                   self.full_mode_ren[ifile,j,0,m] = self.eigcorr_modes[0,gap_index0[j,0],self.band_numbers[0]-1,m] 
                   self.full_mode_ren[ifile,j,1,m] = self.eigcorr_modes[0,gap_index0[j,1],self.band_numbers[1]-1,m] 
                   self.full_mode_ren[ifile,j,2,m] = self.eigcorr_modes[0,gap_index0[j,1],self.band_numbers[1]-1,m]-self.eigcorr_modes[0,gap_index0[j,0],self.band_numbers[0]-1,m] 
           else:
               for m in range(self.nmodes):
                   self.full_mode_ren[ifile,0,m] = self.eigcorr_modes[0,gap_index0[0],self.band_numbers[0]-1,m] 
                   self.full_mode_ren[ifile,1,m] = self.eigcorr_modes[0,gap_index0[1],self.band_numbers[1]-1,m] 
                   self.full_mode_ren[ifile,2,m] = self.eigcorr_modes[0,gap_index0[1],self.band_numbers[1]-1,m] - self.eigcorr_modes[0,gap_index0[0],self.band_numbers[0]-1,m] 

           print(gap_index0[0,0],self.band_numbers[0]-1, gap_index0[0,1],self.band_numbers[1]-1)
           print(zpr_file)
           print(np.sum(self.full_mode_ren[ifile,0,0,:]), np.sum(self.full_mode_ren[ifile,0,1,:]))

           # Plot by mode number...
           x = range(self.nmodes)+np.ones((self.nmodes))
    
           if self.split:
               #conduction band
               _arr[0][0].plot(x, self.full_mode_ren[ifile,0,1,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
               _arr2[0][0].plot(x, self.full_mode_ren[ifile,1,1,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])

               #Valence band
               _arr[1][0].plot(x, self.full_mode_ren[ifile,0,0,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
               _arr2[1][0].plot(x, self.full_mode_ren[ifile,1,0,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])

               #Gap correction
               _arr[2][0].plot(x, self.full_mode_ren[ifile,0,2,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
               _arr2[2][0].plot(x, self.full_mode_ren[ifile,1,2,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
   
   
           else:
               _arr[0][0].plot(x,self.full_mode_ren[ifile,1,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
               _arr[1][0].plot(x,self.full_mode_ren[ifile,0,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
               _arr[2][0].plot(x,self.full_mode_ren[ifile,2,:], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])


       for i in range(3):
           self.set_vrefs(_arr[i][0], x, 0.)

           if self.split:
               self.set_vrefs(_arr2[i][0], x, 0.)

       
       _arr[2,0].set_xlabel('mode number')
       self.set_yaxis(_arr[0][0], 'CB ren ({})'.format(self.gap_units))
       self.set_yaxis(_arr[1][0], 'VB ren ({})'.format(self.gap_units))
       self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))


       if self.split:
           self.set_xaxis(_arr2[2][0], self.temp)
           self.set_yaxis(_arr2[0][0], 'CB ren ({})'.format(self.gap_units))
           self.set_yaxis(_arr2[1][0], 'VB ren ({})'.format(self.gap_units))
           self.set_yaxis(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units))

           self.set_title(_arr[0][0], self.title[0])
           self.set_title(_arr2[0][0], self.title[1])

       else:
           self.set_title(_arr[0][0], self.main_title)

##

       self.set_legend_gap(_arr[2][0])
       self.set_main_title(fig) 

       if self.split:
           self.set_legend_gap(_arr2[2][0])

           self.set_title(_arr[0][0], self.title[0])
           self.set_title(_arr2[0][0], self.title[1])

           self.set_main_title(fig2)

       if self.split:
           self.save_figure_split(fig,fig2)
       else:
           self.save_figure(fig)

       plt.show()
 



#   plot valence band and conduction band renormalisation, at gap kpoint
    def plot_pgap_separate(self):

        file_qty = len(self.zpr_fnames)

        # Define figure
        fig, _arr = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
            plt.subplots_adjust(hspace=0.05, top=0.95) 

        if not self.band_numbers:
            raise Exception('Must provide valence and conduction band numbers vis band_numbers')


        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
            
            self.gridsize = self.zpr.gridsize
            self.gap_units = self.zpr.gap_energy_units

            if self.units == 'eV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev
            if self.units == 'meV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev*1000
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev*1000
               

            if self.split:
                self.loc0 = self.zpr.unperturbed_gap_location_split
                self.gap_location = self.zpr.gap_location_split
                self.gap_ren = self.zpr.gap_renormalization_split

                gap_index0 = np.zeros((2),dtype = int)
                gap_index0[0] = self.find_gap_index(self.loc0[0])
                gap_index0[1] = self.find_gap_index(self.loc0[1])
                
                
            else:
                self.loc0 = self.zpr.unperturbed_gap_location
                self.gap_location = self.zpr.gap_location
                self.gap_ren = self.zpr.gap_renormalization
                gap_index0 = self.find_gap_index(self.loc0)


            if ifile==0:
                if self.split:
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2,2))
                    self.ref_temp = self.temp
                else:
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2))
                    self.ref_temp = self.temp
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

            if self.split:
                self.full_gap_ren[ifile,:,:] = self.gap_ren
 #               print(np.shape(self.gap_ren))
 #               print(np.shape(self.full_gap_ren))
            else:
                self.full_gap_ren[ifile,:] = self.gap_ren

            for T in range(self.ntemp):
                if self.split:
                    for a in range(2):
                        gap_index = self.find_gap_index(self.gap_location[T,a])
                        self.full_band_ren[ifile, T,0,a] = self.eigcorr[0,gap_index,self.band_numbers[0]-1,T]-self.eig0[0,gap_index0[a],self.band_numbers[0]-1] 
                        self.full_band_ren[ifile, T,1,a] = self.eigcorr[0,gap_index,self.band_numbers[1]-1,T]-self.eig0[0,gap_index0[a],self.band_numbers[1]-1]

                        diff = (self.full_band_ren[ifile,T,1,a]-self.full_band_ren[ifile,T,0,a]) - self.full_gap_ren[ifile,T,a]
#                        if abs(diff)>1E-6:
#                            print(ifile,T,a)
#                            print(self.full_gap_ren[ifile,T,a],self.full_band_ren[ifile,T,0,a],self.full_band_ren[ifile,T,1,a],self.full_band_ren[ifile,T,1,a]-self.full_band_ren[ifile,T,0,a])
#                            print(self.gap_location[T,a], self.kpoints[gap_index])
                else:
                    gap_index = self.find_gap_index(self.gap_location[T])
                    self.full_band_ren[ifile, T,0] = self.eigcorr[0,gap_index,self.band_numbers[0]-1,T]-self.eig0[0,gap_index0,self.band_numbers[0]-1]
                    self.full_band_ren[ifile, T,1] = self.eigcorr[0,gap_index,self.band_numbers[1]-1,T]-self.eig0[0,gap_index0,self.band_numbers[1]-1]
                    diff = (self.full_band_ren[ifile,T,1]-self.full_band_ren[ifile,T,0]) - self.full_gap_ren[ifile,T]
#                    if abs(diff)>1E-6:
#                        print(ifile,T)
#                        print(self.full_gap_ren[ifile,T],self.full_band_ren[ifile,T,0],self.full_band_ren[ifile,T,1],self.full_band_ren[ifile,T,1]-self.full_band_ren[ifile,T,0])
#                        print(self.gap_location[T], self.kpoints[gap_index])

#                print(ifile, T, gap_index)

        if self.crit_pressure is not None:
            crit_index = self.find_temp_index()
        else:
            crit_index = None
#        print(crit_index)

        for T in range(self.ntemp):

            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    #figure1
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,0], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,0], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    self.set_legend_pgap2(_arr[2][0])
                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,0], marker='d', linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,0], marker='d', linewidth=1.5, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='d', linewidth=1.5, color=self.color[T])

                    #figure2
                    _arr2[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,1], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,1], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    self.set_legend_pgap2(_arr2[2][0])
                    _arr2[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,1], marker='d', linewidth=1.5, color=self.color[T])
                    _arr2[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,1], marker='d', linewidth=1.5, color=self.color[T])
                    _arr2[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='d', linewidth=1.5, color=self.color[T])

                else:
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    self.set_legend_pgap2(_arr[2][0])
                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1], marker='d', linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0], marker='d', linewidth=1.5, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='d', linewidth=1.5, color=self.color[T])


            else:
                if self.split:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)

                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)

                else:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)


        for i in range(3):
            self.set_vrefs(_arr[i][0], self.pressure, 0.)
            self.set_hrefs(self.ylims, _arr[i][0], self.crit_pressure,'black')

            if self.split:
                self.set_vrefs(_arr2[i][0], self.pressure, 0.)
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'black')

        
        self.set_xaxis(_arr[2][0], self.pressure)
        self.set_yaxis(_arr[0][0], 'CB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))


        if self.split:
            self.set_xaxis(_arr2[2][0], self.pressure)
            self.set_yaxis(_arr2[0][0], 'CB ren ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[1][0], 'VB ren ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units))

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])

        else:
            self.set_title(_arr[0][0], self.main_title)

        plt.show()
        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure(fig)
 
#   plot valence band and conduction band renormalisation, at VB max and CB min
    def plot_pgap_separate_indirect(self):

        file_qty = len(self.zpr_fnames)

        # Define figure
        fig, _arr = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
            plt.subplots_adjust(hspace=0.05, top=0.95) 

        if not self.band_numbers:
            raise Exception('Must provide valence and conduction band numbers vis band_numbers')


        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
            
#            self.gridsize = self.zpr.gridsize
            self.gap_units = self.zpr.gap_energy_units

            if self.units == 'eV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev
            if self.units == 'meV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev*1000
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev*1000
               
            # Set side of phase transition for split2 option
            if self.split2:
                if self.pressure[ifile] < self.crit_pressure:
                    side=0
                else:
                    side=1

            # Retrieve data from file
            if self.split:
                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
#                self.gap_location = self.zpr.gap_location_split
                self.gap_ren = self.zpr.indirect_gap_ren_split
                self.band_ren = self.zpr.indirect_gap_ren_band_split

                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                for a in range(2):
                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

            if self.split2:
                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
#                self.gap_location = self.zpr.indirect_gap_location_split
                self.gap_ren = self.zpr.indirect_gap_ren_split
                self.band_ren = self.zpr.indirect_gap_ren_band_split # T left/right vb/cb

                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                for a in range(2):
                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

 
            else:
#                self.loc0 = self.zpr.unperturbed_gap_location
#                self.gap_location = self.zpr.gap_location
                self.gap_ren = self.zpr.indirect_gap_ren
                self.band_ren = self.zpr.indirect_gap_ren_band
                gap_index0 = np.zeros((2),dtype=int)
                for a in range(2):
                    gap_index0[a] = self.find_gap_index(self.loc0[a])


            # Initialize plotting arrays
            if ifile==0:
                if self.split: # distinct data for left and rifht gap
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2,2)) # file, temp, vbcb, lr 
                    self.ref_temp = self.temp
                if self.split2: # single array, left gap for trovoal phase and right gap for topol phase
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2))
                    self.ref_temp = self.temp
                else:
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2))
                    self.ref_temp = self.temp
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

            if self.split:
                self.full_gap_ren[ifile,:,:] = self.gap_ren
     #           print(np.shape(self.gap_ren))
     #           print(np.shape(self.full_gap_ren))
            if self.split2:
                self.full_gap_ren[ifile,:] = self.gap_ren[:,side]
            else:
                self.full_gap_ren[ifile,:] = self.gap_ren

            for T in range(self.ntemp):
                if self.unpert:
                    if self.split:
                        for a in range(2):
                            self.full_band_ren[ifile, T,0,a] = self.eigcorr[0,gap_index0[0,a],self.band_numbers[0]-1,T]-self.eig0[0,gap_index0[0,a],self.band_numbers[0]-1] 
                            self.full_band_ren[ifile, T,1,a] = self.eigcorr[0,gap_index0[1,a],self.band_numbers[1]-1,T]-self.eig0[0,gap_index0[1,a],self.band_numbers[1]-1]
                    if self.split2:
                            self.full_band_ren[ifile, T,0] = self.eigcorr[0,gap_index0[0,0],self.band_numbers[0]-1,T]-self.eig0[0,gap_index0[0,0],self.band_numbers[0]-1] 
                            self.full_band_ren[ifile, T,1] = self.eigcorr[0,gap_index0[1,0],self.band_numbers[1]-1,T]-self.eig0[0,gap_index0[1,0],self.band_numbers[1]-1]

                    else:
                            self.full_band_ren[ifile, T,0] = self.eigcorr[0,gap_index0[0],self.band_numbers[0]-1,T]-self.eig0[0,gap_index0[0],self.band_numbers[0]-1] 
                            self.full_band_ren[ifile, T,1] = self.eigcorr[0,gap_index0[1],self.band_numbers[1]-1,T]-self.eig0[0,gap_index0[1],self.band_numbers[1]-1]

                else:
                    if self.split:
                        for a in range(2):
                            self.full_band_ren[ifile, T,0,a] = self.band_ren[T,a,0]  
                            self.full_band_ren[ifile, T,1,a] = self.band_ren[T,a,1]

                            diff = (self.full_band_ren[ifile,T,1,a]-self.full_band_ren[ifile,T,0,a]) - self.full_gap_ren[ifile,T,a]
    #                        if abs(diff)>1E-6:
    #                            print(ifile,T,a)
    #                            print(self.full_gap_ren[ifile,T,a],self.full_band_ren[ifile,T,0,a],self.full_band_ren[ifile,T,1,a],self.full_band_ren[ifile,T,1,a]-self.full_band_ren[ifile,T,0,a])
    #                            print(self.gap_location[T,a], self.kpoints[gap_index])
                    if self.split2:
                        self.full_band_ren[ifile, T, :] = self.band_ren[T,side,:]
                    else:
                        self.full_band_ren[ifile, T,0] = self.band_ren[T,0]
                        self.full_band_ren[ifile, T,1] = self.band_ren[T,1] 
                        diff = (self.full_band_ren[ifile,T,1]-self.full_band_ren[ifile,T,0]) - self.full_gap_ren[ifile,T]
    #                    if abs(diff)>1E-6:
    #                        print(ifile,T)
    #                        print(self.full_gap_ren[ifile,T],self.full_band_ren[ifile,T,0],self.full_band_ren[ifile,T,1],self.full_band_ren[ifile,T,1]-self.full_band_ren[ifile,T,0])
    #                        print(self.gap_location[T], self.kpoints[gap_index])

    #                print(ifile, T, gap_index)

        if self.crit_pressure is not None:
            crit_index = self.find_temp_index()
        else:
            crit_index = None
#        print(crit_index)

        for T in range(self.ntemp):

            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    #figure1
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    _arr[2][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.27,1.62), ncol = 1, fontsize=20)

                  #  self.set_legend_pgap2(_arr[2][0])
                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T])

                    #figure2
                    _arr2[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    self.set_legend_pgap2(_arr2[2][0])
                    _arr2[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                    _arr2[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                    _arr2[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T])

                else:
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    _arr[2][0].legend(numpoints=1, loc= 'center right', bbox_to_anchor=(1.27,1.62), ncol=1, fontsize=20)

                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', markersize=8, linewidth=2.0, color=self.color[T])


            else:
                if self.split:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    _arr[1][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.02,1.0), ncol = 1, fontsize=26)


#                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)

                else:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)


        limms = [[-50.,40.],[-2.,50.],[-70.,40.]]
        for i in range(3):
#            lims = _arr[i,0].get_ylim() 
            ylims = limms[i]    
            self.set_vrefs(_arr[i][0], self.pressure, 0.)
            self.set_hrefs(ylims, _arr[i][0], self.crit_pressure,'black')
            self.set_hrefs(ylims, _arr[i][0], self.crit_pressure+0.2,'black')
            _arr[i][0].fill([2.08,2.28,2.28,2.08],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
            _arr[i,0].set_ylim(limms[i])


            if self.split:
                self.set_vrefs(_arr2[i][0], self.pressure, 0.)
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'black')
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure+0.2,'black')
                lims = _arr[i,0].get_ylim() 
                _arr2[i][0].fill([2.08,2.28,2.28,2.08],[lims[0]-5,lims[0]-5,lims[1]+5,lims[1]+5],'gray',alpha=0.2)


        
        self.set_xaxis(_arr[2][0], self.pressure)
        self.set_yaxis(_arr[0][0], 'CB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))

        fig.subplots_adjust(left=0.11,bottom=0.08,right=0.81,top=0.95,wspace=0.2,hspace=0.12)

        #_arr[1,0].text(self.explicit_pressures[8],300, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=30)
        #_arr[1,0].text(self.explicit_pressures[35],300, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=30)

        #_arr[1,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=28, weight='bold', color='gray')



        if self.split:
            self.set_xaxis(_arr2[2][0], self.pressure)
            self.set_yaxis(_arr2[0][0], 'CB ren ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[1][0], 'VB ren ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units))

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])

        else:
            self.set_title(_arr[0][0], self.main_title)

        plt.show()
        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure(fig)

    def plot_split_contribution(self):
        # Plot total contribution to TDR, splitted into VB and CB contributions. Can also be used with individual modes.
        file_qty = len(self.zpr_fnames)

        # Define figure
        fig, _arr = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
            plt.subplots_adjust(hspace=0.05, top=0.95) 

        if not self.band_numbers:
            raise Exception('Must provide valence and conduction band numbers vis band_numbers')


        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
            self.nmodes = self.zpr.nmodes
            self.nqpt = self.zpr.nqpt
            
            #self.gridsize = self.zpr.gridsize
            self.gap_units = self.zpr.gap_energy_units

            if self.units == 'eV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev
            if self.units == 'meV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev*1000
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev*1000
               

            if self.split:
                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
#                self.gap_location = self.zpr.gap_location_split
                self.gap_ren = self.zpr.indirect_gap_ren_split
                self.band_ren = self.zpr.indirect_gap_ren_band_split
                print(self.band_ren[:,0,:])

                self.loc = self.zpr.indirect_gap_location_split #temp lr vbcb cart
#                print(self.loc[:,1,1])

                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                gap_index = np.zeros((self.ntemp,2,2),dtype=int) #( temp vbcb lr)
                for a in range(2):
                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

                    for t in range(self.ntemp):
                        gap_index[t,0,a] = self.find_gap_index(self.loc[t,a,0]) # temp vbcb lr
                        gap_index[t,1,a] = self.find_gap_index(self.loc[t,a,1])
#                print(gap_index[:,1,0])
                
#                print(gap_index[:,:,0])
            else:
#                self.loc0 = self.zpr.unperturbed_gap_location
#                self.gap_location = self.zpr.gap_location
                self.gap_ren = self.zpr.indirect_gap_ren
                self.band_ren = self.zpr.indirect_gap_ren_band
                gap_index0 = np.zeros((2),dtype=int)
                for a in range(2):
                    gap_index0[a] = self.find_gap_index(self.loc0[a])

            if self.split_contribution:
                self.fan_g2 = self.zpr.fan_g2
                self.ddw_g2 = self.zpr.ddw_g2
                self.deltaE_ddw = self.zpr.deltaE_ddw
                self.fan_occterm = self.zpr.fan_occterm
                self.qpt_weight = self.zpr.qpt_weight
                self.ddw_tdep = self.zpr.ddw_tdep



            if ifile==0:
                if self.split:
                    if self.modes:
                        self.full_gap_ren = np.zeros((file_qty, self.ntemp,2,2,self.nmodes)) # file temp lr vbcontr/cbcontr mode
                        self.full_band_ren = np.zeros((file_qty,self.ntemp,2,2,2,self.nmodes)) # file, temp, vbcb, lr vbcontr/cbcontr mode
                        self.ref_temp = self.temp
                    else:
                        self.full_gap_ren = np.zeros((file_qty, self.ntemp,2,2)) # file temp lr vbcontr/cbcontr
                        self.full_band_ren = np.zeros((file_qty,self.ntemp,2,2,2)) # file, temp, vbcb, lr vbcontr/cbcontr
                        self.ref_temp = self.temp
                else:
                    if self.modes:
                        print('modes not done yet')
                    else:
                        self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
                        self.full_band_ren = np.zeros((file_qty,self.ntemp,2,2))
                        self.ref_temp = self.temp
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

            # Copy function from print_coupling_info / get_se_modes... 
            # sum on qpts and modes (add a keyword NOT to sum on modes)
            # plot same as before but loop of dim 2 (vb contr, CB contr)
            # in the same type of raph as before

#            if self.split:
            split_se = self.get_split_se()

            if self.modes:
                for t,n,m,d,v in itt.product(range(self.ntemp), range(2), range(2), range(2),range(self.nmodes)):

                    self.full_band_ren[ifile,t,n,d,m,v] = split_se[0,gap_index[t,n,d],n,m,v,t]

            else:
                for t,n,m,d in itt.product(range(self.ntemp), range(2), range(2), range(2)):

                    self.full_band_ren[ifile,t,n,d,m] = split_se[0,gap_index[t,n,d],n,m,t]

            if self.modes:
                for t,m,d,v in itt.product(range(self.ntemp), range(2), range(2),range(self.nmodes)):
                    self.full_gap_ren[ifile,t,d,m,v] = self.full_band_ren[ifile,t,1,d,m,v] - self.full_band_ren[ifile,t,0,d,m,v]
            else:
                for t,m,d in itt.product(range(self.ntemp), range(2), range(2)):
                    self.full_gap_ren[ifile,t,d,m] = self.full_band_ren[ifile,t,1,d,m] - self.full_band_ren[ifile,t,0,d,m]
     #           print(np.shape(self.gap_ren))
     #           print(np.shape(self.full_gap_ren))
           # else:
#                self.full_gap_ren[ifile,:] = self.get_split_se()
            print(self.pressure[ifile])
            print(self.full_band_ren[ifile,:,1,0,0])
            print(self.full_band_ren[ifile,:,1,0,1])
            print(self.full_band_ren[ifile,:,1,0,0]+self.full_band_ren[ifile,:,1,0,1])


        if self.crit_pressure is not None:
            crit_index = self.find_temp_index()
        else:
            crit_index = None
#        print(crit_index)

        style = ['dashed', 'dashdot']

        for T in range(self.ntemp):

            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    #figure1
                    for d in range(2):
                        _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,0,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                        _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,0,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                        _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')

                        self.set_legend_pgap2(_arr[2][0])
                        _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,0,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                        _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,0,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                        _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T])

                        #figure2
                        _arr2[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,1,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                        _arr2[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,1,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                        _arr2[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')

                        self.set_legend_pgap2(_arr2[2][0])
                        _arr2[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,1,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                        _arr2[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,1,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
                        _arr2[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T])

                else:
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,d], linestyle=style[d], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    self.set_legend_pgap2(_arr[2][0])
                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,d], linestyle=style[d], marker='d', linewidth=2.0, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,d], linestyle=style[d], marker='d', linewidth=2.0, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,d], linestyle=style[d], marker='d', linewidth=2.0, color=self.color[T])


            else:
                if self.split:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)

                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)

                else:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)


        for i in range(3):
            self.set_vrefs(_arr[i][0], self.pressure, 0.)
            self.set_hrefs(self.ylims, _arr[i][0], self.crit_pressure,'black')
            self.set_hrefs(self.ylims, _arr[i][0], self.crit_pressure+0.2,'black')

            lims = _arr[i,0].get_ylim() 
            _arr[i][0].fill([2.08,2.28,2.28,2.08],[lims[0],lims[0],lims[1],lims[1]],'gray',alpha=0.2)
            _arr[i,0].set_ylim(lims)


            if self.split:
                self.set_vrefs(_arr2[i][0], self.pressure, 0.)
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'black')
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure+0.2,'black')
                lims = _arr[i,0].get_ylim() 
                _arr2[i][0].fill([2.08,2.28,2.28,2.08],[lims[0]-5,lims[0]-5,lims[1]+5,lims[1]+5],'gray',alpha=0.2)


        
        self.set_xaxis(_arr[2][0], self.pressure)
        self.set_yaxis(_arr[0][0], 'CB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))


        if self.split:
            self.set_xaxis(_arr2[2][0], self.pressure)
            self.set_yaxis(_arr2[0][0], 'CB ren ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[1][0], 'VB ren ({})'.format(self.gap_units))
            self.set_yaxis(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units))

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])

        else:
            self.set_title(_arr[0][0], self.main_title)

        plt.show()
        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure(fig)


    def get_split_se(self):

        # Reconstruct self-energy, splitted into VB and CB contributions

        # Fan
        fan2 = np.zeros((self.nsppol, self.nkpt, 2, 2, self.nmodes,self.nqpt), dtype= complex)
        fan2.real = self.fan_g2[:,:,:,:,:,:,0]
        fan2.imag = self.fan_g2[:,:,:,:,:,:,1]

        fan_occ = np.zeros((self.nsppol, self.nkpt, 2, 2, self.nmodes,self.ntemp, self.nqpt), dtype= complex)
        fan_occ.real = self.fan_occterm[:,:,:,:,:,:,:,0]
        fan_occ.imag = self.fan_occterm[:,:,:,:,:,:,:,1]


        if self.modes:
            fan = np.einsum('sknmvq,sknmvtq->sknmvtq',fan2,fan_occ)
        else:
            fan = np.einsum('sknmvq,sknmvtq->sknmtq',fan2,fan_occ)

        del fan2, fan_occ

        #Debye-Waller
        ddw2 = np.zeros((self.nsppol, self.nkpt, 2, 2, self.nmodes,self.nqpt), dtype= complex)
        ddw2.real = self.ddw_g2[:,:,:,:,:,:,0]
        ddw2.imag = self.ddw_g2[:,:,:,:,:,:,1]

        deno_ddw = np.zeros((self.nsppol,self.nkpt, 2, 2, self.nqpt), dtype= complex)
        deno_ddw.real = self.deltaE_ddw[:,:,:,:,:,0]
        deno_ddw.imag = self.deltaE_ddw[:,:,:,:,:,1]

        ddw_occterm = np.einsum('sknmq,vtq->sknmvtq',1./deno_ddw,self.ddw_tdep)

        if self.modes:
            ddw = np.einsum('sknmvq,sknmvtq->sknmvtq',ddw2, ddw_occterm)
        else:
            ddw = np.einsum('sknmvq,sknmvtq->sknmtq',ddw2, ddw_occterm)

        del ddw2, deno_ddw, ddw_occterm

        total = fan-ddw
        if self.modes:
            se = np.zeros((self.nsppol,self.nkpt,2,2,self.nmodes,self.ntemp))
        else:
            se = np.zeros((self.nsppol,self.nkpt,2,2,self.ntemp))

        
        #add qpt weight and sum on qpoints 
        if self.modes:
            for q in range(self.nqpt):
#                se += self.qpt_weight[q]*(total[:,:,:,:,:,:,q].real)

                if self.units=='meV':
                    se += self.qpt_weight[q]*(total[:,:,:,:,:,:,q])*cst.ha_to_ev*1000
                else:
                    se += self.qpt_weight[q]*(total[:,:,:,:,:,:,q])*cst.ha_to_ev

        else:
            for q in range(self.nqpt):
                se_q = self.qpt_weight[q]*(total[:,:,:,:,:,q].real)

                se_q = self.make_average(se_q)

#                if self.units=='meV':
#                    se += self.qpt_weight[q]*(total[:,:,:,:,:,q].real)*cst.ha_to_ev*1000
                se += se_q
#                else:
#                    se += self.qpt_weight[q]*(total[:,:,:,:,:,q].real)*cst.ha_to_ev

        return se*cst.ha_to_ev


    def make_average(self,arr):

         """ 
         Average a quantity over degenerated states.
         Does not work with spin yet.
     
         Arguments
         ---------
     
         arr: numpy.ndarray(..., nkpt, nband)
             An array of any dimension, of which the two last indicies are
             the kpoint and the band.
     
         Returns
         -------
     
         arr: numpy.ndarray(..., nkpt, nband)
             The array with the values of degenerated bands averaged.
     
         """
 
         if not self.degen:
             self.get_degen()
 
         nkpt, nband = arr.shape[-2:]
 
         for ikpt in range(nkpt):
             for group in self.degen[ikpt]:
                 average = copy(arr[...,ikpt,group[0][1]])
                 for ispin, iband in group[1:]:
                     average += arr[...,ikpt,iband]
 
                 average /= len(group)
                 for ispin, iband in group:
                     arr[...,ikpt,iband] = average
 
         return arr

    def get_degen(self):
#    """
#     Compute the degeneracy of the bands.
#  
#       Returns
#        -------
#         degen: 2D list (nkpt, )
#             For each k-point, contains a list of groups of (s, n) tuples
#             which are degenerated.
#             For example, if there is only a single spin (spin unpolarized case)
#             and two k-points, and at the first k-point there are two triplets,
#             then the output is
#                 [[[(0,1), (0,2), (0,3)],  [(0,4), (0,5), (0,6)]], []]
#    
#     """
        nspin, nkpt, nband = self.eig0.shape
    
        degen = list()
        for ikpt in range(nkpt):
    
            kpt_degen = list()
            group = list()
            last_ispin, last_iband, last_eig = 0, 0, -float('inf')
    
            for sbe in self.iter_spin_band_eig(ikpt):
                ispin, iband, eig = sbe
    
                if np.isclose(last_eig, eig, rtol=1e-12, atol=1e-5):
                    if not group:
                        group.append((last_ispin, last_iband))
                    group.append((ispin, iband))
    
                else:
                    if group:
                        kpt_degen.append(group)
                        group = list()
    
                last_ispin, last_iband, last_eig = ispin, iband, eig
    
            degen.append(kpt_degen)
    
        self.degen = degen
    
        return degen

    def iter_spin_band_eig(self, ikpt):
#        """
#         Iterator over spin index, band index, and eigenvalues at one k-point.
#         Yields tuples (ispin, iband, eig) in order of increasing eig.
#         """
         nspin, nkpt, nband = self.eig0.shape
         sbe = [(ispin, 0, self.eig0[ispin,ikpt,0]) for ispin in range(nspin)]
         cmp_sbe = lambda sbe1, sbe2: cmp(sbe1[2], sbe2[2])
         while sbe:
             min_sbe = sorted(sbe, cmp=cmp_sbe)[0]
             yield min_sbe
 
             i = sbe.index(min_sbe)
             ispin, iband, eig = min_sbe
 
             if iband == nband - 1:
                 del sbe[i]
             else:
                 sbe[i] = (ispin, iband+1, self.eig0[ispin, ikpt, iband+1])
 


    def plot_self_energy(self):
        
        plot_qty = len(self.bands_to_print)
        file_qty = len(self.zpr_fnames)

        if file_qty > len(self.color):
            raise Exception('color array must have at least {} entries, but it has only {}'.format(file_qty,len(self.color)))

        if self.units == 'eV':
            self.fermi = self.fermi*cst.ha_to_ev
        elif self.units == 'meV':
            self.fermi = self.fermi*1000

        fig, arr = plt.subplots( 1, plot_qty, squeeze = False, sharey = True,  figsize = self.figsize)

        # Read and treat all input files
        for ifile, zpr_file in enumerate(self.zpr_fnames):
                
            # Define file class
            self.zpr = ZPRfile(zpr_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.zpr.nsppol
            self.nkpt = self.zpr.nkpt
            self.max_band = self.zpr.max_band
            self.ntemp = self.zpr.ntemp
            self.temp = self.zpr.temp
            self.kpoints = self.zpr.kpoints
            self.nfreq = self.zpr.nfreq
    
            self.omega_se = self.zpr.omega_se
            self.self_energy = self.zpr.self_energy
            self.spectral_function = self.zpr.spectral_function
            self.broadening = self.zpr.broadening
            self.smearing = self.zpr.smearing

            if self.units is 'eV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev
                self.self_energy = self.self_energy*cst.ha_to_ev
                self.omega_se = self.omega_se*cst.ha_to_ev

            elif self.units is 'meV':
                self.eig0 = self.zpr.eig0*cst.ha_to_ev*1000
                self.eigcorr = self.zpr.eigcorr*cst.ha_to_ev*1000
                self.ylims = self.ylims*1000
                self.self_energy = self.self_energy*cst.ha_to_ev*1000
                self.omega_se = self.omega_se*cst.ha_to_ev*1000


            for iplot in range(plot_qty):
    
                arr[0,iplot].plot(self.omega_se, self.self_energy[0,self.point_for_se, self.bands_to_print[iplot]-1, :, 0], color=self.color[ifile],label=self.title[ifile], linewidth=1.5)

                self.set_xaxis(arr[0,iplot],self.omega_se)
                self.set_vrefs(arr[0,iplot],self.omega_se, 0.)



        if not self.ylims:

            self.ylims = arr[0,iplot].get_ylim()


        # Set axis properties
        self.set_yaxis(arr[0,0],'Re[Self-energy] ({})'.format(self.units))

        for i in range(plot_qty):
            self.set_legend_se(arr[0,i])
            self.set_hrefs(self.ylims, arr[0,i],0.,'black')

        if plot_qty == 2:
            self.set_title(arr[0,0],'Valence band')
            self.set_title(arr[0,1], 'Conduction band')

        self.set_main_title(fig)
        self.save_figure(fig)

        plt.show()

    def plot_spectral_function(self):
        print('spectral')

    def find_gap_index(self, loc):

        loclst = list(loc)
        lst = [list(x) for x in self.kpoints]

        return lst.index(loclst)

    def find_temp_index(self):

        lst = list(self.temp) 

        if self.vbcb:
            if self.temp_to_print in lst:
                return lst.index(self.temp_to_print)
            else:
                return None
        if self.renormalization:
            if self.tmax in lst:
                return lst.index(self.tmax)
            else:
                return None
        if self.pgap:
            lst = list(self.pressure)
            if self.crit_pressure in lst:
                return lst.index(self.crit_pressure)
            else:
                for i in range(len(self.pressure)-1):
                    if self.pressure[i]<self.crit_pressure and self.crit_pressure < self.pressure[i+1]:
                        return i
                return None

    def find_temp_index2(self, pressure):

        lst = list(pressure) 
        if self.crit_pressure in lst:
            return lst.index(self.crit_pressure)
        else:
            for i in range(len(pressure)-1):
                if pressure[i]<self.crit_pressure and self.crit_pressure < pressure[i+1]:
                    return i
            return None

    def find_pressure_index(self, p):

        lst = list(self.explicit_pressures) 
        if p in lst:
            return lst.index(p)
        else:
            for i in range(len(self.explicit_pressures)-1):
                if self.explicit_pressures[i]<p and p< self.explicit_pressures[i+1]:
#                    print("in between")
                    return i
            return None




    def set_xaxis(self,f,x):

        if self.renormalization:
            f.set_xlim((0,len(x)-1))
            if self.xticks is not None:
                plt.setp(f.get_xticklabels(), fontsize=18,weight='bold')
                plt.setp(f,xticks=self.xticks[0], xticklabels=self.xticks[1])
    
                if self.xticks_alignment is not None:
                    for itick, tick in enumerate(f.xaxis.get_majorticklabels()):
                        tick.set_horizontalalignment(self.xticks_alignment[itick]) 

        if self.senergy:
            f.set_xlabel('w-e^0 ({})'.format(self.units), fontsize=18)
            if self.xlims:
                f.set_xlim(self.xlims[0], self.xlims[-1])
            else:
                f.set_xlim(self.omega_se[0], self.omega_se[-1])

        if self.gap:
            f.set_xlabel('Temperature (K)', fontsize=18)
            f.xaxis.set_tick_params(labelsize=16)
            if not self.xlims:
                self.xlims = [min(x), max(x)]
            f.set_xlim(self.xlims)

        if self.vbcb:
            f.set_xlim((0,len(x)-1))
            if self.xticks is not None:
                plt.tick_params(axis='x',labelsize = 18)
#                plt.setp(f.get_xticklabels(), fontsize='large')
                plt.setp(f,xticks=self.xticks[0], xticklabels=self.xticks[1])
    
                if self.xticks_alignment is not None:
                    for itick, tick in enumerate(f.xaxis.get_majorticklabels()):
                        tick.set_horizontalalignment(self.xticks_alignment[itick]) 

        if self.pgap:
            lims = (min(self.pressure)-0.1, max(self.pressure)+0.1)
            f.set_xlim(lims)
            plt.setp(f.get_xticklabels(), fontsize=20, weight='bold')
            if self.separate_bands:
                f.set_xlabel('Pressure (GPa)', fontsize=20)
            else:
                f.set_xlabel('Pressure (GPa)',fontsize=20)

    def set_yaxis(self,f, lab):

        if self.vbcb:
            f.set_ylabel(lab, fontsize=14)
            
        else:
            f.set_ylabel(lab, fontsize = 20)


            if self.ylims is not None:
                f.set_ylim(self.ylims) 
            else:
                lims = f.get_ylim()
                f.set_ylim(lims)

        if self.yticks is not None:
            plt.setp(f.get_yticklabels(), fontsize=20)
            f.set_yticks(np.arange(self.yticks[0], self.yticks[1]+self.yticks[2], self.yticks[2]))
        else:
            plt.setp(f.get_yticklabels(), fontsize=20, weight='bold')
        if self.yminorticks:
            f.yaxis.set_minor_locator(AutoMinorLocator(self.yminorticks+1))

    def set_ylimits_vbcb(self,f,lims):

        f.set_ylim(lims)

            
    def adjust_ylimits(self, f, a, t):

        #Get current ylims
        lims = f.get_ylim()

        # Correct if necessary
        if lims[0] < self.ylimits[a][0]:
            self.ylimits[a][0] = lims[0]
        if lims[1] > self.ylimits[a][1]:
            self.ylimits[a][1] = lims[1]


    def set_hrefs(self, lims, f, val,col):

        print(lims)
        if not lims:
            lims = f.get_ylim()
#        print(lims[0], lims[1])
        # Set vertical line at reference point (omega-eps^0=0 for self energy, unperturbed gap energy for gap, 0 for gap renormalization)
        y = np.linspace(lims[0], lims[1], 10)
        zer = val*np.ones(len(y))

#        if self.vbcb == True:
#            f.plot(zer,y,'k:')
#        else:
        f.plot(zer,y,col)


    def set_vrefs(self,f,x,val):

        if self.pgap:
            lims = (min(self.pressure)-0.1, max(self.pressure)+0.1)
            x = np.linspace(lims[0], lims[1], 50)
            zer = val*np.ones((len(x)))

        else:
            # Set constant horizontal line at reference value
            zer = val*np.ones(len(x))
    #        if self.vbcb==True:
    #            f.plot(x,zer,'k')
    #        else:
        f.plot(x, zer, color='black', linestyle='solid')


    def set_legend_pgap(self,f):

        box = f.get_position()
        f.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
        f.legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.27), ncol = self.ntemp+1, fontsize=20)

    def set_legend_pgap2(self,f):

        box = f.get_position()
        if self.split:
#            f.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
            f.legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.40), ncol = self.ntemp+1, fontsize=22)
#           f.legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.40), ncol = 3, fontsize=22)
        else:
 #           f.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
            f.legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.40), ncol = self.ntemp+1, fontsize=22)

    def set_legend_ren(self, f, a):

        labels = ['unperturbed']

        if self.subplots:
            if self.temp[a] <= self.tmax:
                labels.append(str(self.temp[a])+'K')
        else:
            for T in self.temp:
                if T <= self.tmax:
                    labels.append(str(T)+'K')
            
            f.legend(labels,loc = 'lower right',bbox_to_anchor=(1.01,-0.15 ),ncol=len(labels),fontsize=16)
#            f.legend(labels,loc = 'center right',bbox_to_anchor=(1.01,-0.15 ),ncol=1,fontsize=16)


    def set_legend_ren2(self, f, lines):

        labels = ['unperturbed']

#        if self.subplots:
#            if self.temp[a] <= self.tmax:
#                labels.append(str(self.temp[a])+'K')
#        else:
        for T in self.temp:
            if T <= self.tmax:
                labels.append(str(T)+'K')
            
        f.legend(lines,labels,loc=9, bbox_to_anchor=(0.5,-0.1),ncol=len(labels))

    def set_legend_se(self,f):

        if self.title:
            f.legend()

    def set_legend_gap(self, f):
        #f.legend(numpoints=1, loc=3, fontsize=16)
        f.legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.10,2.0), ncol = 1, fontsize=12)


    def set_legend_vbcb(self, f):
#        box = f.get_position()
#        f.set_position([box.x0, box.y0 + box.height*0.1, box.width, box.height*0.9])
#        f.legend(loc='upper center',bbox_to_anchor=(0.5,-0.5))
        f.legend()

    def set_title(self, f, j):

        if self.vbcb:
            f.set_title(j, fontsize=16) 
        elif self.title:
            f.set_title(j, fontsize=24)
        else:
            f.set_title(j, fontsize=24)

    def set_main_title(self,g):
        # Set main title
        if self.main_title:
            if self.senergy:
                t = self.main_title+'\nKpt : '+str(self.kpoints[self.point_for_se])
            else:
                t = self.main_title
            g.suptitle(t, fontsize=18)

    def save_figure(self,g):

        create_directory('figures/')

        if self.savefile:
            g.savefig('figures/'+self.savefile+'.png')

    def save_figure2(self,g):

        create_directory('figures/')

        g.savefig('figures/unperturbed_energy.png')

    def save_figure_split(self,g1,g2):

        create_directory('figures/')

        if self.savefile:
            g1.savefig('figures/{}_lgap.png'.format(self.savefile))
            g2.savefig('figures/{}_rgap.png'.format(self.savefile))


##########################
# Create a directory if it does not exist
def create_directory(fname):

    dirname = os.path.dirname(fname)
    if not dirname:
        return
    if not os.path.exists(dirname):
        os.system('mkdir -p ' + dirname)

# Main function

def plotter(
        #Input file
        zpr_fnames = list(),
        gap_fname = None,
        rootname = 'zpr.png',

        # Parameters
        nsppol = None,
        nkpt = None,
        max_band = None,
        eig0 = None,
        eigcorr = None,
        ntemp = None,
        temp = None,
        units = 'eV',
        fermi = None,
        fermi_td = None,
        tmax = None,
        scissor = None,
        kpoints = None,
        pressure = None,
        crit_pressure = None,

        # Options
        figsize = (15,9),
        color = None,
        labels = None,
        linestyle = None,
        ylims = None,
        yticks = None,
        yminorticks = None,
        xlims = None,
        xticks = None,
        xticks_alignment = None,
        bands_to_print = None,
        band_numbers = None,
        point_for_se = None,
        temp_to_print = 0.,
        subplots = False,
        main_title = None,
        title = None,
        savefile = None,

        # Type of plot
        renormalization = True,
        spectral = False,
        senergy = False,
        broad = False,
        gap = False,
        pgap = False,
        vbcb = False,
        separate_bands = False,
        split = False,
        split2 = False,
        follow = False,
        indirect = False,
        unpert = False,
        split_contribution = False,
        modes = False,
        verbose = False,
        zero_gap_value = None,
        zero_gap_units = 'eV',
        experimental_data = None,

        **kwargs):

        zpr_plot = ZPR_plotter(
            zpr_fnames = zpr_fnames,
            rootname = rootname,
            gap_fname = gap_fname,

            nsppol = nsppol,
            nkpt = nkpt,
            max_band = max_band,
            kpoints = kpoints,
            ntemp = ntemp,
            temp = temp,
            eig0 = eig0,
            eigcorr = eigcorr,
            units = units,

            color = color,
            labels = labels,
            linestyle = linestyle,
            fermi = fermi,
            fermi_td = fermi_td,
            tmax = tmax,
            scissor = scissor,
            pressure = pressure,
            crit_pressure = crit_pressure,

            xlims = xlims,
            xticks = xticks,
            xticks_alignment = xticks_alignment,
            ylims = ylims,
            yticks = yticks,
            yminorticks = yminorticks,
            figsize = figsize,

            bands_to_print = bands_to_print,
            band_numbers = band_numbers,
            point_for_se = point_for_se,
            temp_to_print = temp_to_print,
            zero_gap_value = zero_gap_value,
            zero_gap_units = zero_gap_units,
            experimental_data = experimental_data,

            subplots = subplots,
            main_title = main_title,
            title = title,
            savefile = savefile,

            renormalization = renormalization,
            senergy = senergy,
            spectral = spectral,
            broad = broad,
            gap = gap,
            pgap = pgap,
            vbcb = vbcb,
            separate_bands = separate_bands,
            split = split,
            split2 = split2,
            follow = follow,
            indirect = indirect,
            unpert = unpert,
            split_contribution = split_contribution,
            modes = modes,
            verbose = verbose,

            **kwargs)
            
    
        # Plot the renormalized bandstructure
        if renormalization:
            zpr_plot.plot_zpr()

        #Plot the self-energy
        if senergy:
            zpr_plot.plot_self_energy()
        
        # Plot the spectral function
        if spectral:
            zpr_plot.plot_spectral_function()

        # Plot the gap as a function of temperature
        if gap:
            if separate_bands:
                zpr_plot.plot_gap_separate()
            else:
                zpr_plot.plot_gap()

        # plot valence band and conduction band ZPR separately
        if vbcb:
            zpr_plot.plot_vbcb()

        
        if modes:
            zpr_plot.plot_mode_decomposition()

        # Plot gap (T) as a function of pressure
        if pgap:
            if indirect:
                if separate_bands:
                    if split_contribution:
                        zpr_plot.plot_split_contribution()
                    else:
                        zpr_plot.plot_pgap_separate_indirect()
                else:
#                    zpr_plot.plot_pgap_indirect_only()
                    zpr_plot.plot_pgap_indirect()

            else:
                if separate_bands:
                    zpr_plot.plot_pgap_separate()
                else:
                    zpr_plot.plot_pgap()

        

        return zpr_plot







