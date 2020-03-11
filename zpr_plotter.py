#! /usr/bin/env python
import matplotlib
matplotlib.use('TKAgg')
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib import cm
import matplotlib.colors as mcolors

import netCDF4 as nc
import numpy as np
import os
import warnings
import itertools as itt
from copy import copy
from scipy.interpolate import interp1d
############################
# this code plots the corrected bandstructure on a symmetrized kpath obtained with zpr_analyser.py

# Inputfile : <rootname>_ZPR.nc
# Output : png image

# V. Brousseau 2017/07
###########################

# Tex settings
rc('text', usetex = True)
#rc('font', family = 'serif',weight = 'bold')
rc('font', family = 'sans-serif',weight = 'bold')
#rc('font',family='sans-serif')
plt.rcParams['axes.unicode_minus'] = False
params = {'text.latex.preamble' : [r'\usepackage{amsmath,amssymb}']}
#r'\usepackage[utf8]{inputenc}',r'\DeclareUnicodeCharacter{2212}{$-$}',
#r'\usepackage{cmbright}']}
#         r'\renewcommand{\familydefault}{\sfdefault}']}
plt.rcParams.update(params)

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

class TEfile(CDFfile):
        
    def __init__(self, *args, **kwargs):

        super(TEfile, self).__init__(*args, **kwargs)

        self.eig0 = None
        self.kpoints = None
        self.temp = None
        self.eigcorr = None
        self.kpoints_spline = None

    # Open the ZPR.nc file and read it
    def read_nc(self, fname = None):
        
        fname = fname if fname else self.fname
        super(TEfile, self).read_nc(fname)

        with nc.Dataset(fname, 'r') as ncdata:

            self.temp = ncdata.variables['temperatures'][:]
            self.correction = ncdata.variables['eigenvalue_corrections'][:,:,:,:]

            self.eig0 = ncdata.variables['bare_eigenvalues'][:,:,:]
            self.kpoints = ncdata.variables['reduced_coordinates_of_kpoints'][:,:]
            self.unperturbed_gap_location = ncdata.variables['unperturbed_gap_location'][:]
            self.unperturbed_gap_energy = ncdata.variables['unperturbed_gap_energy'][:]
            self.gap_location = ncdata.variables['gap_location'][:,:]
            self.gap_energy = ncdata.variables['gap_energy'][:]
            self.gap_energy_units = ncdata.variables['gap_energy'].getncattr('units')
            self.gap_renormalization = ncdata.variables['direct_gap_te_renormalization'][:]
#            self.band_energy = ncdata.variables['band_energy'][:,:]
#            self.band_renormalization = ncdata.variables['band_renormalization'][:,:]


#            self.unperturbed_gap_location_split = ncdata.variables['unperturbed_gap_location_split'][:,:]
            #self.unperturbed_gap_energy_split = ncdata.variables['unperturbed_gap_energy_split'][:,:]
 #           self.gap_location_split = ncdata.variables['gap_location_split'][:,:,:]
 #           self.gap_energy_split = ncdata.variables['gap_energy_split'][:,:]
 #           self.gap_energy_units_split = ncdata.variables['gap_energy_split'].getncattr('units')
 #           self.gap_renormalization_split = ncdata.variables['gap_renormalization_split'][:,:]
#            self.band_energy_split = ncdata.variables['band_energy_split'][:,:,:]
#            self.band_renormalization_split = ncdata.variables['band_renormalization_split'][:,:,:]

            self.unperturbed_indirect_gap_location = ncdata.variables['unperturbed_indirect_gap_location'][:,:]
            self.indirect_gap_location = ncdata.variables['indirect_gap_location'][:,:,:]
            self.unperturbed_indirect_gap_energy = ncdata.variables['unperturbed_indirect_gap_energy'][:]
            self.indirect_gap_energy = ncdata.variables['indirect_gap_energy'][:]
            self.indirect_gap_ren = ncdata.variables['indirect_gap_te_renormalization'][:]
#            self.indirect_gap_energy_band = ncdata.variables['indirect_gap_energy_band'][:,:]
            self.indirect_gap_ren_band = ncdata.variables['te_renorm_indirect_band'][:,:]


            
#            self.unperturbed_indirect_gap_location_split = ncdata.variables['unperturbed_indirect_gap_location_split'][:,:,:]
#            self.indirect_gap_location_split = ncdata.variables['indirect_gap_location_split'][:,:,:,:]
#            self.unperturbed_indirect_gap_energy_split = ncdata.variables['unperturbed_indirect_gap_energy_split'][:,:]
#            self.indirect_gap_energy_split = ncdata.variables['indirect_gap_energy_split'][:,:]
#            self.indirect_gap_ren_split = ncdata.variables['indirect_gap_renormalization_split'][:,:]
#            self.indirect_gap_energy_band_split = ncdata.variables['indirect_gap_energy_band_split'][:,:,:]
#            self.indirect_gap_ren_band_split = ncdata.variables['indirect_gap_renormalization_band_split'][:,:,:]


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
        self.kpoints_spline = None

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

            self.fan_occ = ncdata.variables['reduced_fan_occ'][:,:,:,:]
            self.fan_unocc = ncdata.variables['reduced_fan_unocc'][:,:,:,:]
            self.ddw_occ = ncdata.variables['reduced_ddw_occ'][:,:,:,:]
            self.ddw_unocc = ncdata.variables['reduced_ddw_unocc'][:,:,:,:]

            status = ncdata.get_variables_by_attributes(name='corrected_eigenvalue_spline_interpolation')
            if status != []:
                self.eigcorr_spline = ncdata.variables['corrected_eigenvalue_spline_interpolation'][:,:,:,:]
            status = ncdata.get_variables_by_attributes(name='eigenvalue_corrections_spline_interpolation')
            if status != []:
                self.td_ren_spline = ncdata.variables['eigenvalue_corrections_spline_interpolation'][:,:,:,:]
            status = ncdata.get_variables_by_attributes(name='reduced_coordinated_of_kpoints_spline')
            if status != []:
                self.kpoints_spline = ncdata.variables['reduced_coordinated_of_kpoints_spline'][:,:]


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
    def nkpt_spline(self):
        return np.int(self.kpoints_spline.shape[0]) if self.kpoints_spline is not None else None

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
    te_fnames = None
    bare_fnames = None
    gap_fname = None
    deltap_fname = None
        
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
    crit_pressure2 = None
    explicit_pressures = None
    valence = None
    spline = False
    spline_factor = 5
    
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
    cond_ylims = None
    val_ylims = None
    gap_ylims = None
    egap_ylims = None
    cond_yticks = None
    val_yticks = None
    gap_yticks = None
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
    cmap_cols = None

    renormalization = True
    spectral = False
    senergy = False
    broad = False
    gap = False
    pgap = False
    te_pgap = False
    total_pgap = False
    vbcb = False
    separate_bands = False
    phase_diagram = False
    split = False
    split2 = False
    follow = False
    indirect = False
    unpert = False
    split_contribution = False
    split_occupied_subspace = False
    modes = False
    verbose = False
    zero_gap_value = None
    experimental_data = None
    expdata = None
    zero_gap_units = 'eV'
    extrapolate_ren = False
    read_correction_from_spline = False
    crossover = False
    diagram_together = False
    save_format = 'png'

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
            cond_ylims = None,
            val_ylims = None,
            gap_ylims = None,
            egap_ylims = None,
            cond_yticks = None,
            val_yticks = None,
            gap_yticks = None,
            xlims = None,
            xticks = None,
            xticks_labels = None,
            xticks_alignment = None,
            cmap_cols = None,
            bands_to_print = None,
            band_numbers = None,
            point_for_se = None,
            temp_to_print = 0.,
            subplots = False,
            main_title = None,
            title = None,
            savefile = None,
            spline_factor = 5,
            spline = False,

            gap_loc0 = None,
            gap_energy0 = None,
            gap_loc = None,
            gap_energy = None,
            gap_ren = None,
            valence = None,


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
            te_pgap = False,
            total_pgap = False,
            separate_bands = False,
            phase_diagram = False,
            split = False,
            split2 = False,
            follow=False,
            indirect = False,
            unpert = False,
            split_contribution = False,
            split_occupied_subspace = False,
            modes = False,
            degen = None,
            verbose = False,
            zero_gap_value = None,
            zero_gap_units = 'eV',
            extrapolate_ren = False,
            read_correction_from_spline = False,

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
            crit_pressure2 = None,
            explicit_pressures = None,
            experimental_data = None,
            expdata = None,
            crossover = False,
            diagram_together = False,

            units = 'eV',
            gap_units = None,
            rootname = 'zpr.png',
            save_format = 'png',

            # Input file
            zpr_fnames = list(),
            te_fnames = None,
            gap_fname = None,
            bare_fnames = None,
            deltap_fname = None,
            
            **kwargs):

        # Set input file
        self.zpr_fnames = zpr_fnames
        self.gap_fname = gap_fname
        self.te_fnames = te_fnames
        self.bare_fnames = bare_fnames
        self.deltap_fname = deltap_fname


        # Define options
        self.figsize = figsize
        self.color = color
        self.labels = labels
        self.linestyle = linestyle
        self.ylims = ylims
        self.yticks = yticks
        self.yminorticks = yminorticks
        self.cond_ylims = cond_ylims
        self.val_ylims = val_ylims
        self.gap_ylims = gap_ylims
        self.egap_ylims = egap_ylims
        self.cond_yticks = cond_yticks
        self.val_yticks = val_yticks
        self.gap_yticks = gap_yticks

        self.xlims = xlims
        self.xticks = xticks
        self.xticks_labels = xticks_labels
        self.xticks_alignment = xticks_alignment
        self.cmap_cols = cmap_cols
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
        self.crit_pressure2 = crit_pressure2
        self.valence = valence
        self.spline = spline
        self.spline_factor = spline_factor
        self.extrapolate_ren = extrapolate_ren
        self.read_correction_from_spline = read_correction_from_spline
        self.crossover = crossover
        self.diagram_together = diagram_together

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
        self.te_pgap = te_pgap
        self.total_pgap = total_pgap
        self.separate_bands = separate_bands
        self.phase_diagram = phase_diagram
        self.split = split
        self.split2 = split2
        self.follow = follow
        self.indirect = indirect
        self.unpert = unpert
        self.split_contribution = split_contribution
        self.split_occupied_subspace = split_occupied_subspace
        self.modes = modes
        self.save_format = save_format

        # Check if input is correct 
        self.check_input()

        self.set_rootname(rootname)

        if self.experimental_data is not None:
            self.expdata = EXPfile(self.experimental_data)
            self.expdata.read_nc()
    
    def check_input(self):
        
        # Checks that all required input is provided
        if not self.zpr_fnames:
            if not self.te_pgap:
                raise Exception('Must provide files for zpr_fnames')


        if sum(map(bool,[self.renormalization, self.senergy, self.spectral])) > 1:
            raise Exception('Please specify ONLY ONE quantity to plot from renormalization, self-energy and spectral function')

        if self.senergy:
            if self.point_for_se is None:
                raise Exception('Must provide a kpoint index to treat in point_for_se')  
    
            if not self.bands_to_print:
                raise Exception('Must provide band numbers for self-energy in bands_to_print') 

        if self.units is not 'meV' and self.units is not 'eV':
            raise Exception('Units must be eV or meV')

        if self.te_pgap:
            if not self.te_fnames:
                raise Exception('Must provide files for te_fnames')


    def read_file(self):
    
        if self.te_pgap:
            f = self.te
        else:
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
            self.spline = True

            kpt_array = np.arange(self.nkpt)
            if self.spline:
                self.nkpt_fine = 5*self.nkpt
                kpt_array_fine = np.linspace(kpt_array[0],kpt_array[-1],self.nkpt_fine)
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


            if self.spline:

                nband = len(self.bands_to_print)

                # Find index of A point 
                kpt_is_a = self.find_gap_index([0.0, 0.0, 0.5])
                if kpt_is_a is None:
                    kpt_is_a = self.find_gap_index([0.0, 0.0,0.53849])

                #First, treat the bare eigenvalues
                self.eig0_int = np.zeros((self.nsppol, self.nkpt_fine,nband))
                eig0 = np.zeros((self.nsppol, self.nkpt,nband))
                eig0[0,:,0] = np.concatenate((self.eig0[0,:kpt_is_a,self.bands_to_print[0]-1],self.eig0[0,kpt_is_a:,self.bands_to_print[1]-1]))
                eig0[0,:,1] = np.concatenate((self.eig0[0,:kpt_is_a,self.bands_to_print[1]-1],self.eig0[0,kpt_is_a:,self.bands_to_print[0]-1]))
                eig0[0,:,2] = np.concatenate((self.eig0[0,:kpt_is_a,self.bands_to_print[2]-1],self.eig0[0,kpt_is_a:,self.bands_to_print[3]-1]))
                eig0[0,:,3] = np.concatenate((self.eig0[0,:kpt_is_a,self.bands_to_print[3]-1],self.eig0[0,kpt_is_a:,self.bands_to_print[2]-1]))

                for iband, band in enumerate(self.bands_to_print): 
                    spl0 = interp1d(kpt_array,eig0[0,:,iband], kind='cubic')
                    self.eig0_int[0,:,iband] = spl0(kpt_array_fine)

                # Then, the Tdep ones
                self.eigcorr_int = np.zeros((self.nsppol, self.nkpt_fine, nband, self.ntemp))
                for t in range(self.ntemp):
                    eig = np.zeros((self.nsppol, self.nkpt,nband))
                    eig[0,:,0] = np.concatenate((self.eigcorr[0,:kpt_is_a,self.bands_to_print[0]-1,t],self.eigcorr[0,kpt_is_a:,self.bands_to_print[1]-1,t]))
                    eig[0,:,1] = np.concatenate((self.eigcorr[0,:kpt_is_a,self.bands_to_print[1]-1,t],self.eigcorr[0,kpt_is_a:,self.bands_to_print[0]-1,t]))
                    eig[0,:,2] = np.concatenate((self.eigcorr[0,:kpt_is_a,self.bands_to_print[2]-1,t],self.eigcorr[0,kpt_is_a:,self.bands_to_print[3]-1,t]))
                    eig[0,:,3] = np.concatenate((self.eigcorr[0,:kpt_is_a,self.bands_to_print[3]-1,t],self.eigcorr[0,kpt_is_a:,self.bands_to_print[2]-1,t]))

                    for iband, band in enumerate(self.bands_to_print): 
                        spl0 = interp1d(kpt_array,eig[0,:,iband], kind='cubic')
                        self.eigcorr_int[0,:,iband,t] = spl0(kpt_array_fine)
                
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

                if self.spline == True:
                    if self.bands_to_print:
                        for jband,iband in enumerate(self.bands_to_print):
                            arr[iplot][ifile].plot(kpt_array_fine, self.eig0_int[0,:,jband]-self.fermi_td[ifile]+scissor_array[iband-1], color='k',linewidth=1.5)
                            
                            if self.subplots:
                                arr[iplot][ifile].plot(kpt_array_fine, self.eigcorr_int[0,:,jband,iplot]-self.fermi_td[ifile]+scissor_array[iband-1], color=self.color[iplot],linewidth=1.5)
                            else:
                                for itemp in range(tmax_index+1):
                                    arr[iplot][ifile].plot(kpt_array_fine, self.eigcorr_int[0,:,jband,itemp]-self.fermi_td[ifile]+scissor_array[iband-1], color = self.color[itemp],linewidth=1.5)
            
                    else:
                        for iband in range(self.max_band):
                            arr[iplot][ifile].plot(kpt_array_fine, self.eig0_int[0,:,iband]-self.fermi_td[ifile]+scissor_array[iband], color='k', label='Static T=0')
                            
                            if self.subplots:
                                arr[iplot][ifile].plot(kpt_array_fine, self.eigcorr_int[0,:,jband,iplot]-self.fermi_td[ifile]+scissor_array[iband], color = self.color[iplot], label=str(self.temp[iplot])+'K')
                            else:
                                for itemp in range(tmax_index+1):
                                    arr[iplot][ifile].plot(kpt_array_fine, self.eigcorr_int[0,:,jband,itemp]-self.fermi_td[ifile]+scissor_array[iband], color = self.color[itemp], label=str(self.temp[itemp])+'K')
     
                else:
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
                            arr[iplot][ifile].plot(kpt_array, self.eig0[0,:,iband]-self.fermi_td[ifile]+scissor_array[iband], color='k', label='Static T=0')
                            
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
                if self.linestyle[ifile] == 'dashed':
                    #conduction band
                    _arr[0][0].plot(self.temp, self.band_ren[:,0,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (4,4))
                    _arr2[0][0].plot(self.temp, self.band_ren[:,1,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (4,4))
                    #Valence band
                    _arr[1][0].plot(self.temp, self.band_ren[:,0,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (4,4))
                    _arr2[1][0].plot(self.temp, self.band_ren[:,1,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (4,4))
                    #Gap correction
                    _arr[2][0].plot(self.temp, self.gap_ren[:,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (4,4))
                    _arr2[2][0].plot(self.temp, self.gap_ren[:,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (4,4))
                elif self.linestyle[ifile] == 'dotted':
                    #conduction band
                    _arr[0][0].plot(self.temp, self.band_ren[:,0,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (2,2))
                    _arr2[0][0].plot(self.temp, self.band_ren[:,1,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (2,2))
                    #Valence band
                    _arr[1][0].plot(self.temp, self.band_ren[:,0,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (2,2))
                    _arr2[1][0].plot(self.temp, self.band_ren[:,1,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (2,2))
                    #Gap correction
                    _arr[2][0].plot(self.temp, self.gap_ren[:,0], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (2,2))
                    _arr2[2][0].plot(self.temp, self.gap_ren[:,1], marker='o', linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile], dashes = (2,2))
                else:
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
                _arr[0][0].plot(self.temp, self.band_ren[:,1], marker='o', markersize=8, linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                _arr[1][0].plot(self.temp, self.band_ren[:,0], marker='o', markersize=8,linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])
                _arr[2][0].plot(self.temp, self.gap_ren, marker='o', markersize=8,linewidth=1.5, color=self.color[ifile], label=self.labels[ifile], linestyle=self.linestyle[ifile])

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
####### FROM HERE
        limms = [[-50.,50.],[-2.,50.],[-90.,40.]]
##########FIX ME : add default data and if condition
        limms = [self.cond_ylims, self.val_ylims, self.gap_ylims]
        for i in range(3):
#            lims = _arr[i,0].get_ylim() 
            ylims = limms[i]    
            self.set_vrefs(_arr[i][0], self.temp, 0.,style='dashed')
            _arr[i,0].set_ylim(limms[i])


            if self.split:
                self.set_vrefs(_arr2[i][0], self.temp, 0.)

        
        self.set_xaxis(_arr[2][0], self.temp)
        self.set_yaxis_separate(_arr[0][0], 'CBM ren ({})'.format(self.gap_units),self.cond_ylims, self.cond_yticks)
        self.set_yaxis_separate(_arr[1][0], 'VBM ren ({})'.format(self.gap_units),self.val_ylims, self.val_yticks)
        self.set_yaxis_separate(_arr[2][0], 'Gap ren ({})'.format(self.gap_units),self.gap_ylims, self.gap_yticks)


        self.set_yaxis(_arr[0][0], 'CBM ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VBM ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))

        fig.subplots_adjust(left=0.11,bottom=0.08,right=0.81,top=0.95,wspace=0.2,hspace=0.12)

        if self.split:
            self.set_xaxis(_arr2[2][0], self.temp)
            self.set_yaxis_separate(_arr2[0][0], 'CB ren ({})'.format(self.gap_units),self.cond_ylims, self.cond_yticks)
            self.set_yaxis_separate(_arr2[1][0], 'VB ren ({})'.format(self.gap_units),self.val_ylims, self.val_yticks)
            self.set_yaxis_separate(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units),self.gap_ylims, self.gap_yticks)

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])

        else:
            if self.main_title:
                self.set_title(_arr0[0][0], self.main_title)

        legend1 = _arr[0][0].legend(loc=9,bbox_to_anchor=(0.5,1.28),fontsize=24, handletextpad=0.4,handlelength=1.0,frameon=True,ncol = file_qty,columnspacing=0.7)
        _arr[0][0].add_artist(legend1)
        
        fig.subplots_adjust(hspace=0.0,top=0.91,right=0.90,bottom=0.10,left=0.12)

        if self.split:
            fig.align_ylabels()
            fig2.align_ylabels()
            self.save_figure_split(fig,fig2)
        else:
            fig.align_ylabels()
            self.save_figure(fig)

###### TO HERE
#
#        for i in range(3):
#            self.set_vrefs(_arr[i][0], self.temp, 0.)
#
#            if self.split:
#                self.set_vrefs(_arr2[i][0], self.temp, 0.)
#
#        
#        self.set_xaxis(_arr[2][0], self.temp)
#        self.set_yaxis(_arr[0][0], 'CB ren ({})'.format(self.gap_units))
#        self.set_yaxis(_arr[1][0], 'VB ren ({})'.format(self.gap_units))
#        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))
#
#
#        if self.split:
#            self.set_xaxis(_arr2[2][0], self.temp)
#            self.set_yaxis(_arr2[0][0], 'CB ren ({})'.format(self.gap_units))
#            self.set_yaxis(_arr2[1][0], 'VB ren ({})'.format(self.gap_units))
#            self.set_yaxis(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units))
#
#            self.set_title(_arr[0][0], self.title[0])
#            self.set_title(_arr2[0][0], self.title[1])
#
#        else:
#            if self.main_title is not None:
#                self.set_title(_arr[0][0], self.main_title)
#
####
#
#        self.set_legend_gap(_arr[2][0])
#        if self.main_title is not None:
#            self.set_main_title(fig) 
#
#        if self.split:
#            self.set_legend_gap(_arr2[2][0])
#
#            self.set_title(_arr[0][0], self.title[0])
#            self.set_title(_arr2[0][0], self.title[1])
#
#            if self.main_title is not None:
#                self.set_main_title(fig2)
#
#        if self.split:
#            self.save_figure_split(fig,fig2)
#        else:
#            self.save_figure(fig)

#        plt.show()
 
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
                       
                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black', label='Static T=0')
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
                       
                        _arr2[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
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
                       
                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='Static T=0')
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
        only=True  # only the Egap(P) for each T. If False, add gap renorm on top subplot
        extrapolate_bands = True # for now, default extrapolates on the TOTAL GAP REN, not bands separately...
        #when exprapolation on bands work properly, remove the other option. This should be the default. 

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

            
        ####### Define critical pressure index, for plotting                                                                                                                                            
        if self.crit_pressure is not None:
            crit_index = self.find_crit_pressure(self.crit_pressure, self.pressure)
        else:
            crit_index = None
    #        print(crit_index)

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

                if self.unpert:
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

                if self.unpert:
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
                self.band_ren = self.zpr.indirect_gap_ren_band

                if self.unpert:
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
                    self.full_band_ren = np.zeros((file_qty, self.ntemp,2))

                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((file_qty))
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

            # Extrapolate both band behavior towards the TPT
            if ifile==0:                                                                                                                                                                        
                if extrapolate_bands:
                    if self.split:
                        raise Exception('split is not implemented yet for extrapolation')
                    elif self.split2:
                        raise Exception('split2 is not implemented yet for extrapolation')
                    else:
                        print('Extrapolating towards TPT')
                        # Adjust these values if needed!
                        pressure1 = np.arange(1.5,2.0,0.1)
                        tmparr = np.arange(2.1,2.6,0.1)
                        pressure2 = np.arange(2.3,3.0,0.1)
                        tmparr2 = np.arange(2.0,2.3,0.1)

                        '''why do I do this?!?'''
                        self.extr_pressure1 = np.concatenate((pressure1, tmparr))
                        self.extr_pressure2 = np.concatenate((tmparr2,pressure2))

                        self.extr_full_band_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp,2)) #pressure, temperature, val/cond
                        self.extr_full_band_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp,2))

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
                self.full_band_ren[ifile,:,:] = self.band_ren
                if not self.gap_fname:
                    self.full_energy0[ifile] = self.energy0
                else:
                    p_index = self.find_pressure_index(self.pressure[ifile])
                    #print(p_index)
                    self.full_energy0[ifile] = self.explicit_gap_energies[p_index]

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

        #### End of loop on files. All data is stored now.

#        if self.crit_pressure is not None:
#            crit_index = self.find_temp_index()
        if self.gap_fname is not None:
#            crit_index2 = self.find_temp_index2(self.explicit_pressures)
            crit_index2 = self.find_crit_pressure(self.crit_pressure2,self.explicit_pressures)
            print('crit_index2',crit_index2)
#        else:
#            crit_index = None


        if self.follow:
            self.full_energy0[-1,:] = -self.full_energy0[-1,:]
            self.full_gap_energy[-1,:,:] = -self.full_gap_energy[-1,:,:]

        # Extrapolate behavior linearly towards TPT... it IS a bit sketchy.
        # Default, extrapolate on gap energy (not renormalization)
        if not extrapolate_bands:
            pressure1 = np.arange(1.5,2.0,0.1)
            tmparr = np.arange(2.1,2.6,0.1) # change 2.2 to 2.6 foir full BZ
            self.extr_pressure2 = np.arange(2.3, 3.0, 0.1)
            self.extr_pressure1 = np.concatenate((pressure1, tmparr))
        # For extrabolate_bands, it was done earlier.

        # linearize unperturbed gap energy to extrapolate inside WSM phase
        y0,y1 = np.polyfit(self.explicit_pressures[15:20], self.explicit_gap_energies[15:20],1)
        tmp2 = y1 + tmparr*y0


        if not extrapolate_bands:
            #for extrapolate bands, most of it was done earlier
            self.extr_full_gap_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp))
            self.extr_full_gap_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp))
            self.extr_full_gap_energy1 = np.zeros((len(self.extr_pressure1), self.ntemp))
            self.extr_full_gap_energy2 = np.zeros((len(self.extr_pressure2), self.ntemp))
            extr_energy01 = np.zeros((len(pressure1)))
            self.extr_full_energy02 = np.zeros((len(self.extr_pressure2)))
        else:
            extr_energy01 = np.zeros((len(pressure1)))
            self.extr_full_energy02 = np.zeros((len(self.extr_pressure2)))
            self.extr_full_gap_energy1 = np.zeros((len(self.extr_pressure1), self.ntemp))
            self.extr_full_gap_energy2 = np.zeros((len(self.extr_pressure2), self.ntemp))


        # Get full_gap_energy0 from explicit data
        ''' add  a  if explicit_energies is not None???'''
        for p, pres in enumerate(pressure1):
            ind = self.find_pressure_index(pres)
            extr_energy01[p] = self.explicit_gap_energies[ind]

        self.extr_full_energy01 = np.concatenate((extr_energy01,tmp2))

        for p, pres in enumerate(self.extr_pressure2):
            ind = self.find_pressure_index(pres)
            self.extr_full_energy02[p] = self.explicit_gap_energies[ind]

        #Store extrapolated phase boundaries (T)
        self.pc1 = np.zeros((self.ntemp))
        self.pc2 = np.zeros((self.ntemp))

        for T in range(self.ntemp):

            if extrapolate_bands:
                 if self.split:
                     raise Exception('split is not implemented yet for extrapolation!!! How did you get here?!?')
                 elif self.split2:
                     raise Exception('split2 is not implemented yet for extrapolation!!! How did you get here?!?')
                 else:

                     #loop on val and cond extremas
                     for b in range(2):
                         #trivial side
                         x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_band_ren[crit_index-1:crit_index+1,T,b], 1)
                         self.extr_full_band_ren1[:,T,b] = x1 + x0*self.extr_pressure1
                         # topol side
                         x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_band_ren[crit_index+1:crit_index+3,T,b], 1)
                         self.extr_full_band_ren2[:,T,b] = x1 + x0*self.extr_pressure2
                    
                     self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_band_ren1[:,T,1]-self.extr_full_band_ren1[:,T,0]
                     '''first, check what the plot is giving. Then, figure out how to get the intersect...'''
                     self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_band_ren2[:,T,1]-self.extr_full_band_ren2[:,T,0]
                     x0,x1 = np.polyfit(self.extr_pressure1,self.extr_full_gap_energy1[:,T],1)
                     print('trivial side : {} GPa'.format(-x1/x0))
                     self.pc1[T] = -x1/x0
                     x0,x1 = np.polyfit(self.extr_pressure2,self.extr_full_gap_energy2[:,T],1)
                     print('topol side : {} GPa'.format(-x1/x0))
                     self.pc2[T] = -x1/x0

            else:
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


                else:
                    #trivial side
                    x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_gap_ren[crit_index-1:crit_index+1,T], 1)
                    print('For T={} K:'.format(self.temp[T]))
                    self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]
                    x0,x1 = np.polyfit(self.extr_pressure1,self.extr_full_gap_energy1[:,T],1)
                    print('trivial side : {} GPa'.format(-x1/x0))
                    self.pc1[T] = -x1/x0
                    # topol side
                    x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_gap_ren[crit_index+1:crit_index+3,T], 1)
                    self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]
                    x0,x1 = np.polyfit(self.extr_pressure2,self.extr_full_gap_energy2[:,T],1)
                    print('topol side : {} GPa'.format(-x1/x0))
                    self.pc2[T] = -x1/x0

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
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='Static T=0')
                            if self.crit_index2 is not None:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='Static T=0')

                            else:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='Static T=0')

                        else:
                            if only:
                                _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='Static T=0')
                            else:
                                _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='Static T=0')
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
                                    _arr[0][0].plot(self.explicit_pressures[crit_index2+2,:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[crit_index2+2:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label='Static T=0')

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
                                _arr2[0][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                            else:
                                arr2[1][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                        else:
                            if only:
                                _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                            else:
                                _arr2[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
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

#########################################
       #-------- crit_index is defined but no split or split2
                else:
                    if T==self.temp[0]:
                        ###### FIX ME : keep explicit_pressures and explicit_gap_energies, just get a new crit index for this array! 
                        if self.gap_fname is not None:
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='Static T=0')
                            if crit_index2 is not None:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[:crit_index2+1], self.explicit_gap_energies[:crit_index2+1],marker='d',markersize=7, color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:crit_index2+1], self.explicit_gap_energies[:crit_index2+1],marker='d',markersize=7, color='black', label='Static T=0')
                            else:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', markersize=7,color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d',markersize=7, color='black', label='Static T=0')

                        else:
                            if only:
                                _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black', label='Static T=0')
                            else:
                                _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black', label='Static T=0')
#                    else:
#                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black')
    
                    if only:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='d',markersize=7, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    else:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='d',markersize=7, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
        #                self.set_legend_gap(_arr[0][0])
                        _arr[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='d',markersize=7, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
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
                                     _arr[0][0].plot(self.explicit_pressures[crit_index2+1:], self.explicit_gap_energies[crit_index2+1:],marker='o',markersize=7, color='black', label=None)
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[crit_index2+1:], self.explicit_gap_energies[crit_index2+1:],marker='o',markersize=7, color='black', label=None)

                    if only:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o',markersize=7, linewidth=2.0, color=self.color[T])
                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        print(self.extr_full_gap_energy2[:,T])
                        print(self.extr_pressure2)

                    else:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o',markersize=7, linewidth=2.0, color=self.color[T])
                        _arr[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o',markersize=7, linewidth=2.0, color=self.color[T])

                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


                    if self.gap_fname is not None:
                        pass
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
                    raise Exception('Only keyword not implemented yet without critical pressure index')
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
        limms = [self.gap_ylims, self.egap_ylims]

        if only:
            ylims = limms[1]    
            self.set_vrefs(_arr[0][0], self.pressure, 0.)
#            if i==0:
            pc1 = self.crit_pressure
            pc2 = self.crit_pressure2
            self.set_hrefs(ylims, _arr[0][0], pc1,'gray')
            self.set_hrefs(ylims, _arr[0][0], pc2,'gray')
            _arr[0][0].fill([pc1,pc2,pc2,pc1],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
#            else:
#            self.set_hrefs(ylims, _arr[0,0], 2.18,'orange')
#            self.set_hrefs(ylims, _arr[0,0], 2.50,'orange')
#            _arr[0,0].fill([2.18,2.50,2.50,2.18],[ylims[0],ylims[0],ylims[1],ylims[1]],'orange',alpha=0.2)
#            self.set_hrefs(ylims, _arr[0,0], 2.68,'magenta')


        else:
            for i in range(2):
    #            lims = _arr[i,0].get_ylim() 
    #            ylims = lims
                ylims = limms[i]    
                self.set_vrefs(_arr[i][0], self.pressure, 0.)
    #            if i==0:
                pc1 = self.crit_pressure
                pc2 = self.crit_pressure2

                self.set_hrefs(ylims, _arr[i][0], pc1,'black')
                self.set_hrefs(ylims, _arr[i][0], pc2,'black')
                _arr[i][0].fill([pc1,pc2,pc2,pc1],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
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


#        if only:
#            _arr[0,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=24) 
#            _arr[0,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=24) 
#
##            _arr[0,0].text(self.explicit_pressures[5],290, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=26)
##            _arr[0,0].text(self.explicit_pressures[38],290, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=26)
#
#
#            _arr[0,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=26, weight='bold', color='gray')
#            _arr[0,0].text(self.explicit_pressures[25],250, r'WSM 300K', fontsize=26, weight='bold', color='orange')
#
#        else:
#            _arr[1,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=24) 
#            _arr[1,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=24) 
#
##            _arr[1,0].text(self.explicit_pressures[5],290, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=26)
##            _arr[1,0].text(self.explicit_pressures[38],290, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=26)
#
#
#            _arr[1,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=26, weight='bold', color='gray')
#            _arr[1,0].text(self.explicit_pressures[25],250, r'WSM 300K', fontsize=26, weight='bold', color='orange')



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

        # Custom text
        if only:
            _arr[0][0].text(0.7,250, r'$\mathbb{Z}_2\!=\!0$',fontsize=24,color='k')
            _arr[0][0].text(2.7,250, r'$\mathbb{Z}_2\!=\!1$',fontsize=24,color='k')
            _arr[0][0].text(1.90,200,r'$\Rightarrow$',fontsize=26,color='#5A5A5A')
            _arr[0][0].text(1.50,210,r'Static',fontsize=26,color='#5A5A5A',weight='bold')
            _arr[0][0].text(1.50,190,r'WSM',fontsize=26,color='#5A5A5A',weight='bold')

            legend_handles = []
            legend_handles.append(Line2D([0],[0],color='k',linewidth=1.5,label=r'0 K Static'))
            for t,T in enumerate(self.ref_temp):
                legend_handles.append(Line2D([0],[0],color=self.color[t],linewidth=1.5,label=r'{:>3.0f} K'.format(T)))

            legend1 = _arr[0][0].legend(handles=legend_handles, loc=9, bbox_to_anchor=(0.5,1.16),fontsize=24,handletextpad=0.4,handlelength=1.0,frameon=True,ncol=len(self.ref_temp)+1,columnspacing=0.7)
            _arr[0][0].add_artist(legend1)

            legend2_handles = []
            legend2_handles.append(Line2D([0],[0],marker='d',markersize=12,markerfacecolor='None',markeredgecolor='k', linestyle='None',label=r'P$_{\text{C1}}$ plane'))
            legend2_handles.append(Line2D([0],[0],marker='o',markersize=12,markerfacecolor='None',markeredgecolor='k', linestyle='None',label=r'P$_{\text{C2}}$ plane'))
            legend2 = _arr[0][0].legend(handles=legend2_handles, loc=1,bbox_to_anchor=(1.0,1.0),fontsize=24, handletextpad=0.4,handlelength=1.4,frameon=True,ncol=1,labelspacing=0.1)
            _arr[0][0].add_artist(legend2)

        else:
            print('Complete legend not implemented yet')


#        plt.show()
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig.align_ylabels()
            fig2.align_ylabels()
            self.save_figure_split(fig,fig2)
        else:
            fig.align_ylabels()
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
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='Static T=0')
                            _arr[0][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='Static T=0')


                        else:
                            _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='Static T=0')
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
                            _arr2[0][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                        else:
                            _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
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
                            _arr[0][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s],marker='d', color='black', label='Static T=0')
                        else:
                            _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s],marker='d', color='black', label='Static T=0')
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
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
            plt.subplots_adjust(hspace=0.05, top=0.95) 

        if not self.band_numbers:
            raise Exception('Must provide valence and conduction band numbers vis band_numbers')

####### Define critical pressure index, for plotting
        if self.crit_pressure is not None:
            crit_index = self.find_crit_pressure(self.crit_pressure, self.pressure)
        else:
            crit_index = None
#        print(crit_index)


        ########################
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

                if self.unpert:
                    gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                    for a in range(2):
                        gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                        gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

            if self.split2:
                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
#                self.gap_location = self.zpr.indirect_gap_location_split
                self.gap_ren = self.zpr.indirect_gap_ren_split
                self.band_ren = self.zpr.indirect_gap_ren_band_split # T left/right vb/cb

                if self.unpert:
                    gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                    for a in range(2):
                        gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                        gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

 
            else:
                self.loc0 = self.zpr.unperturbed_indirect_gap_location
#                self.gap_location = self.zpr.gap_location
                self.gap_ren = self.zpr.indirect_gap_ren
                self.band_ren = self.zpr.indirect_gap_ren_band
                if self.unpert:
                    gap_index0 = np.zeros((2),dtype=int)
                    for a in range(2):
                        gap_index0[a] = self.find_gap_index(self.loc0[a])


            # Initialize plotting arrays
            if ifile==0:
                if self.split: # distinct data for left and right gap
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2,2)) # file, temp, vbcb, lr 
                    self.ref_temp = self.temp
                if self.split2: # single array, left gap for trivial phase and right gap for topol phase
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

                # Extrapolate both band behavior towards the TPT
            if ifile==0:
                if self.extrapolate_ren:
                    if self.split:
                        raise Exception('split is not implemented yet for extrapolation')
                    elif self.split2:
                        raise Exception('split2 is not implemented yet for extrapolation')
                    else:
                        print('Extrapolating towards TPT')
                        # Adjust these values if needed!
                        pressure1 = np.arange(1.5,2.0,0.1)
                        tmparr = np.arange(2.1,2.6,0.1) 
                        pressure2 = np.arange(2.3,3.0,0.1)
                        tmparr2 = np.arange(2.0,2.3,0.1)

                        '''why do I do this?!?'''
                        self.extr_pressure1 = np.concatenate((pressure1, tmparr))
                        self.extr_pressure2 = np.concatenate((pressure2,tmparr2))

                        self.extr_full_gap_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp,2)) #pressure, temperature, val/cond
                        self.extr_full_gap_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp,2))
                        #self.extr_full_gap_ren_spline1 = np.zeros((len(self.extr_pressure1), self.ntemp,2)) #pressure, temperature, val/cond
                        #self.extr_full_gap_ren_spline2 = np.zeros((len(self.extr_pressure2), self.ntemp,2))


#############################################
            # Define full band ren
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
                        if abs(diff)>1E-6:
                            print('difference between full_band_ren and full_gap_ren from input file:')
                            print(ifile,T)
                            print(self.full_gap_ren[ifile,T],self.full_band_ren[ifile,T,0],self.full_band_ren[ifile,T,1],self.full_band_ren[ifile,T,1]-self.full_band_ren[ifile,T,0])
                            print(self.gap_location[T], self.kpoints[gap_index])

    # End loop on files. All data is stored now.
        if self.extrapolate_ren:
                if self.split:
                    raise Exception('split is not implemented yet for extrapolation!!! How did you get here?!?')
                elif self.split2:
                    raise Exception('split2 is not implemented yet for extrapolation!!! How did you get here?!?')
                else:

                    for T in range(self.ntemp):

                        #loop on val and cond extremas
                        print(self.ref_temp[T])
                        for b in range(2):
                            #trivial side
                            x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_band_ren[crit_index-1:crit_index+1,T,b], 1)
                            self.extr_full_gap_ren1[:,T,b] = x1 + x0*self.extr_pressure1
                            print(x0,x1)
        #                        self.extr_full_gap_energy1[:,T,b] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T,b]
        #                        x0,x1 = np.polyfit(self.extr_pressure1,self.extr_full_gap_energy1[:,T],1)
        #                        print('trivial side : {} GPa'.format(-x1/x0))
        #                        self.pc1[T] = -x1/x0
                            #ftriv = interp1d(self.pressure[:crit_index+1], self.full_band_ren[:crit_index+1,T,b], kind='slinear',bounds_error=False, fill_value='extrapolate')
                            #self.extr_full_gap_ren_spline1[:,T,b] = ftriv(self.extr_pressure1)

                            # topol side
                            x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_band_ren[crit_index+1:crit_index+3,T,b], 1)
                            self.extr_full_gap_ren2[:,T,b] = x1 + x0*self.extr_pressure2
                            print(x0,x1)

        #                        self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]
        #                        x0,x1 = np.polyfit(self.extr_pressure2,self.extr_full_gap_energy2[:,T],1)
        #                        print('topol side : {} GPa'.format(-x1/x0))
        #                        self.pc2[T] = -x1/x0
                            #ftopo = interp1d(self.pressure[crit_index+1:], self.full_band_ren[crit_index+1:,T,b], kind='slinear',bounds_error=False, fill_value='extrapolate')
                            #self.extr_full_gap_ren_spline2[:,T,b] = ftopo(self.extr_pressure2)



#                print(ifile, T, gap_index)

#############################################
        # Start plotting
        for T in range(self.ntemp):

            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    #figure1
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    _arr[2][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.27,1.62), ncol = 1, fontsize=20)

                  #  self.set_legend_pgap2(_arr[2][0])
                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])

                    #figure2
                    _arr2[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    self.set_legend_pgap2(_arr2[2][0])
                    _arr2[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr2[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr2[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])

                else:
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=r'{} K'.format(self.ref_temp[T]))
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                  #  _arr[2][0].legend(numpoints=1, loc= 'center right', bbox_to_anchor=(1.27,1.62), ncol=1, fontsize=20)

                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', markersize=8, linewidth=1.5, color=self.color[T])

                    if self.extrapolate_ren:
                        # plot extrapolated data
                        _arr[0,0].plot(self.extr_pressure1, self.extr_full_gap_ren1[:,T,1], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure1, self.extr_full_gap_ren1[:,T,0], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[2,0].plot(self.extr_pressure1, self.extr_full_gap_ren1[:,T,1]-self.extr_full_gap_ren1[:,T,0], linestyle='--', linewidth=2.0, label = None, color = self.color[T])

                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T,1], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T,0], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[2,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T,1]-self.extr_full_gap_ren2[:,T,0], linestyle='--', linewidth=2.0, label = None, color = self.color[T])

                            ## This is just to test spline interpolation. It gives shit when extrapolating. 
#                        _arr[0,0].plot(self.extr_pressure1, self.extr_full_gap_ren_spline1[:,T,1], linestyle='dotted', linewidth=2.0, label = None, color = 'k')
#                        _arr[1,0].plot(self.extr_pressure1, self.extr_full_gap_ren_spline1[:,T,0], linestyle='dotted', linewidth=2.0, label = None, color = 'k')
#                        _arr[2,0].plot(self.extr_pressure1, self.extr_full_gap_ren_spline1[:,T,1]-self.extr_full_gap_ren_spline1[:,T,0], linestyle='--', linewidth=2.0, label = None, color = 'k')
#
#                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_ren_spline2[:,T,1], linestyle='dotted', linewidth=2.0, label = None, color = 'k')
#                        _arr[1,0].plot(self.extr_pressure2, self.extr_full_gap_ren_spline2[:,T,0], linestyle='dotted', linewidth=2.0, label = None, color = 'k')
#                        _arr[2,0].plot(self.extr_pressure2, self.extr_full_gap_ren_spline2[:,T,1]-self.extr_full_gap_ren_spline2[:,T,0], linestyle='--', linewidth=2.0, label = None, color = 'k')
 
 
            else:
                if self.split:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    _arr[1][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.02,1.0), ncol = 1, fontsize=26)


#                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)

                else:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)


########### Continue with final plotting options

        limms = [[-50.,50.],[-2.,50.],[-90.,40.]]
##########FIX ME : add default data and if condition
        limms = [self.cond_ylims, self.val_ylims, self.gap_ylims]

        for i in range(3):
            _arr[i,0].set_rasterization_zorder(0)
#            lims = _arr[i,0].get_ylim() 
            ylims = limms[i]    
            self.set_vrefs(_arr[i][0], self.pressure, 0.,style='dashed')
            self.set_hrefs(ylims, _arr[i][0], self.crit_pressure,'gray')
            self.set_hrefs(ylims, _arr[i][0], self.crit_pressure+0.2,'gray')
            _arr[i][0].fill([2.08,2.28,2.28,2.08],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2,zorder=-1, rasterized=True)
            _arr[i,0].set_ylim(limms[i])


            if self.split:
                self.set_vrefs(_arr2[i][0], self.pressure, 0.)
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'gray')
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure+0.2,'gray')
                lims = _arr[i,0].get_ylim() 
                _arr2[i][0].fill([2.08,2.28,2.28,2.08],[lims[0]-5,lims[0]-5,lims[1]+5,lims[1]+5],'gray',alpha=0.2,zorder=-1, rasterized=True)


        
        self.set_xaxis(_arr[2][0], self.pressure)
        self.set_yaxis_separate(_arr[0][0], 'CBM ren ({})'.format(self.gap_units),self.cond_ylims, self.cond_yticks)
        self.set_yaxis_separate(_arr[1][0], 'VBM ren ({})'.format(self.gap_units),self.val_ylims, self.val_yticks)
        self.set_yaxis_separate(_arr[2][0], 'Gap ren ({})'.format(self.gap_units),self.gap_ylims, self.gap_yticks)


        self.set_yaxis(_arr[0][0], 'CBM ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VBM ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))

        fig.subplots_adjust(left=0.11,bottom=0.08,right=0.81,top=0.95,wspace=0.2,hspace=0.12)

        #_arr[1,0].text(self.explicit_pressures[8],300, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=30)
        #_arr[1,0].text(self.explicit_pressures[35],300, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=30)

        #_arr[1,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=28, weight='bold', color='gray')



        if self.split:
            self.set_xaxis(_arr2[2][0], self.pressure)
            self.set_yaxis_separate(_arr2[0][0], 'CBM ren ({})'.format(self.gap_units),self.cond_ylims, self.cond_yticks)
            self.set_yaxis_separate(_arr2[1][0], 'VBM ren ({})'.format(self.gap_units),self.val_ylims, self.val_yticks)
            self.set_yaxis_separate(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units),self.gap_ylims, self.gap_yticks)

            self.set_title(_arr[0][0], self.title[0])
            self.set_title(_arr2[0][0], self.title[1])

        else:
            if self.main_title:
                self.set_title(_arr0[0][0], self.main_title)

        # Custum stuff
        _arr[2][0].text(0.7,-82, r'$\mathbb{Z}_2\!=\!0$', fontsize=32,color='k')
        _arr[2][0].text(2.7, -82, r'$\mathbb{Z}_2\!=\!1$',fontsize = 32,color='k')
        _arr[2][0].text(1.85,-60,r'$\Rightarrow$',fontsize=26,color='#5A5A5A')
        _arr[2][0].text(1.35,-60,r'WSM',fontsize=26,color='#5A5A5A',weight='bold')

        legend_handles = [] 
        for t, temp in enumerate(self.ref_temp):
            legend_handles.append(Line2D([0],[0],color=self.color[t],linewidth=2.5, label=r'{:>3.0f} K'.format(self.ref_temp[t])))
#                Line2D([0],[0],color='b',marker='o',markersize=8,linestyle='None',label=r'P$_{\text{C2}}$ plane')]
        #legend_handles.append('')
        #legend_handles.append('')

        legend1 = _arr[0][0].legend(handles=legend_handles, loc=9,bbox_to_anchor=(0.5,1.28),fontsize=24, handletextpad=0.4,handlelength=1.0,frameon=True,ncol = len(self.ref_temp),columnspacing=0.7)
        _arr[0][0].add_artist(legend1)
        
        legend2_handles=[]
        legend2_handles.append(Line2D([0],[0],marker='d',markersize=12,markerfacecolor='None', markeredgecolor='k', linestyle='None',label=r'P$_{\text{C1}}$ plane'))
        legend2_handles.append(Line2D([0],[0],marker='o',markersize=12,markerfacecolor='None', markeredgecolor='k', linestyle='None',label=r'P$_{\text{C2}}$ plane'))
        legend2 = _arr[0][0].legend(handles=legend2_handles, loc=1,bbox_to_anchor=(1.01,1.02),fontsize=24, handletextpad=0.4,handlelength=1.4,frameon=True,ncol = 1,labelspacing=0.1)
        _arr[0][0].add_artist(legend2)

        fig.subplots_adjust(hspace=0.0,top=0.90,right=0.95,left=0.12,bottom=0.09)

        if self.split:
            fig.align_ylabels()
            fig2.align_ylabels()
            self.save_figure_split(fig,fig2)
        else:
            fig.align_ylabels()
            self.save_figure(fig)

        #plt.show()

    def plot_splitted_subspaces(self):

        file_qty = len(self.zpr_fnames)

        # Define figure
        fig, _arr = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
        fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
        fig3, _arr3 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)

        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            raise Exception('Not implemented for split yet.')
#        if self.split:
#            fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
#            plt.subplots_adjust(hspace=0.05, top=0.95) 

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
               
#            # Set side of phase transition for split2 option
#            if self.split2:
#                if self.pressure[ifile] < self.crit_pressure:
#                    side=0
#                else:
#                    side=1

            # Retrieve data from file
#            if self.split:
#                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
##                self.gap_location = self.zpr.gap_location_split
#                self.gap_ren = self.zpr.indirect_gap_ren_split
#                self.band_ren = self.zpr.indirect_gap_ren_band_split
#
#                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
#                for a in range(2):
#                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
#                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])
#
#            if self.split2:
#                self.loc0 = self.zpr.unperturbed_indirect_gap_location_split
##                self.gap_location = self.zpr.indirect_gap_location_split
#                self.gap_ren = self.zpr.indirect_gap_ren_split
#                self.band_ren = self.zpr.indirect_gap_ren_band_split # T left/right vb/cb
#
#                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
#                for a in range(2):
#                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
#                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

 
#            else:
            self.loc0 = self.zpr.unperturbed_indirect_gap_location
            self.loc = self.zpr.indirect_gap_location

            self.gap_ren = self.zpr.indirect_gap_ren
            self.band_ren = self.zpr.indirect_gap_ren_band
            gap_index0 = np.zeros((2),dtype=int)
            gap_index = np.zeros((self.ntemp,2),dtype=int)
            for a in range(2):
                gap_index0[a] = self.find_gap_index(self.loc0[a])
                
                for t in range(self.ntemp):
                    gap_index[t,a] = self.find_gap_index(self.loc[t,a])


            # Initialize plotting arrays
            if ifile==0:
#                if self.split: # distinct data for left and rifht gap
#                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
#                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2,2)) # file, temp, vbcb, lr 
#                    self.ref_temp = self.temp
#                if self.split2: # single array, left gap for trovoal phase and right gap for topol phase
#                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
#                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2))
#                    self.ref_temp = self.temp
#                else:
                self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                self.full_band_ren = np.zeros((file_qty,self.ntemp,2))
                self.ref_temp = self.temp
                self.fan_occ = np.zeros((file_qty,self.ntemp,2))
                self.fan_unocc = np.zeros((file_qty,self.ntemp,2))
                self.ddw_occ = np.zeros((file_qty,self.ntemp,2))
                self.ddw_unocc = np.zeros((file_qty,self.ntemp,2))

            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

#            if self.split:
#                self.full_gap_ren[ifile,:,:] = self.gap_ren
#     #           print(np.shape(self.gap_ren))
#     #           print(np.shape(self.full_gap_ren))
#            if self.split2:
#                self.full_gap_ren[ifile,:] = self.gap_ren[:,side]
#            else:
            self.full_gap_ren[ifile,:] = self.gap_ren

            for T in range(self.ntemp):
                if self.unpert: ### FIX ME : WHAT DOES THIS DO???? This looks like the renorm(T) at the bare gap location...
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
#                    if self.split:
#                        for a in range(2):
#                            self.full_band_ren[ifile, T,0,a] = self.band_ren[T,a,0]  
#                            self.full_band_ren[ifile, T,1,a] = self.band_ren[T,a,1]
#
#                            diff = (self.full_band_ren[ifile,T,1,a]-self.full_band_ren[ifile,T,0,a]) - self.full_gap_ren[ifile,T,a]
#    #                        if abs(diff)>1E-6:
#    #                            print(ifile,T,a)
#    #                            print(self.full_gap_ren[ifile,T,a],self.full_band_ren[ifile,T,0,a],self.full_band_ren[ifile,T,1,a],self.full_band_ren[ifile,T,1,a]-self.full_band_ren[ifile,T,0,a])
#    #                            print(self.gap_location[T,a], self.kpoints[gap_index])
#                    if self.split2:
#                        self.full_band_ren[ifile, T, :] = self.band_ren[T,side,:]
#                    else:
                    for a in range(2): 
                        self.full_band_ren[ifile, T,a] = self.band_ren[T,a]
                        self.fan_occ[ifile,T,a] = self.zpr.fan_occ[0,gap_index[T,a],self.band_numbers[a]-1,T]*cst.ha_to_ev*1000
                        self.fan_unocc[ifile,T,a] = self.zpr.fan_unocc[0,gap_index[T,a],self.band_numbers[a]-1,T]*cst.ha_to_ev*1000
                        self.ddw_occ[ifile,T,a] = self.zpr.ddw_occ[0,gap_index[T,a],self.band_numbers[a]-1,T]*cst.ha_to_ev*1000
                        self.ddw_unocc[ifile,T,a] = self.zpr.ddw_unocc[0,gap_index[T,a],self.band_numbers[a]-1,T]*cst.ha_to_ev*1000

                    diff = (self.full_band_ren[ifile,T,1]-self.full_band_ren[ifile,T,0]) - self.full_gap_ren[ifile,T]
#                    if abs(diff)>1E-6:
    #                        print(ifile,T)
    #                        print(self.full_gap_ren[ifile,T],self.full_band_ren[ifile,T,0],self.full_band_ren[ifile,T,1],self.full_band_ren[ifile,T,1]-self.full_band_ren[ifile,T,0])
    #                        print(self.gap_location[T], self.kpoints[gap_index])

    #                print(ifile, T, gap_index)

            if self.verbose:

                print('Valence band')
                print(self.fan_occ[ifile,:,0])
                print(self.fan_unocc[ifile,:,0])
                print(self.ddw_occ[ifile,:,0])
                print(self.ddw_unocc[ifile,:,0])

                print('full band ren : {}'.format(self.full_band_ren[ifile,:,0]))
                print('sum : {}'.format(self.fan_occ[ifile,:,0]+self.fan_unocc[ifile,:,0]+self.ddw_occ[ifile,:,0]+self.ddw_unocc[ifile,:,0]))
                print('Conduction band')
                print(self.fan_occ[ifile,:,1])
                print(self.fan_unocc[ifile,:,1])
                print(self.ddw_occ[ifile,:,1])
                print(self.ddw_unocc[ifile,:,1])


                print('full band ren : {}'.format(self.full_band_ren[ifile,:,1]))
                print('sum : {}'.format(self.fan_occ[ifile,:,1]+self.fan_unocc[ifile,:,1]+self.ddw_occ[ifile,:,1]+self.ddw_unocc[ifile,:,1]))

        if self.crit_pressure is not None:
            crit_index = self.find_crit_index()
        else:
            crit_index = None
#        print(crit_index)

        temp_index = self.find_temp_index()
        print(temp_index)

        # Print pressure-dependent results, for chosen temperature
        print('Valence band')
        print('    Fan occ {}'.format(self.fan_occ[:,temp_index,0]))
        print('    Fan unocc {}'.format(self.fan_unocc[:,temp_index,0]))
        print('    DW occ {}'.format(self.ddw_occ[:,temp_index,0]))
        print('    DW unocc {}'.format(self.ddw_unocc[:,temp_index,0]))

        print('\n    Fan total {}'.format(self.fan_occ[:,temp_index,0],self.fan_unocc[:,temp_index,0]))
        print('    DW total {}'.format(self.ddw_occ[:,temp_index,0],self.ddw_unocc[:,temp_index,0]))
        print('\n    Occ total {}'.format(self.fan_occ[:,temp_index,0],self.ddw_occ[:,temp_index,0]))
        print('    Unocc total {}'.format(self.fan_unocc[:,temp_index,0],self.ddw_unocc[:,temp_index,0]))


        print('full band ren : {}'.format(self.full_band_ren[:,temp_index,0]))
        print('sum : {}'.format(self.fan_occ[:,temp_index,0]+self.fan_unocc[:,temp_index,0]+self.ddw_occ[:,temp_index,0]+self.ddw_unocc[:,temp_index,0]))
        print('\n\nConduction band')
        print('    Fan occ {}'.format(self.fan_occ[:,temp_index,1]))
        print('    Fan unocc {}'.format(self.fan_unocc[:,temp_index,1]))
        print('    DW occ {}'.format(self.ddw_occ[:,temp_index,1]))
        print('    DW unocc {}'.format(self.ddw_unocc[:,temp_index,1]))

        print('\n    Fan total {}'.format(self.fan_occ[:,temp_index,1],self.fan_unocc[:,temp_index,1]))
        print('    DW total {}'.format(self.ddw_occ[:,temp_index,1],self.ddw_unocc[:,temp_index,1]))
        print('\n    Occ total {}'.format(self.fan_occ[:,temp_index,1],self.ddw_occ[:,temp_index,1]))
        print('    Unocc total {}'.format(self.fan_unocc[:,temp_index,1],self.ddw_unocc[:,temp_index,1]))

        print('\n\nfull band ren : {}'.format(self.full_band_ren[:,temp_index,1]))
        print('sum : {}'.format(self.fan_occ[:,temp_index,1]+self.fan_unocc[:,temp_index,1]+self.ddw_occ[:,temp_index,1]+self.ddw_unocc[:,temp_index,1]))


 
#        for T in range(self.ntemp):
        for T in [temp_index]:

            if crit_index is not None:
                s = crit_index+1

#                if self.split:
#                    #figure1
#                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#
#                    _arr[2][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.27,1.62), ncol = 1, fontsize=20)
#
#                  #  self.set_legend_pgap2(_arr[2][0])
#                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
#                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
#                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
#
#                    #figure2
#                    _arr2[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr2[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr2[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
#
#                    self.set_legend_pgap2(_arr2[2][0])
#                    _arr2[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
#                    _arr2[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
#                    _arr2[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='o', markersize=8, linewidth=2.0, color=self.color[T])
#
#                else:
                _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')


                _arr[0][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='red', label='Fan occ')
                _arr[1][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='red', label='Fan occ')
                _arr[2][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,1]-self.fan_occ[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='red', label='Fan occ')

                _arr[0][0].plot(self.pressure[0:s], self.fan_unocc[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='blue', label='Fan unocc')
                _arr[1][0].plot(self.pressure[0:s], self.fan_unocc[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='blue', label='Fan unocc')
                _arr[2][0].plot(self.pressure[0:s], self.fan_unocc[0:s,T,1]-self.fan_unocc[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='blue', label='Fan unocc')

                _arr[0][0].plot(self.pressure[0:s], self.ddw_occ[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='green', label='DW occ')
                _arr[1][0].plot(self.pressure[0:s], self.ddw_occ[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='green', label='DW occ')
                _arr[2][0].plot(self.pressure[0:s], self.ddw_occ[0:s,T,1]-self.ddw_occ[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='green', label='DW occ')

                _arr[0][0].plot(self.pressure[0:s], self.ddw_unocc[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='orange', label='DW unocc')
                _arr[1][0].plot(self.pressure[0:s], self.ddw_unocc[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='orange', label='DW unocc')
                _arr[2][0].plot(self.pressure[0:s], self.ddw_unocc[0:s,T,1]-self.ddw_unocc[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='orange', label='DW unocc')



                _arr[2][0].legend(numpoints=1, loc= 'center right', bbox_to_anchor=(1.27,1.62), ncol=1, fontsize=20)

                # Fig 2 : Fan and DW total
                _arr2[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr2[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr2[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')

                _arr2[0][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,1]+self.fan_unocc[0:s,T,1] , marker='o', markersize=8, linewidth=2.0, color='red', label='Fan')
                _arr2[1][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,0]+self.fan_unocc[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='red', label='Fan')
                _arr2[2][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,1]+self.fan_unocc[0:s,T,1]-(self.fan_occ[0:s,T,0]+self.fan_unocc[0:s,T,0]), marker='o', markersize=8, linewidth=2.0, color='red', label='Fan')

                _arr2[0][0].plot(self.pressure[0:s], self.ddw_occ[0:s,T,1]+self.ddw_unocc[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='blue', label='DW')
                _arr2[1][0].plot(self.pressure[0:s], self.ddw_occ[0:s,T,0]+self.ddw_unocc[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='blue', label='DW')
                _arr2[2][0].plot(self.pressure[0:s], self.ddw_occ[0:s,T,1]+self.ddw_unocc[0:s,T,1]-(self.ddw_occ[0:s,T,0]+self.ddw_unocc[0:s,T,0]), marker='o', markersize=8, linewidth=2.0, color='blue', label='DW')

                _arr2[2][0].legend(numpoints=1, loc= 'center right', bbox_to_anchor=(1.17,1.62), ncol=1, fontsize=20)

                # Fig 3 : Total contribution from occupied/unoccupied bands
                _arr3[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr3[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr3[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')

                _arr3[0][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,1]+self.ddw_occ[0:s,T,1] , marker='o', markersize=8, linewidth=2.0, color='red', label='occ')
                _arr3[1][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,0]+self.ddw_occ[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='red', label='occ')
                _arr3[2][0].plot(self.pressure[0:s], self.fan_occ[0:s,T,1]+self.ddw_occ[0:s,T,1]-(self.fan_occ[0:s,T,0]+self.ddw_occ[0:s,T,0]), marker='o', markersize=8, linewidth=2.0, color='red', label='occ')

                _arr3[0][0].plot(self.pressure[0:s], self.fan_unocc[0:s,T,1]+self.ddw_unocc[0:s,T,1], marker='o', markersize=8, linewidth=2.0, color='blue', label='unocc')
                _arr3[1][0].plot(self.pressure[0:s], self.fan_unocc[0:s,T,0]+self.ddw_unocc[0:s,T,0], marker='o', markersize=8, linewidth=2.0, color='blue', label='unocc')
                _arr3[2][0].plot(self.pressure[0:s], self.fan_unocc[0:s,T,1]+self.ddw_unocc[0:s,T,1]-(self.fan_unocc[0:s,T,0]+self.ddw_unocc[0:s,T,0]), marker='o', markersize=8, linewidth=2.0, color='blue', label='unocc')

                _arr3[2][0].legend(numpoints=1, loc= 'center right', bbox_to_anchor=(1.17,1.62), ncol=1, fontsize=20)


            ####################################
            ####################################
                    ### After Pc

                _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')


                _arr[0][0].plot(self.pressure[s:], self.fan_occ[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='red', label='Fan occ')
                _arr[1][0].plot(self.pressure[s:], self.fan_occ[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='red', label='Fan occ')
                _arr[2][0].plot(self.pressure[s:], self.fan_occ[s:,T,1]-self.fan_occ[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='red', label='Fan occ')

                _arr[0][0].plot(self.pressure[s:], self.fan_unocc[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='blue', label='Fan unocc')
                _arr[1][0].plot(self.pressure[s:], self.fan_unocc[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='blue', label='Fan unocc')
                _arr[2][0].plot(self.pressure[s:], self.fan_unocc[s:,T,1]-self.fan_unocc[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='blue', label='Fan unocc')

                _arr[0][0].plot(self.pressure[s:], self.ddw_occ[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='green', label='DW occ')
                _arr[1][0].plot(self.pressure[s:], self.ddw_occ[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='green', label='DW occ')
                _arr[2][0].plot(self.pressure[s:], self.ddw_occ[s:,T,1]-self.ddw_occ[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='green', label='DW occ')

                _arr[0][0].plot(self.pressure[s:], self.ddw_unocc[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='orange', label='DW unocc')
                _arr[1][0].plot(self.pressure[s:], self.ddw_unocc[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='orange', label='DW unocc')
                _arr[2][0].plot(self.pressure[s:], self.ddw_unocc[s:,T,1]-self.ddw_unocc[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='orange', label='DW unocc')


                # Fan and DW total
                _arr2[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr2[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr2[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')

                _arr2[0][0].plot(self.pressure[s:], self.fan_occ[s:,T,1]+self.fan_unocc[s:,T,1] , marker='o', markersize=8, linewidth=2.0, color='red', label='Fan')
                _arr2[1][0].plot(self.pressure[s:], self.fan_occ[s:,T,0]+self.fan_unocc[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='red', label='Fan')
                _arr2[2][0].plot(self.pressure[s:], self.fan_occ[s:,T,1]+self.fan_unocc[s:,T,1]-(self.fan_occ[s:,T,0]+self.fan_unocc[s:,T,0]), marker='o', markersize=8, linewidth=2.0, color='red', label='Fan')

                _arr2[0][0].plot(self.pressure[s:], self.ddw_occ[s:,T,1]+self.ddw_unocc[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='blue', label='DW')
                _arr2[1][0].plot(self.pressure[s:], self.ddw_occ[s:,T,0]+self.ddw_unocc[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='blue', label='DW')
                _arr2[2][0].plot(self.pressure[s:], self.ddw_occ[s:,T,1]+self.ddw_unocc[s:,T,1]-(self.ddw_occ[s:,T,0]+self.ddw_unocc[s:,T,0]), marker='o', markersize=8, linewidth=2.0, color='blue', label='DW')

                # Fig 3 : Total contribution from occupied/unoccupied bands
                _arr3[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr3[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')
                _arr3[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', markersize=8, linewidth=2.0, color='black', label='Full')

                _arr3[0][0].plot(self.pressure[s:], self.fan_occ[s:,T,1]+self.ddw_occ[s:,T,1] , marker='o', markersize=8, linewidth=2.0, color='red', label='occ')
                _arr3[1][0].plot(self.pressure[s:], self.fan_occ[s:,T,0]+self.ddw_occ[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='red', label='occ')
                _arr3[2][0].plot(self.pressure[s:], self.fan_occ[s:,T,1]+self.ddw_occ[s:,T,1]-(self.fan_occ[s:,T,0]+self.ddw_occ[s:,T,0]), marker='o', markersize=8, linewidth=2.0, color='red', label='occ')

                _arr3[0][0].plot(self.pressure[s:], self.fan_unocc[s:,T,1]+self.ddw_unocc[s:,T,1], marker='o', markersize=8, linewidth=2.0, color='blue', label='unocc')
                _arr3[1][0].plot(self.pressure[s:], self.fan_unocc[s:,T,0]+self.ddw_unocc[s:,T,0], marker='o', markersize=8, linewidth=2.0, color='blue', label='unocc')
                _arr3[2][0].plot(self.pressure[s:], self.fan_unocc[s:,T,1]+self.ddw_unocc[s:,T,1]-(self.fan_unocc[s:,T,0]+self.ddw_unocc[s:,T,0]), marker='o', markersize=8, linewidth=2.0, color='blue', label='unocc')

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


        # set vrefs at 0
        for i in range(3):
            self.set_vrefs(_arr[i][0], self.pressure, 0.)
            self.set_vrefs(_arr2[i][0], self.pressure, 0.)
            self.set_vrefs(_arr3[i][0], self.pressure, 0.)

#        limms = [[-50.,50.],[-2.,50.],[-90.,40.]]
##########FIX ME : add default data and if condition
        if self.cond_ylims and self.val_ylims:
            if self.gap_ylims:
                limms = [self.cond_ylims, self.val_ylims, self.gap_ylims]
                for i in range(3):
                    lims = _arr[i,0].get_ylim() 
                    ylims = limms[i]    
                    self.set_hrefs(ylims, _arr[i][0], self.crit_pressure,'black')
                    self.set_hrefs(ylims, _arr[i][0], self.crit_pressure+0.2,'black')
                    _arr[i][0].fill([2.08,2.28,2.28,2.08],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
                    _arr[i,0].set_ylim(limms[i])


#            if self.split:
#                self.set_vrefs(_arr2[i][0], self.pressure, 0.)
#                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'black')
#                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure+0.2,'black')
#                lims = _arr[i,0].get_ylim() 
#                _arr2[i][0].fill([2.08,2.28,2.28,2.08],[lims[0]-5,lims[0]-5,lims[1]+5,lims[1]+5],'gray',alpha=0.2)


        
        self.set_xaxis(_arr[2][0], self.pressure)
        self.set_yaxis(_arr[0][0], 'CB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))

        self.set_xaxis(_arr2[2][0], self.pressure)
        self.set_yaxis(_arr2[0][0], 'CB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr2[1][0], 'VB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units))

        self.set_xaxis(_arr3[2][0], self.pressure)
        self.set_yaxis(_arr3[0][0], 'CB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr3[1][0], 'VB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr3[2][0], 'Gap ren ({})'.format(self.gap_units))

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
            if self.main_title:
                self.set_title(_arr[0][0], self.main_title)


        if self.split:
            self.save_figure_split(fig,fig2)
        else:
            self.save_figure_subspaces(fig,fig2,fig3)

        plt.show()

    def plot_te_pgap_separate_indirect(self):

        file_qty = len(self.te_fnames)

        '''FIX ME : add extrapolate as an external option'''
        extrapolate = False
        # Define figure
        fig, _arr = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            raise Exception('split not yet implemented for te_pgap')
            fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
            plt.subplots_adjust(hspace=0.05, top=0.95) 

        if not self.band_numbers:
            raise Exception('Must provide valence and conduction band numbers vis band_numbers')


        # Read and treat all input files
        for ifile, te_file in enumerate(self.te_fnames):
                
            print(te_file)
            # Define file class
            self.te = TEfile(te_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.te.nsppol
            self.nkpt = self.te.nkpt
            self.max_band = self.te.max_band
            self.ntemp = self.te.ntemp
            self.temp = self.te.temp
            self.kpoints = self.te.kpoints
           ################### HERE ############################ 
#            self.gridsize = self.zpr.gridsize
            self.gap_units = self.te.gap_energy_units

            input_units = self.te.gap_energy_units

            if input_units == self.units:
                pass
            else:
                if self.units == 'eV':
                    self.eig0 = self.te.eig0*cst.ha_to_ev
                    self.eigcorr = self.te.eigcorr*cst.ha_to_ev
                if self.units == 'meV':
                    self.eig0 = self.te.eig0*cst.ha_to_ev*1000
                    self.eigcorr = self.te.eigcorr*cst.ha_to_ev*1000

            # Set side of phase transition for split2 option
            if self.split2:
                if self.pressure[ifile] < self.crit_pressure:
                    side=0
                else:
                    side=1

            # Retrieve data from file
            if self.split:
                self.loc0 = self.te.unperturbed_indirect_gap_location_split
#                self.gap_location = self.te.gap_location_split
                self.gap_ren = self.te.indirect_gap_ren_split
                self.band_ren = self.te.indirect_gap_ren_band_split

                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                for a in range(2):
                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

            if self.split2:
                self.loc0 = self.te.unperturbed_indirect_gap_location_split
#                self.gap_location = self.te.indirect_gap_location_split
                self.gap_ren = self.te.indirect_gap_ren_split
                self.band_ren = self.te.indirect_gap_ren_band_split # T left/right vb/cb

                gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                for a in range(2):
                    gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                    gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

 
            else:
                self.loc0 = self.te.unperturbed_indirect_gap_location
#                self.gap_location = self.te.gap_location
                self.gap_ren = self.te.indirect_gap_ren
                self.band_ren = self.te.indirect_gap_ren_band
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
                        if abs(diff)>1E-6:
                            print('Coherence warning!!! for the following file/temperature, the difference of band ren is not the same as gap ren') 
                            print(ifile,T)
                            print(self.full_gap_ren[ifile,T],self.full_band_ren[ifile,T,0],self.full_band_ren[ifile,T,1],self.full_band_ren[ifile,T,1]-self.full_band_ren[ifile,T,0])
    #                        print(self.gap_location[T], self.kpoints[gap_index])

    #                print(ifile, T, gap_index)

    ##### Finished reading TE data #####

    ### From here, start plotting ####

        if self.crit_pressure is not None:
            crit_index = self.find_pressure_index(self.crit_pressure)
        else:
            crit_index = None

        print(crit_index)
        for T in range(self.ntemp):

            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    #figure1
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    _arr[2][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.27,1.62), ncol = 1, fontsize=20)

                  #  self.set_legend_pgap2(_arr[2][0])
                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])

                    #figure2
                    _arr2[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    self.set_legend_pgap2(_arr2[2][0])
                    _arr2[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr2[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr2[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])

                else:
                    _arr[0][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=r'{} K'.format(self.ref_temp[T]))
                    _arr[1][0].plot(self.pressure[0:s], self.full_band_ren[0:s,T,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                  #  _arr[2][0].legend(numpoints=1, loc= 'center right', bbox_to_anchor=(1.27,1.62), ncol=1, fontsize=20)

                    _arr[0][0].plot(self.pressure[s:], self.full_band_ren[s:,T,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_band_ren[s:,T,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', markersize=8, linewidth=1.5, color=self.color[T])


            else:
                if self.split:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    _arr[1][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.02,1.0), ncol = 1, fontsize=26)


#                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)

                else:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)


        limms = [[-50.,50.],[-2.,50.],[-90.,40.]]
##########FIX ME : add default data and if condition
        limms = [self.cond_ylims, self.val_ylims, self.gap_ylims]
        for i in range(3):
            _arr[i,0].set_rasterization_zorder(0)
#            lims = _arr[i,0].get_ylim() 
            ylims = limms[i]    
            self.set_vrefs(_arr[i][0], self.pressure, 0.,style='dashed')
            self.set_hrefs(ylims, _arr[i][0], self.crit_pressure,'gray')
            self.set_hrefs(ylims, _arr[i][0], self.crit_pressure+0.2,'gray')
            _arr[i][0].fill([2.08,2.28,2.28,2.08],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2,zorder=-1,rasterized=True)
            _arr[i,0].set_ylim(limms[i])


            if self.split:
                self.set_vrefs(_arr2[i][0], self.pressure, 0.)
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'gray')
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure+0.2,'gray')
                lims = _arr[i,0].get_ylim() 
                _arr2[i][0].fill([2.08,2.28,2.28,2.08],[lims[0]-5,lims[0]-5,lims[1]+5,lims[1]+5],'gray',alpha=0.2,zorder=-1,rasterized=True)


        
        self.set_xaxis(_arr[2][0], self.pressure)
        self.set_yaxis_separate(_arr[0][0], 'CBM ren ({})'.format(self.gap_units),self.cond_ylims, self.cond_yticks)
        self.set_yaxis_separate(_arr[1][0], 'VBM ren ({})'.format(self.gap_units),self.val_ylims, self.val_yticks)
        self.set_yaxis_separate(_arr[2][0], 'Gap ren ({})'.format(self.gap_units),self.gap_ylims, self.gap_yticks)


        self.set_yaxis(_arr[0][0], 'CBM ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VBM ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))

        # Custum stuff
        _arr[2][0].text(0.7,-48, r'$\mathbb{Z}_2\!=\!0$', fontsize=32,color='k')
        _arr[2][0].text(3.8, -48, r'$\mathbb{Z}_2\!=\!1$',fontsize = 32,color='k')
        _arr[2][0].text(1.85,-25,r'$\Rightarrow$',fontsize=26,color='#5A5A5A')
        _arr[2][0].text(1.35,-25,r'WSM',fontsize=26,color='#5A5A5A',weight='bold')

#        if self.main_title:
#            self.set_title(_arr[0][0], self.main_title)

        fig.subplots_adjust(left=0.11,bottom=0.08,right=0.81,top=0.95,wspace=0.2,hspace=0.12)

        #_arr[1,0].text(self.explicit_pressures[8],300, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=30)
        #_arr[1,0].text(self.explicit_pressures[35],300, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=30)

        #_arr[1,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=28, weight='bold', color='gray')



#        if self.split:
#            self.set_xaxis(_arr2[2][0], self.pressure)
#            self.set_yaxis_separate(_arr2[0][0], 'CB ren ({})'.format(self.gap_units),self.cond_ylims, self.cond_yticks)
#            self.set_yaxis_separate(_arr2[1][0], 'VB ren ({})'.format(self.gap_units),self.val_ylims, self.val_yticks)
#            self.set_yaxis_separate(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units),self.gap_ylims, self.gap_yticks)
#
#            self.set_title(_arr[0][0], self.title[0])
#            self.set_title(_arr2[0][0], self.title[1])
#
#        else:
#            if self.main_title:
#                self.set_title(_arr0[0][0], self.main_title)
#
#        # Custum stuff
#        _arr[2][0].text(0.7,-82, r'$\mathbb{Z}_2\!=\!0$', fontsize=24,color='k')
#        _arr[2][0].text(2.7, -82, r'$\mathbb{Z}_2\!=\!1$',fontsize = 24,color='k')
#        _arr[2][0].text(1.85,-60,r'$\Rightarrow$',fontsize=20,color='#5A5A5A')
#        _arr[2][0].text(1.35,-60,r'WSM',fontsize=20,color='#5A5A5A',weight='bold')
#
        legend_handles = [] 
        for t, temp in enumerate(self.ref_temp):
            legend_handles.append(Line2D([0],[0],color=self.color[t],linewidth=2.5, label=r'{:>3.0f} K'.format(self.ref_temp[t])))
#                Line2D([0],[0],color='b',marker='o',markersize=8,linestyle='None',label=r'P$_{\text{C2}}$ plane')]
        #legend_handles.append('')
        #legend_handles.append('')

        legend1 = _arr[0][0].legend(handles=legend_handles, loc=9,bbox_to_anchor=(0.5,1.28),fontsize=24, handletextpad=0.4,handlelength=1.0,frameon=True,ncol = len(self.ref_temp),columnspacing=0.7)
        _arr[0][0].add_artist(legend1)
        
        legend2_handles=[]
        legend2_handles.append(Line2D([0],[0],marker='d',markersize=12,markerfacecolor='None', markeredgecolor='k', linestyle='None',label=r'P$_{\text{C1}}$ plane'))
        legend2_handles.append(Line2D([0],[0],marker='o',markersize=12,markerfacecolor='None', markeredgecolor='k', linestyle='None',label=r'P$_{\text{C2}}$ plane'))
        legend2 = _arr[1][0].legend(handles=legend2_handles, loc=1,bbox_to_anchor=(1.01,0.38),fontsize=24, handletextpad=0.4,handlelength=1.4,frameon=True,ncol = 1,labelspacing=0.1)

        _arr[1][0].add_artist(legend2)

        fig.subplots_adjust(hspace=0.0,bottom=0.09,top=0.90,right=0.95,left=0.14)

        if self.split:
            fig.align_ylabels()
            fig2.align_ylabels()
            self.save_figure_split(fig,fig2)
        else:
            fig.align_ylabels()
            self.save_figure(fig)

        #plt.show()

    def plot_te_pgap_indirect(self):


        file_qty = len(self.te_fnames)

        # Define figure
        '''FIX ME : add only as an external option'''
        only=True  # only the Egap(P) for each T. If False, add gap renorm on top subplot
        extrapolate_bands = True ## Make this the default when everything is done

        if only:
            fig, _arr = plt.subplots(1,1, figsize=self.figsize, squeeze=False, sharex=True)
        else:
            fig, _arr = plt.subplots(2,1, figsize=self.figsize, squeeze=False, sharex=True)
        #plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            raise Exception('split not yet implemented for te_pgap')
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
        for ifile, te_file in enumerate(self.te_fnames):
                
            # Define file class
            self.te = TEfile(te_file, read=False)

            # Read input file
            self.read_file()

            # Set parameters for this file
            self.nsppol = self.te.nsppol
            self.nkpt = self.te.nkpt
            self.max_band = self.te.max_band
            self.ntemp = self.te.ntemp
            self.temp = self.te.temp
            self.kpoints = self.te.kpoints
            
            #self.gridsize = self.te.gridsize
            self.gap_units = self.te.gap_energy_units

            input_units = self.te.gap_energy_units
            
            if input_units == self.units:
                pass
            else:
                if self.units == 'eV':
                    self.eig0 = self.te.eig0*cst.ha_to_ev
                    self.eigcorr = self.te.eigcorr*cst.ha_to_ev
                if self.units == 'meV':
                    self.eig0 = self.te.eig0*cst.ha_to_ev*1000
                    self.eigcorr = self.te.eigcorr*cst.ha_to_ev*1000


            if self.split :
                self.energy0 = self.te.unperturbed_indirect_gap_energy_split
                self.loc0 = self.te.unperturbed_indirect_gap_location_split
                self.gap_location = self.te.indirect_gap_location_split
                self.gap_energy = self.te.indirect_gap_energy_split
                self.gap_ren = self.te.indirect_gap_ren_split
                self.band_energy = self.te.indirect_gap_energy_band_split

                if self.unpert:
                    gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                    for a in range(2):
                        gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                        gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

            elif self.split2:
                if self.pressure[ifile] > self.crit_pressure:
                    side = 1    
                else:
                    side = 0     

                self.energy0 = self.te.unperturbed_indirect_gap_energy_split
                self.loc0 = self.te.unperturbed_indirect_gap_location_split
                self.gap_location = self.te.indirect_gap_location_split
                self.gap_energy = self.te.indirect_gap_energy_split
                self.gap_ren = self.te.indirect_gap_ren_split
                self.band_energy = self.te.indirect_gap_energy_band_split

                if self.unpert:
                    gap_index0 = np.zeros((2,2),dtype = int) # (vbcb, lr)
                    for a in range(2):
                        gap_index0[0,a] = self.find_gap_index(self.loc0[a,0])
                        gap_index0[1,a] = self.find_gap_index(self.loc0[a,1])

 

            else:
                self.energy0 = self.te.unperturbed_indirect_gap_energy
                self.loc0 = self.te.unperturbed_indirect_gap_location
                self.gap_location = self.te.indirect_gap_location
                self.gap_energy = self.te.indirect_gap_energy
                self.gap_ren = self.te.indirect_gap_ren
                self.band_ren = self.te.indirect_gap_ren_band

                if self.unpert:
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
                    self.full_band_ren = np.zeros((file_qty, self.ntemp,2))

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
                self.full_band_ren[ifile,:,:] = self.band_ren
                if not self.gap_fname:
                    self.full_energy0[ifile] = self.energy0
                else:
                    p_index = self.find_pressure_index(self.pressure[ifile])
                    #print(p_index)
                    self.full_energy0[ifile] = self.explicit_gap_energies[p_index]

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

        ### End of loop on files. All data is stored now.

        if self.crit_pressure is not None:
            crit_index = self.find_crit_pressure(self.crit_pressure,self.pressure)
            if self.gap_fname is not None:
                crit_index2 = self.find_crit_pressure(self.crit_pressure2,self.explicit_pressures)
                print('crit_index2',crit_index2)
        else:
            crit_index = None


        if self.follow:
            self.full_energy0[-1,:] = -self.full_energy0[-1,:]
            self.full_gap_energy[-1,:,:] = -self.full_gap_energy[-1,:,:]

    # Extrapolate behavior linearly towards TPT... it IS a bit sketchy.
        pressure1 = np.arange(1.5,2.0,0.1)
        tmparr = np.arange(2.1,2.6,0.1) # change 2.2 to 2.6 foir full BZ
        pressure2 = np.arange(2.3,3.0,0.1)
        tmparr2 = np.arange(2.0,2.3,0.1)
#        self.extr_pressure2 = np.arange(2.3, 3.0, 0.1)

        # linearize unperturbed gap energy to extrapolate inside WSM phase
        y0,y1 = np.polyfit(self.explicit_pressures[15:20], self.explicit_gap_energies[15:20],1)
        tmp2 = y1 + tmparr*y0
        self.extr_pressure1 = np.concatenate((pressure1, tmparr))

        y0,y1 = np.polyfit(self.explicit_pressures[21:29], self.explicit_gap_energies[21:29],1)
        tmp22 = y1 + tmparr2*y0
        self.extr_pressure2 = np.concatenate((tmparr2,pressure2))

        self.extr_full_gap_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp))
        self.extr_full_gap_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp))
        self.extr_full_gap_energy1 = np.zeros((len(self.extr_pressure1), self.ntemp))
        self.extr_full_gap_energy2 = np.zeros((len(self.extr_pressure2), self.ntemp))
        self.extr_full_band_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp,2)) #pressure, temperature, val/cond
        self.extr_full_band_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp,2)) 
        extr_energy01 = np.zeros((len(pressure1)))
        extr_energy02 = np.zeros((len(pressure2)))
#        self.extr_full_energy02 = np.zeros((len(self.extr_pressure2)))

        # Get full_gap_energy0 (explicit bare gap energies) from explicit data
        for p, pres in enumerate(pressure1):
            ind = self.find_pressure_index(pres)
            extr_energy01[p] = self.explicit_gap_energies[ind]

        self.extr_full_energy01 = np.concatenate((extr_energy01,tmp2))

        for p, pres in enumerate(pressure2):
            ind = self.find_pressure_index(pres)
            extr_energy02[p] = self.explicit_gap_energies[ind]

        self.extr_full_energy02 = np.concatenate((tmp22,extr_energy02))

        #Store extrapolated phase boundaries (T)
        self.pc1 = np.zeros((self.ntemp))
        self.pc2 = np.zeros((self.ntemp))

        for T in range(self.ntemp):

            if self.split:
                if extrapolate_bands:
                    raise Exception('Extrapolate_bands not implemented for split yet')
                else:
                    #trivial side
                    x0,x1 = np.polyfit(self.pressure[0:crit_index+1], self.full_gap_ren[0:crit_index+1,T,0], 1)
                    self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]

                    # topol side
                    x0,x1 =  np.polyfit(self.pressure[crit_index+1:], self.full_gap_ren[crit_index+1:,T,0], 1)
                    self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]

            elif self.split2:
                if extrapolate_bands:
                    raise Exception('Extrapolate_bands not implemented for split yet')
                else:
                    #trivial side
                    x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_gap_ren[crit_index-1:crit_index+1,T], 1)
                    self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]

                    # topol side
                    x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_gap_ren[crit_index+1:crit_index+3,T], 1)
                    self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]


            else:
                if extrapolate_bands:
                    #loop on band extrema
                    print('\n\nFor T={} K:'.format(self.temp[T]))

                    for b in range(2):
                        #trivial side
                        x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_band_ren[crit_index-1:crit_index+1,T,b], 1)
                        self.extr_full_band_ren1[:,T,b] = x1 + x0*self.extr_pressure1

                        # topol side
                        x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_band_ren[crit_index+1:crit_index+3,T,b], 1)
                        self.extr_full_band_ren2[:,T,b] = x1 + x0*self.extr_pressure2
                        
                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_band_ren1[:,T,1] - self.extr_full_band_ren1[:,T,0]
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_band_ren2[:,T,1] - self.extr_full_band_ren2[:,T,0]


#                    print('For T={} K:'.format(self.temp[T]))
                    x0,x1 = np.polyfit(self.extr_pressure1,self.extr_full_gap_energy1[:,T],1)
                    print('trivial side : {} GPa'.format(-x1/x0))
                    self.pc1[T] = -x1/x0
                    x0,x1 = np.polyfit(self.extr_pressure2,self.extr_full_gap_energy2[:,T],1)
                    print('topol side : {} GPa'.format(-x1/x0))
                    self.pc2[T] = -x1/x0

                else:
                    #trivial side
                    x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_gap_ren[crit_index-1:crit_index+1,T], 1)
                    print('For T={} K:'.format(self.temp[T]))
                    self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]
                    x0,x1 = np.polyfit(self.extr_pressure1,self.extr_full_gap_energy1[:,T],1)
                    print('trivial side : {} GPa'.format(-x1/x0))
                    self.pc1[T] = -x1/x0
                    # topol side
                    x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_gap_ren[crit_index+1:crit_index+3,T], 1)
                    self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]
                    x0,x1 = np.polyfit(self.extr_pressure2,self.extr_full_gap_energy2[:,T],1)
                    print('topol side : {} GPa'.format(-x1/x0))
                    self.pc2[T] = -x1/x0

        #print(self.extr_pressure1)
        #print(self.extr_full_gap_energy1)
        #print(self.extr_pressure2)
        #print('before plotting: extr_full_gap_energy2',self.extr_full_gap_energy2)

        for T in range(self.ntemp):
            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    # figure 1
                    if T==self.temp[0]:
                        ###### FIX ME : keep explicit_pressures and explicit_gap_energies, just get a new crit index for this array! 
                        if self.gap_fname is not None:
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='Static T=0')
                            if self.crit_index2 is not None:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='Static T=0')

                            else:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='Static T=0')

                        else:
                            if only:
                                _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='Static T=0')
                            else:
                                _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='Static T=0')
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
                                    _arr[0][0].plot(self.explicit_pressures[crit_index2+2,:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[crit_index2+2:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label='Static T=0')

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
                                _arr2[0][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                            else:
                                arr2[1][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                        else:
                            if only:
                                _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                            else:
                                _arr2[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
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

#########################################
       #-------- crit_index but no split or split2
                else:
                    if T==self.temp[0]:
                        ###### FIX ME : keep explicit_pressures and explicit_gap_energies, just get a new crit index for this array! 
                        if self.gap_fname is not None:
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='Static T=0')
                            if crit_index2 is not None:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[:crit_index2+1], self.explicit_gap_energies[:crit_index2+1],marker='d',markersize=7, color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:crit_index2+1], self.explicit_gap_energies[:crit_index2+1],marker='d',markersize=7, color='black', label='Static T=0')
                            else:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', markersize=7,color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d',markersize=7, color='black', label='Static T=0')

                        else:
                            if only:
                                _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black', label='Static T=0')
                            else:
                                _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black', label='Static T=0')
#                    else:
#                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black')
    
                    if only:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='d',markersize=7, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    else:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='d',markersize=7, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
        #                self.set_legend_gap(_arr[0][0])
                        _arr[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='d',markersize=7, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
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
                                     _arr[0][0].plot(self.explicit_pressures[crit_index2+1:], self.explicit_gap_energies[crit_index2+1:],marker='o',markersize=7, color='black', label=None)
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[crit_index2+1:], self.explicit_gap_energies[crit_index2+1:],marker='o',markersize=7, color='black', label=None)

                    if only:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o',markersize=7, linewidth=2.0, color=self.color[T])
                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])

                    else:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o',markersize=7, linewidth=2.0, color=self.color[T])
                        _arr[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o',markersize=7, linewidth=2.0, color=self.color[T])

                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


                    if self.gap_fname is not None:
                        pass
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
                    raise Exception('Only keyword not implemented yet without critical pressure index')
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
        if self.gap_ylims is not None:
            limms[0] = self.gap_ylims
        if self.egap_ylims is not None:
            limms[1] = self.egap_ylims

        if only:
            ylims = limms[1]    
            self.set_vrefs(_arr[0][0], self.pressure, 0.)
#            if i==0:
            if self.crit_pressure:
                pc1 = self.crit_pressure
                self.set_hrefs(ylims, _arr[0][0], pc1,'gray')

                if self.crit_pressure2:
                    pc2 = self.crit_pressure2
                    self.set_hrefs(ylims, _arr[0][0], pc2,'gray')
                    _arr[0][0].fill([pc1,pc2,pc2,pc1],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
                else:
                    print('No value was given for PC2')
            else:
                print('No value was given for PC1')
#            self.set_hrefs(ylims, _arr[0,0], 2.18,'orange')
#            self.set_hrefs(ylims, _arr[0,0], 2.50,'orange')
#            _arr[0,0].fill([2.18,2.50,2.50,2.18],[ylims[0],ylims[0],ylims[1],ylims[1]],'orange',alpha=0.2)
#            self.set_hrefs(ylims, _arr[0,0], 2.68,'magenta')


        else:
            for i in range(2):
    #            lims = _arr[i,0].get_ylim() 
    #            ylims = lims
                ylims = limms[i]    
                self.set_vrefs(_arr[i][0], self.pressure, 0.)
    #            if i==0:
                pc1 = self.crit_pressure
                pc2 = self.crit_pressure2

                self.set_hrefs(ylims, _arr[i][0], pc1,'black')
                self.set_hrefs(ylims, _arr[i][0], pc2,'black')
                _arr[i][0].fill([pc1,pc2,pc2,pc1],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
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


#        if only:
#            _arr[0,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=24) 
#            _arr[0,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=24) 
#
##            _arr[0,0].text(self.explicit_pressures[5],290, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=26)
##            _arr[0,0].text(self.explicit_pressures[38],290, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=26)
#
#
#            _arr[0,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=26, weight='bold', color='gray')
#            _arr[0,0].text(self.explicit_pressures[25],250, r'WSM 300K', fontsize=26, weight='bold', color='orange')
#
#        else:
#            _arr[1,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=24) 
#            _arr[1,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=24) 
#
##            _arr[1,0].text(self.explicit_pressures[5],290, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=26)
##            _arr[1,0].text(self.explicit_pressures[38],290, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=26)
#
#
#            _arr[1,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=26, weight='bold', color='gray')
#            _arr[1,0].text(self.explicit_pressures[25],250, r'WSM 300K', fontsize=26, weight='bold', color='orange')



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

        # Custom text
        if only:
            _arr[0][0].text(0.7,250, r'$\mathbb{Z}_2\!=\!0$',fontsize=24,color='k')
            _arr[0][0].text(2.7,250, r'$\mathbb{Z}_2\!=\!1$',fontsize=24,color='k')
            _arr[0][0].text(1.90,200,r'$\Rightarrow$',fontsize=20,color='#5A5A5A')
            _arr[0][0].text(1.50,210,r'Static',fontsize=20,color='#5A5A5A',weight='bold')
            _arr[0][0].text(1.50,190,r'WSM',fontsize=20,color='#5A5A5A',weight='bold')

            legend_handles = []
            legend_handles.append(Line2D([0],[0],color='k',linewidth=1.5,label=r'0 K Static'))
            for t,T in enumerate(self.ref_temp):
                legend_handles.append(Line2D([0],[0],color=self.color[t],linewidth=1.5,label=r'{:>3.0f} K'.format(T)))

            legend1 = _arr[0][0].legend(handles=legend_handles, loc=9, bbox_to_anchor=(0.5,1.15),fontsize=20,handletextpad=0.4,handlelength=1.4,frameon=True,ncol=len(self.ref_temp)+1,columnspacing=1)
            _arr[0][0].add_artist(legend1)

            legend2_handles = []
            legend2_handles.append(Line2D([0],[0],marker='d',markersize=8,markerfacecolor='None',markeredgecolor='k', linestyle='None',label=r'P$_{\text{C1}}$ plane'))
            legend2_handles.append(Line2D([0],[0],marker='o',markersize=8,markerfacecolor='None',markeredgecolor='k', linestyle='None',label=r'P$_{\text{C2}}$ plane'))
            legend2 = _arr[0][0].legend(handles=legend2_handles, loc=1,bbox_to_anchor=(1.0,1.0),fontsize=16, handletextpad=0.4,handlelength=1.4,frameon=True,ncol=1,labelspacing=0.1,borderpad=0.2)
            _arr[0][0].add_artist(legend2)

        else:
            print('Complete legend not implemented yet')


#        plt.show()
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig.align_ylabels()
            fig2.align_ylabels()
            self.save_figure_split(fig,fig2)
        else:
            fig.align_ylabels()
            self.save_figure(fig)


    def plot_total_pgap_separate_indirect(self):

        if not self.zpr_fnames and self.te_fnames:
            raise Exception('For total pgap, both te_fnames and zpr_fnames must be provided')
        if not self.valence:
            raise Exception('For total pgap, must provide valence band index (index starting at 1)')

        file_qty = len(self.zpr_fnames)

        if len(self.te_fnames) != file_qty:
            raise Exception('zpr_fnames and te_fnames must contain the same number of files, but zpr_fnames has {} and te_fnames has {}.'.format(file_qty, len(te_fnames)))

        # Define figure
        fig, _arr = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            raise Exception('split not yet implemented for total_pgap')
            fig2, _arr2 = plt.subplots(3,1, figsize=self.figsize, squeeze=False, sharex=True)
            plt.subplots_adjust(hspace=0.05, top=0.95) 

        if not self.band_numbers:
            raise Exception('Must provide valence and conduction band numbers vis band_numbers')


        # Read and treat all input files
        for ifile in range(file_qty):
                
            # Define file class
            zpr_file = ZPRfile(self.zpr_fnames[ifile], read=False)
            zpr_file.read_nc()
            te_file = TEfile(self.te_fnames[ifile], read=False)
            te_file.read_nc()


            # Set parameters for this file
            self.nsppol = zpr_file.nsppol
            self.nkpt = zpr_file.nkpt
            self.max_band = te_file.max_band
            self.ntemp = zpr_file.ntemp
            self.temp = zpr_file.temp
            self.kpoints = zpr_file.kpoints
            if self.units == 'eV':
                self.eig0 = zpr_file.eig0*cst.ha_to_ev
                zpr_file.correction = zpr_file.correction*cst.ha_to_ev
            elif self.units == 'meV':
                self.eig0 = zpr_file.eig0*cst.ha_to_ev*1000
                zpr_file.correction = zpr_file.correction*cst.ha_to_ev*1E3

            # Checkups, put in a separate function...
            if te_file.nsppol != self.nsppol:
                raise Exception('zpr and te files (index {}) do not have the same nsppol.'.format(ifile))
            if te_file.nkpt != self.nkpt:
                raise Exception('zpr and te files (index {}) do not have the same nkpt.'.format(ifile))
            if zpr_file.max_band != self.max_band:
                zpr_max_band = zpr_file.max_band
            else:
                zpr_max_band = None
                raise Exception('zpr and te files (index {}) do not have the same max_band.'.format(ifile))
            if te_file.ntemp != self.ntemp:
                raise Exception('zpr and te files (index {}) do not have the same ntemp.'.format(ifile))
            if not np.array_equal(te_file.temp, self.temp):
                raise Exception('zpr and te files (index {}) do not have the same temperature array.'.format(ifile))
            if not np.array_equal(te_file.kpoints,self.kpoints):
                raise Exception('zpr and te files (index {}) do not have the same kpoint array.'.format(ifile))
            if not np.allclose(te_file.eig0, self.eig0[0,:,:self.max_band]):
                raise Exception('zpr and te files (index {}) do not have the same bare eigenvalues.'.format(ifile))


           ############################################### 
#            self.gridsize = self.zpr.gridsize
            self.gap_units = zpr_file.gap_energy_units

            input_units = te_file.gap_energy_units

            if input_units == self.units:
                pass
            else:
                if self.units == 'eV':
                    if input_units == 'Hartree':
                        te_file.eig0 = te_file.eig0*cst.ha_to_ev
                        te_file.correction = te_file.correction*cst.ha_to_ev
                    elif input_units == 'meV':
                         te_file.eig0 = te_file.eig0*1E-3
                         te_file.correction = te_file.correction*1E-3
                    else:
                        raise Exception('What kind of energy units are you using?!?')
                if self.units == 'meV':
                    if input_units == 'Hartree':
                        te_file.eig0 = te_file.eig0*cst.ha_to_ev*1000
                        te_file.correction = te_file.correction*cst.ha_to_ev*1000
                    elif input_units == 'eV':
                        te_file.eig0 = te_file.eig0*1E3
                        te_file.correction = te_file.correction*1E3
                    else:
                        raise Exception('What kind of energy units are you using?!?')

            # Set side of phase transition for split2 option
            if self.split2:
                if self.pressure[ifile] < self.crit_pressure:
                    side=0
                else:
                    side=1

            # Retrieve data from files
            self.correction = np.zeros((self.nsppol,self.nkpt, self.max_band, self.ntemp))
            self.eigcorr = np.zeros((self.nsppol,self.nkpt, self.max_band, self.ntemp))

            for t in range(self.ntemp):
                if self.max_band is not None:
                    self.correction[:,:,:,t] = te_file.correction[:,:,:,t] + zpr_file.correction[:,:,:self.max_band,t]
                    self.eigcorr[:,:,:,t] = self.eig0[:,:,:self.max_band] + self.correction[:,:,:,t]
                else:
                    self.correction = te_file.correction + zpr_file.correction
                    self.eigcorr = self.eig0 + self.correction
            
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
                    self.full_gap_ren_band = np.zeros((file_qty,self.ntemp,2))
                    self.full_gap_loc = np.zeros((file_qty, self.ntemp+1,2))
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp+1))
                    self.full_gap_energy_band = np.zeros((file_qty, self.ntemp+1,2))
                    self.ref_temp = self.temp
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

#            if self.split:
#                self.full_gap_ren[ifile,:,:] = self.gap_ren
#     #           print(np.shape(self.gap_ren))
#     #           print(np.shape(self.full_gap_ren))
#            if self.split2:
#                self.full_gap_ren[ifile,:] = self.gap_ren[:,side]
#            else:
#
#                self.full_gap_ren[ifile,:] = self.gap_ren

            if self.split:
                raise Exception('Split not ready yet! How come did you get this far?!?')

            if self.indirect:
                self.full_gap_loc[ifile,:,:], self.full_gap_energy[ifile,:], self.full_gap_ren[ifile,:], self.full_gap_energy_band[ifile,:,:],self.full_gap_ren_band[ifile,:,:] = self.get_indirect_gap_info(self.eig0, self.eigcorr,False)
            else:
                raise Exception('Direct gap not yet implemented, get to work!')

        print(self.full_gap_ren_band[:,-1,0])
        print(self.full_gap_ren_band[:,-1,1])
        print(self.full_gap_ren[:,-1])

        if self.crit_pressure is not None:
            crit_index = self.find_temp_index()
        else:
            crit_index = None

        for T in range(self.ntemp):

            if crit_index is not None:
                s = crit_index+1

                if self.split:
                    #figure1
                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren_band[0:s,T,1,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure[0:s], self.full_gap_ren_band[0:s,T,0,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,0], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    _arr[2][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.27,1.62), ncol = 1, fontsize=20)

                  #  self.set_legend_pgap2(_arr[2][0])
                    _arr[0][0].plot(self.pressure[s:], self.full_gap_ren_band[s:,T,1,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_gap_ren_band[s:,T,0,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,0], marker='o', markersize=8, linewidth=1.5, color=self.color[T])

                    #figure2
                    _arr2[0][0].plot(self.pressure[0:s], self.full_gap_ren_band[0:s,T,1,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[1][0].plot(self.pressure[0:s], self.full_gap_ren_band[0:s,T,0,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr2[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T,1], marker='d', markersize=8, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                    self.set_legend_pgap2(_arr2[2][0])
                    _arr2[0][0].plot(self.pressure[s:], self.full_gap_ren_band[s:,T,1,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr2[1][0].plot(self.pressure[s:], self.full_gap_ren_band[s:,T,0,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])
                    _arr2[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T,1], marker='o', markersize=8, linewidth=1.5, color=self.color[T])

                else:
                    _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren_band[0:s,T,1], marker='d', markersize=7, linewidth=1.5, color=self.color[T], label=r'{} K'.format(self.ref_temp[T]))
                    _arr[1][0].plot(self.pressure[0:s], self.full_gap_ren_band[0:s,T,0], marker='d', markersize=7, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[2][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='d', markersize=7, linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')

                  #  _arr[2][0].legend(numpoints=1, loc= 'center right', bbox_to_anchor=(1.27,1.62), ncol=1, fontsize=20)

                    _arr[0][0].plot(self.pressure[s:], self.full_gap_ren_band[s:,T,1], marker='o', markersize=7, linewidth=1.5, color=self.color[T])
                    _arr[1][0].plot(self.pressure[s:], self.full_gap_ren_band[s:,T,0], marker='o', markersize=7, linewidth=1.5, color=self.color[T])
                    _arr[2][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o', markersize=7, linewidth=1.5, color=self.color[T])


            else:
                if self.split:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)
                    _arr[1][0].legend(numpoints = 1, loc = 'center right', bbox_to_anchor=(1.02,1.0), ncol = 1, fontsize=26)


#                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
#                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)

                else:
                    _arr[0][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].plot(self.pressure, self.full_gap_ren[:,T], marker='d', linewidth=1.5, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    _arr[1][0].legend(numpoints = 1, loc = 'lower center', bbox_to_anchor=(0.5,-0.25), ncol = self.ntemp, fontsize=16)


        limms = [[-50.,50.],[-2.,50.],[-90.,40.]]
##########FIX ME : add default data and if condition
        if self.cond_ylims is not None:
            limms[0] = self.cond_ylims
        if self.val_ylims is not None:
            limms[1] = self.val_ylims
        if self.gap_ylims is not None:
            limms[2] = self.gap_ylims

        for i in range(3):
#            lims = _arr[i,0].get_ylim() 
            ylims = limms[i]    
            _arr[i,0].set_ylim(limms[i])
            self.set_vrefs(_arr[i][0], self.pressure, 0.,style='dashed')
            if self.crit_pressure:
                pc1 = self.crit_pressure
                self.set_hrefs(ylims, _arr[i][0], self.crit_pressure,'gray')
                if self.crit_pressure2:
                    pc2 = self.crit_pressure2
                    self.set_hrefs(ylims, _arr[i][0], self.crit_pressure+0.2,'gray')
                    _arr[i][0].fill([pc1,pc2,pc2,pc1],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2)
                else:
                    print('No value provided for PC2')
            else:
                print('No value provided for PC1')


            if self.split:
                self.set_vrefs(_arr2[i][0], self.pressure, 0.)
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure,'gray')
                self.set_hrefs(self.ylims, _arr2[i][0], self.crit_pressure+0.2,'gray')
                lims = _arr[i,0].get_ylim() 
                _arr2[i][0].fill([2.08,2.28,2.28,2.08],[lims[0]-5,lims[0]-5,lims[1]+5,lims[1]+5],'gray',alpha=0.2)


        
        self.set_xaxis(_arr[2][0], self.pressure)
        self.set_yaxis_separate(_arr[0][0], 'CB ren ({})'.format(self.gap_units),self.cond_ylims, self.cond_yticks)
        self.set_yaxis_separate(_arr[1][0], 'VB ren ({})'.format(self.gap_units),self.val_ylims, self.val_yticks)
        self.set_yaxis_separate(_arr[2][0], 'Gap ren ({})'.format(self.gap_units),self.gap_ylims, self.gap_yticks)


        self.set_yaxis(_arr[0][0], 'CB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[1][0], 'VB ren ({})'.format(self.gap_units))
        self.set_yaxis(_arr[2][0], 'Gap ren ({})'.format(self.gap_units))

#        if self.main_title:
#            self.set_title(_arr[0][0], self.main_title)

        fig.subplots_adjust(left=0.11,bottom=0.08,right=0.81,top=0.95,wspace=0.2,hspace=0.12)

        #_arr[1,0].text(self.explicit_pressures[8],300, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=30)
        #_arr[1,0].text(self.explicit_pressures[35],300, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=30)

        #_arr[1,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=28, weight='bold', color='gray')



#        if self.split:
#            self.set_xaxis(_arr2[2][0], self.pressure)
#            self.set_yaxis_separate(_arr2[0][0], 'CB ren ({})'.format(self.gap_units),self.cond_ylims, self.cond_yticks)
#            self.set_yaxis_separate(_arr2[1][0], 'VB ren ({})'.format(self.gap_units),self.val_ylims, self.val_yticks)
#            self.set_yaxis_separate(_arr2[2][0], 'Gap ren ({})'.format(self.gap_units),self.gap_ylims, self.gap_yticks)
#
#            self.set_title(_arr[0][0], self.title[0])
#            self.set_title(_arr2[0][0], self.title[1])
#
#        else:
#            if self.main_title:
#                self.set_title(_arr0[0][0], self.main_title)
#
#        # Custum stuff
#        _arr[2][0].text(0.7,-82, r'$\mathbb{Z}_2\!=\!0$', fontsize=24,color='k')
#        _arr[2][0].text(2.7, -82, r'$\mathbb{Z}_2\!=\!1$',fontsize = 24,color='k')
#        _arr[2][0].text(1.85,-60,r'$\Rightarrow$',fontsize=20,color='#5A5A5A')
#        _arr[2][0].text(1.35,-60,r'WSM',fontsize=20,color='#5A5A5A',weight='bold')
#
        legend_handles = [] 
        for t, temp in enumerate(self.ref_temp):
            legend_handles.append(Line2D([0],[0],color=self.color[t],linewidth=1.5, label=r'{:>3.0f} K'.format(self.ref_temp[t])))
#                Line2D([0],[0],color='b',marker='o',markersize=8,linestyle='None',label=r'P$_{\text{C2}}$ plane')]
        #legend_handles.append('')
        #legend_handles.append('')

        legend1 = _arr[0][0].legend(handles=legend_handles, loc=9,bbox_to_anchor=(0.5,1.3),fontsize=20, handletextpad=0.4,handlelength=1.4,frameon=True,ncol = len(self.ref_temp),columnspacing=1)
        _arr[0][0].add_artist(legend1)
        
        legend2_handles=[]
        legend2_handles.append(Line2D([0],[0],marker='d',markersize=8,markerfacecolor='None', markeredgecolor='k', linestyle='None',label=r'P$_{\text{C1}}$ plane'))
        legend2_handles.append(Line2D([0],[0],marker='o',markersize=8,markerfacecolor='None', markeredgecolor='k', linestyle='None',label=r'P$_{\text{C2}}$ plane'))
        legend2 = _arr[1][0].legend(handles=legend2_handles, loc=1,bbox_to_anchor=(1.0,0.28),fontsize=16, handletextpad=0.4,handlelength=1.4,frameon=True,ncol = 1,labelspacing=0.1,borderpad=0.2)
        _arr[1][0].add_artist(legend2)

        fig.subplots_adjust(hspace=0.0,top=0.90,right=0.95)

        if self.split:
            fig.align_ylabels()
            fig2.align_ylabels()
            self.save_figure_split(fig,fig2)
        else:
            fig.align_ylabels()
            self.save_figure(fig)

    def plot_total_pgap_indirect(self):

        # Define figure
        only=True  # only the Egap(P) for each T. If False, add gap renorm on top subplot
        extrapolate_bands = True

        if only:
            if self.diagram_together:
                fig, _arr = plt.subplots(1,2, figsize=(self.figsize[0]*2, self.figsize[1]), squeeze=False, sharex=True)
            else:
                fig, _arr = plt.subplots(1,1, figsize=self.figsize, squeeze=False, sharex=True)
        else:
            fig, _arr = plt.subplots(2,1, figsize=self.figsize, squeeze=False, sharex=True)
        #plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            raise Exception('split not yet implemented for te_pgap')
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

            self.explicit_pressures = self.exgap.explicit_pressures


        if not self.zpr_fnames and self.te_fnames:
            raise Exception('For total pgap, both te_fnames and zpr_fnames must be provided')
        if not self.valence:
            raise Exception('For total pgap, must provide valence band index (index starting at 1)')

        file_qty = len(self.zpr_fnames)

        if len(self.te_fnames) != file_qty:
            raise Exception('zpr_fnames and te_fnames must contain the same number of files, but zpr_fnames has {} and te_fnames has {}.'.format(file_qty, len(te_fnames)))

        if not self.band_numbers:
            raise Exception('Must provide valence and conduction band numbers vis band_numbers')

        if self.bare_fnames is not None:
            if len(self.bare_fnames) != file_qty:
                raise Exception('bare_fnames and zpr_fnames must contain the same number of files, but bare_fnames has {} and zpr_fnames has {}.'.format(len(self.bare_fnames), file_qty))


        figt,arrt = plt.subplots(5,6,figsize=(20,15))
        #figB, arrB = plt.subplots(1,1,squeeze=True,figsize=self.figsize)
        #colB = col = ['slategray','purple','blue','green','yellow','red']
        # Read and treat all input files
        for ifile in range(file_qty):
                
            #print('file #',ifile+1)
            # Define file class
            zpr_file = ZPRfile(self.zpr_fnames[ifile], read=False)
            zpr_file.read_nc()
            te_file = TEfile(self.te_fnames[ifile], read=False)
            te_file.read_nc()
            if self.bare_fnames is not None:
                bare_file = EIGfile(self.bare_fnames[ifile], read=False)
                bare_file.read_nc()


            # Set parameters for this file
            self.nsppol = zpr_file.nsppol
            self.max_band = te_file.max_band
            self.ntemp = zpr_file.ntemp
            self.temp = zpr_file.temp
            if self.read_correction_from_spline:
                self.kpoints = zpr_file.kpoints_spline
                zpr_correction = zpr_file.td_ren_spline
                self.nkpt = zpr_file.nkpt_spline
            else:
                self.kpoints = zpr_file.kpoints
                zpr_correction = zpr_file.correction
                self.nkpt = zpr_file.nkpt

            if self.units == 'eV':
                zpr_eig0 = zpr_file.eig0*cst.ha_to_ev
                zpr_correction = zpr_correction*cst.ha_to_ev
            elif self.units == 'meV':
                zpr_eig0 = te_file.eig0*cst.ha_to_ev*1000
                zpr_correction = zpr_correction*cst.ha_to_ev*1E3

            #print('EPI: ', zpr_file.indirect_gap_ren)
            #if ifile==4:
            #    print('Gap location, EPI:',zpr_file.indirect_gap_location[0,:,:])
            #    print('Gap location, TE:',te_file.indirect_gap_location[0,:,:])

            # Checkups, put in a separate function...
            if te_file.nsppol != self.nsppol:
                raise Exception('zpr and te files (index {}) do not have the same nsppol.'.format(ifile))
            if te_file.nkpt != self.nkpt:
                raise Exception('zpr and te files (index {}) do not have the same nkpt. Did you want to read EPI correction from spline interpolation?'.format(ifile))
            if zpr_file.max_band != self.max_band:
                zpr_max_band = zpr_file.max_band
            else:
                zpr_max_band = None
                raise Exception('zpr and te files (index {}) do not have the same max_band.'.format(ifile))
            if te_file.ntemp != self.ntemp:
                raise Exception('zpr and te files (index {}) do not have the same ntemp.'.format(ifile))
            if not np.array_equal(te_file.temp, self.temp):
                raise Exception('zpr and te files (index {}) do not have the same temperature array.'.format(ifile))
            if not np.allclose(te_file.kpoints,self.kpoints):
                raise Exception('zpr and te files (index {}) do not have the same kpoint array.\n zpr has {} kpts and te has {}.'.format(ifile, self.kpoints.shape[0],te_file.kpoints.shape[0]))
            if not self.read_correction_from_spline:
                if not np.allclose(te_file.eig0, zpr_eig0[0,:,:self.max_band]):
                    raise Exception('zpr and te files (index {}) do not have the same bare eigenvalues.'.format(ifile))


           ############################################### 
#            self.gridsize = self.zpr.gridsize
            self.gap_units = zpr_file.gap_energy_units

            input_units = te_file.gap_energy_units

            if input_units == self.units:
                pass
            else:
                if self.units == 'eV':
                    if input_units == 'Hartree':
                        te_file.eig0 = te_file.eig0*cst.ha_to_ev
                        te_file.correction = te_file.correction*cst.ha_to_ev
                    elif input_units == 'meV':
                         te_file.eig0 = te_file.eig0*1E-3
                         te_file.correction = te_file.correction*1E-3
                    else:
                        raise Exception('What kind of energy units are you using?!?')
                if self.units == 'meV':
                    if input_units == 'Hartree':
                        te_file.eig0 = te_file.eig0*cst.ha_to_ev*1000
                        te_file.correction = te_file.correction*cst.ha_to_ev*1000
                    elif input_units == 'eV':
                        te_file.eig0 = te_file.eig0*1E3
                        te_file.correction = te_file.correction*1E3
                    else:
                        raise Exception('What kind of energy units are you using?!?')

            #print('TE :', te_file.indirect_gap_ren)#, te_file.indirect_gap_ren[5])
#            print('sum:', te_file.indirect_gap_ren[0]+zpr_file.indirect_gap_ren[0])
            if self.bare_fnames is not None:
                #print('Defining bare eigenvalues from bare_fnames')
                self.eig0 = bare_file.eig0
                if self.units == 'eV':
                    self.eig0 = self.eig0*cst.ha_to_ev
                elif self.units == 'meV':
                    self.eig0 = self.eig0*cst.ha_to_ev*1E3
            elif self.read_correction_from_spline:
                self.eig0 = te_file.eig0
            else:
                self.eig0 = zpr_file.eig0

            # Set side of phase transition for split2 option
            if self.split2:
                if self.pressure[ifile] < self.crit_pressure:
                    side=0
                else:
                    side=1

            # Retrieve data from files
            self.correction = np.zeros((self.nsppol,self.nkpt, self.max_band, self.ntemp))
            self.eigcorr = np.zeros((self.nsppol,self.nkpt, self.max_band, self.ntemp))

#            print(te_file.correction[0,:,27,-1])
            for t in range(self.ntemp):
                if self.max_band is not None:
                    self.correction[:,:,:,t] = te_file.correction[:,:,:,t] + zpr_correction[:,:,:self.max_band,t]
                    self.eigcorr[:,:,:,t] = self.eig0[:,:,:self.max_band] + self.correction[:,:,:,t]
                else:
                    self.correction = te_file.correction + zpr_correction
                    self.eigcorr = self.eig0 + self.correction
#            print(self.correction[0,:,27,-1])
            ###here###
            
        #    print(te_file.correction[0,440:480,28,0])
            
            # Initialize plotting arrays
            if ifile==0:
                if self.split: # distinct data for left and rifht gap
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2,2)) # file, temp, vbcb, lr 
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((file_qty))

                if self.split2: # single array, left gap for trovoal phase and right gap for topol phase
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.full_band_ren = np.zeros((file_qty,self.ntemp,2))
                    self.ref_temp = self.temp
                    self.full_energy0 = np.zeros((file_qty))

                else:
                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
                    self.full_gap_ren_band = np.zeros((file_qty,self.ntemp,2))
                    self.full_gap_loc = np.zeros((file_qty, self.ntemp+1,2))
                    self.full_gap_energy = np.zeros((file_qty, self.ntemp+1))
                    self.full_gap_energy_band = np.zeros((file_qty, self.ntemp+1,2))
                    self.full_energy0 = np.zeros((file_qty))
                    self.gap_loc_coordinates = np.zeros((file_qty, self.ntemp+1,2,3))
                    self.rashba_energy = np.zeros((file_qty, self.ntemp+1, 2)) # P,T,VB/CB
                    full_energy = np.zeros((file_qty,self.ntemp))
                    self.ref_temp = self.temp
            else:
                if np.array_equal(self.temp,self.ref_temp) == False:
                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))

            if self.ntemp > len(self.color):
                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

#            if self.split:
#                self.full_gap_ren[ifile,:,:] = self.gap_ren
#     #           print(np.shape(self.gap_ren))
#     #           print(np.shape(self.full_gap_ren))
#            if self.split2:
#                self.full_gap_ren[ifile,:] = self.gap_ren[:,side]
#            else:
#
#                self.full_gap_ren[ifile,:] = self.gap_ren

            if self.split:
                raise Exception('Split not ready yet! How come did you get this far?!?')

            if self.indirect:
                self.full_gap_loc[ifile,:,:], self.full_gap_energy[ifile,:], self.full_gap_ren[ifile,:], self.full_gap_energy_band[ifile,:,:],self.full_gap_ren_band[ifile,:,:] = self.get_indirect_gap_info(self.eig0, self.eigcorr,False)
                for tt in range(self.ntemp+1):
                    self.gap_loc_coordinates[ifile,tt,0,:] = self.kpoints[int(self.full_gap_loc[ifile,tt,0])]
                    self.gap_loc_coordinates[ifile,tt,1,:] = self.kpoints[int(self.full_gap_loc[ifile,tt,1])]

                self.rashba_energy[ifile,:,:] = self.get_rashba_info(ifile)

            else:
                raise Exception('Direct gap not yet implemented, get to work!')

            if not self.gap_fname:
                self.full_energy0[ifile] = self.full_gap_energy[ifile,0]
            else:
                # find P in explicit pressures array
                p_index = self.find_pressure_index(self.pressure[ifile])
                self.full_energy0[ifile] = self.explicit_gap_energies[p_index]

            # Plotting for test purposes
            if ifile == 8:
                tarr = [0,1,2,3,4,5]
                #typ='eig'
                kptarr = range(len(self.kpoints))
#                figt, arrt = plt.subplots(2,5)
#                arrt.plot(kptarr,self.eig0[0,:,27:29],'k')
                
                for ti in tarr:
                    #if typ == 'corr':
                    arrt[3,ti].plot(kptarr,zpr_correction[0,:,27,ti],'b',label='EPI')
                    arrt[3,ti].plot(kptarr,te_file.correction[0,:,27,ti], 'r', label='TE')
                    arrt[3,ti].plot(kptarr,self.correction[0,:,27,ti], 'g', label='EPI+TE')

                    arrt[1,ti].plot(kptarr,zpr_correction[0,:,28,ti],'b',label='EPI')
                    arrt[1,ti].plot(kptarr,te_file.correction[0,:,28,ti], 'r', label='TE')
                    arrt[1,ti].plot(kptarr,self.correction[0,:,28,ti], 'g', label='EPI+TE')
                    
                    # total gap correction
                    arrt[4,ti].plot(kptarr,zpr_correction[0,:,28,ti]-zpr_correction[0,:,27,ti],'b',label='EPI')
                    arrt[4,ti].plot(kptarr,te_file.correction[0,:,28,ti]-te_file.correction[0,:,27,ti], 'r', label='TE')
                    arrt[4,ti].plot(kptarr,self.correction[0,:,28,ti]-self.correction[0,:,27,ti], 'g', label='EPI+TE')

                    #elif typ == 'eig':
                    arrt[2,ti].plot(kptarr, self.eig0[0,:,27],'k', label='bare')
                    arrt[2,ti].plot(kptarr,self.eig0[0,:,27]+zpr_correction[0,:,27,ti],'b',label='EPI')
                    arrt[2,ti].plot(kptarr,self.eig0[0,:,27]+te_file.correction[0,:,27,ti], 'r', label='TE')
                    arrt[2,ti].plot(kptarr,self.eigcorr[0,:,27,ti], 'g', label='EPI+TE')

                    arrt[0,ti].plot(kptarr, self.eig0[0,:,28],'k', label='bare')
                    arrt[0,ti].plot(kptarr,self.eig0[0,:,28]+zpr_correction[0,:,28,ti],'b',label='EPI')
                    arrt[0,ti].plot(kptarr,self.eig0[0,:,28]+te_file.correction[0,:,28,ti], 'r', label='TE')
                    arrt[0,ti].plot(kptarr,self.eigcorr[0,:,28,ti], 'g', label='EPI+TE')


                    
                    v1 = self.find_gap_index(zpr_file.indirect_gap_location[ti,0,:])
                    c1 = self.find_gap_index(zpr_file.indirect_gap_location[ti,1,:])
                    v2 = self.find_gap_index(te_file.indirect_gap_location[ti,0,:])
                    c2 = self.find_gap_index(te_file.indirect_gap_location[ti,1,:])
                    v3 = int(self.full_gap_loc[ifile,ti+1,0] )
                    c3 = int(self.full_gap_loc[ifile,ti+1,1])

                    bare_vbm = np.argmax(self.eig0[0,:,27])
                    bare_cbm = np.argmin(self.eig0[0,:,28])

                    #if typ == 'corr':
                    arrt[3,ti].plot(v1,zpr_correction[0,v1,27,ti], 'bx')
                    arrt[3,ti].plot(v2,te_file.correction[0,v2,27,ti], 'rx')
                    arrt[3,ti].plot(v3,self.correction[0,v3,27,ti], 'gx')
                    arrt[1,ti].plot(c1,zpr_correction[0,c1,28,ti], 'bx')
                    arrt[1,ti].plot(c2,te_file.correction[0,c2,28,ti], 'rx')
                    arrt[1,ti].plot(c3,self.correction[0,c3,28,ti], 'gx')
                    #elif typ == 'eig': 
                    arrt[2,ti].plot(v1,self.eig0[0,v1,27]+zpr_correction[0,v1,27,ti], 'bx')
                    arrt[2,ti].plot(v2,self.eig0[0,v2,27]+te_file.correction[0,v2,27,ti], 'rx')
                    arrt[2,ti].plot(v3,self.eigcorr[0,v3,27,ti], 'gx')
                    arrt[0,ti].plot(c1,self.eig0[0,c1,28]+zpr_correction[0,c1,28,ti], 'bx')
                    arrt[0,ti].plot(c2,self.eig0[0,c2,28]+te_file.correction[0,c2,28,ti], 'rx')
                    arrt[0,ti].plot(c3,self.eigcorr[0,c3,28,ti], 'gx')
                    #arrt[1,ti].plot(bare_vbm*np.ones((10)), np.linspace(arrt[1,ti].get_ylim()[0], arrt[1,ti].get_ylim()[1], 10), 'k:')
                    #arrt[0,ti].plot(bare_cbm*np.ones((10)), np.linspace(arrt[0,ti].get_ylim()[0], arrt[0,ti].get_ylim()[1], 10), 'k:')



#                    arrt[0,ti].set_xlim(100,200)
#                    arrt[1,ti].set_xlim(100,200)
                    arrt[0,ti].set_xlim(0,260)
                    arrt[2,ti].set_xlim(0,260)
                    arrt[1,ti].set_xlim(0,260)
                    arrt[3,ti].set_xlim(0,260)
                    arrt[4,ti].set_xlim(0,260)


                    arrt[0,ti].set_ylim(1100,2000)
                    arrt[2,ti].set_ylim(1000,1500)
                    arrt[2,ti].plot(bare_vbm*np.ones((10)), np.linspace(arrt[2,ti].get_ylim()[0], arrt[2,ti].get_ylim()[1], 10), 'k:')
                    arrt[0,ti].plot(bare_cbm*np.ones((10)), np.linspace(arrt[0,ti].get_ylim()[0], arrt[0,ti].get_ylim()[1], 10), 'k:')
                    arrt[3,ti].plot(bare_vbm*np.ones((10)), np.linspace(arrt[3,ti].get_ylim()[0], arrt[3,ti].get_ylim()[1], 10), 'k:')
                    arrt[1,ti].plot(bare_cbm*np.ones((10)), np.linspace(arrt[1,ti].get_ylim()[0], arrt[1,ti].get_ylim()[1], 10), 'k:')
                    arrt[2,ti].plot(v3*np.ones((10)), np.linspace(arrt[2,ti].get_ylim()[0], arrt[2,ti].get_ylim()[1], 10), 'g:')
                    arrt[0,ti].plot(c3*np.ones((10)), np.linspace(arrt[0,ti].get_ylim()[0], arrt[0,ti].get_ylim()[1], 10), 'g:')
                    arrt[3,ti].plot(v3*np.ones((10)), np.linspace(arrt[3,ti].get_ylim()[0], arrt[3,ti].get_ylim()[1], 10), 'g:')
                    arrt[1,ti].plot(c3*np.ones((10)), np.linspace(arrt[1,ti].get_ylim()[0], arrt[1,ti].get_ylim()[1], 10), 'g:')
                    arrt[4,ti].plot(c3*np.ones((10)), np.linspace(arrt[4,ti].get_ylim()[0], arrt[4,ti].get_ylim()[1], 10), ':', color='purple')
                    arrt[4,ti].plot(v3*np.ones((10)), np.linspace(arrt[4,ti].get_ylim()[0], arrt[4,ti].get_ylim()[1], 10), ':', color='pink')




                    for a1,a2 in itt.product(range(5),range(6)):
                        arrt[a1,a2].plot(np.linspace(0,260,10), np.zeros((10)), 'k')

    #                arrt.set_ylim(1400,1600)
                    arrt[0,ti].set_title('{} K'.format(self.temp[ti]))
                    
                figt.suptitle('{}GPa'.format(self.pressure[ifile]))
                arrt[0,0].set_ylabel('CBM eigcorr')
                arrt[2,0].set_ylabel('VBM eigcorr')
                arrt[1,0].set_ylabel('CBM corr')
                arrt[3,0].set_ylabel('VBM corr')
                arrt[4,0].set_ylabel('Gap corr')
                #if typ == 'corr':
                #    plt.savefig('totalcorr_0p5gpa_{}K.png'.format(self.temp[ti]))
                #    os.system('open totalcorr_0p5gpa_{}K.png'.format(self.temp[ti]))
                #elif typ == 'eig':
                plt.savefig('shift2d_1gpa.png'.format(self.temp[ti]))
                os.system('open shift2d_1gpa.png'.format(self.temp[ti]))

####### unchanged stuff strts here

#                  # Initialize full arrays
#            if ifile==0:
#                if self.split:
#                    self.full_gap_energy = np.zeros((file_qty, self.ntemp,2))
#                    self.full_gap_ren = np.zeros((file_qty, self.ntemp,2))
#                    self.ref_temp = self.temp
#                    self.full_energy0 = np.zeros((len(self.explicit_pressures),2))
#
#                if self.split2: # single array, left gap for trovoal phase and right gap for topol phase
#                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
#                    self.full_gap_energy = np.zeros((file_qty, self.ntemp))
#                    self.ref_temp = self.temp
#                    self.full_energy0 = np.zeros((file_qty))
#
#                else:
#                    self.full_gap_energy = np.zeros((file_qty, self.ntemp))
#                    self.full_gap_ren = np.zeros((file_qty, self.ntemp))
#                    self.ref_temp = self.temp
#                    self.full_energy0 = np.zeros((file_qty))
#            else:
#                if np.array_equal(self.temp,self.ref_temp) == False:
#                    raise Exception('All files must have the same temperature array! Please correct file #{}.'.format(ifile+1))
#
#            if self.ntemp > len(self.color):
#                raise Exception('Color vector is not long enough! Please provide {} color list.'.format(self.ntemp))

#            print(self.gap_energy)
#            if self.split:
#                self.full_gap_energy[ifile,:,:] = self.gap_energy
#                self.full_gap_ren[ifile,:,:] = self.gap_ren
#                if not self.gap_fname:
#                    self.full_energy0[ifile,:] = self.energy0
#                else:
#                    # find P in explicit pressures array
#                    p_index = self.find_pressure_index(self.pressure[ifile])
#                    self.full_energy0[ifile,:] = self.explicit_gap_energies[p_index]
#
#            if self.split2:
#                self.full_gap_ren[ifile,:] = self.gap_ren[:,side]
#                self.full_gap_energy[ifile,:] = self.gap_energy[:,side]
#                print(self.gap_fname)
#                if not self.gap_fname:
#                    self.full_energy0[ifile] = self.energy0[ifile,side]
#                else:
#                    # find P in explicit pressures array
#                    p_index = self.find_pressure_index(self.pressure[ifile])
#                    self.full_energy0[ifile] = self.explicit_gap_energies[p_index]
#            else:
#                self.full_gap_energy[ifile,:] = self.gap_energy
#                self.full_gap_ren[ifile,:] = self.gap_ren
#                if not self.gap_fname:
#                    self.full_energy0[ifile] = self.energy0
#                else:
#                    p_index = self.find_pressure_index(self.pressure[ifile])
#                    print(p_index)
#                    self.full_energy0[ifile] = self.explicit_gap_energies[p_index]

            if self.gap_fname is not None:
                if self.split:
                    for a in range(2):
                #        self.full_energy0[ifile,a] = self.explicit_gap_energies[ifile]
                        self.full_gap_energy[ifile,:,a] =  self.full_gap_ren[ifile,:,a] + np.ones(self.ntemp)*self.full_energy0[ifile,a]
                elif self.split2:
                    full_energy[ifile,:] = self.full_gap_ren[ifile,:] + np.ones(self.ntemp)*self.full_energy0[ifile]
                else:
                 #   self.full_energy0[ifile] = self.explicit_gap_energies[ifile]
                    full_energy[ifile,:] =  self.full_gap_ren[ifile,:] + np.ones(self.ntemp)*self.full_energy0[ifile]


#            if ifile == 4:
#                kptarr = range(len(self.kpoints))
#                bbands = [27,28]
#                for j in bbands:
#                    arrB.plot(kptarr, self.eig0[0,:,j], 'k')
#                    arrB.plot(kptarr, self.eig0[0,:,j] + zpr_correction[0,:,j,0], 'bx')
#                    arrB.plot(kptarr, self.eig0[0,:,j] + te_file.correction[0,:,j,0],'rx')
#                    arrB.plot(kptarr, self.eigcorr[0,:,j,0], 'g')
#                    #arrB.set_ylim(1925,2080)
#                    #arrB.set_xlim(440,480)
#                plt.savefig('test_total.png')
#                os.system('open test_total.png')
#        print('conduction')
#        print(self.full_gap_ren_band[:,:,1])
#        print('valence')
#        print(self.full_gap_ren_band[:,:,0])
#        print('gap')
#        print(self.full_gap_ren[:,:])

#        print(self.kpoints.shape)
#        print(self.full_gap_loc[4,1,:])
#        print('Gap location EPI+TE:',self.kpoints[int(self.full_gap_loc[4,1,0]),:], self.kpoints[int(self.full_gap_loc[4,1,1]),:])

        self.full_gap_energy = full_energy
        if self.crit_pressure is not None:
            crit_index = self.find_crit_pressure(self.crit_pressure, self.pressure)
            if self.gap_fname is not None:
                crit_index2 = self.find_crit_pressure(self.crit_pressure2, self.explicit_pressures)
                print('crit_index2',crit_index2)
        else:
            crit_index = None


        if self.follow:
            self.full_energy0[-1,:] = -self.full_energy0[-1,:]
            self.full_gap_energy[-1,:,:] = -self.full_gap_energy[-1,:,:]

        #print('EPI+TE :')
        #print(self.full_gap_ren[:,0])

    # Extrapolate behavior linearly towards TPT... it IS a bit sketchy.
        pressure1 = np.arange(1.5,2.0,0.1)
        tmparr = np.arange(2.1,2.6,0.1) # change 2.2 to 2.6 foir full BZ
        pressure2 = np.arange(2.3,3.0,0.1)
        tmparr2 = np.arange(2.0,2.3,0.1)
#        self.extr_pressure2 = np.arange(2.3, 3.0, 0.1)

        # linearize unperturbed gap energy to extrapolate inside WSM phase
        y0,y1 = np.polyfit(self.explicit_pressures[15:20], self.explicit_gap_energies[15:20],1)
        tmp2 = y1 + tmparr*y0
        self.extr_pressure1 = np.concatenate((pressure1, tmparr))

        y0,y1 = np.polyfit(self.explicit_pressures[21:29], self.explicit_gap_energies[21:29],1)
        tmp22 = y1 + tmparr2*y0
        self.extr_pressure2 = np.concatenate((tmparr2,pressure2))

        self.extr_full_gap_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp))
        self.extr_full_gap_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp))
        self.extr_full_gap_energy1 = np.zeros((len(self.extr_pressure1), self.ntemp))
        self.extr_full_gap_energy2 = np.zeros((len(self.extr_pressure2), self.ntemp))
        self.extr_full_band_ren1 = np.zeros((len(self.extr_pressure1), self.ntemp,2)) #pressure, temperature, val/cond
        self.extr_full_band_ren2 = np.zeros((len(self.extr_pressure2), self.ntemp,2))    
        extr_energy01 = np.zeros((len(pressure1)))
        extr_energy02 = np.zeros((len(pressure2)))
#        self.extr_full_energy02 = np.zeros((len(self.extr_pressure2)))

        # Get full_gap_energy0 (explicit bare gap energies) from explicit data
        for p, pres in enumerate(pressure1):
            ind = self.find_pressure_index(pres)
            extr_energy01[p] = self.explicit_gap_energies[ind]

        self.extr_full_energy01 = np.concatenate((extr_energy01,tmp2))

        for p, pres in enumerate(pressure2):
            ind = self.find_pressure_index(pres)
            extr_energy02[p] = self.explicit_gap_energies[ind]

        self.extr_full_energy02 = np.concatenate((tmp22,extr_energy02))

        #Store extrapolated phase boundaries (T)
        self.pc1 = np.zeros((self.ntemp))
        self.pc2 = np.zeros((self.ntemp))

        for T in range(self.ntemp):

            if self.split:
                if extrapolate_bands:
                    raise Exception('Extrapolate_bands not implemented for split yet')
                else:
                    #trivial side
                    x0,x1 = np.polyfit(self.pressure[0:crit_index+1], self.full_gap_ren[0:crit_index+1,T,0], 1)
                    self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]

                    # topol side
                    x0,x1 =  np.polyfit(self.pressure[crit_index+1:], self.full_gap_ren[crit_index+1:,T,0], 1)
                    self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]

            elif self.split2:
                if extrapolate_bands:
                    raise Exception('Extrapolate_bands not implemented for split yet')
                else:
                    #trivial side
                    x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_gap_ren[crit_index-1:crit_index+1,T], 1)
                    self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]

                    # topol side
                    x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_gap_ren[crit_index+1:crit_index+3,T], 1)
                    self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]


            else:
                if extrapolate_bands:
                    # Loop on band extrema:
                    print('\n\nFor T={} K:'.format(self.temp[T]))

                    for b in range(2):
                        #trivial side
                        x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_gap_ren_band[crit_index-1:crit_index+1,T,b], 1)
                        self.extr_full_band_ren1[:,T,b] = x1 + x0*self.extr_pressure1

                        # topol side
                        x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_gap_ren_band[crit_index+1:crit_index+3,T,b], 1)
                        self.extr_full_band_ren2[:,T,b] = x1 + x0*self.extr_pressure2

                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_band_ren1[:,T,1] - self.extr_full_band_ren1[:,T,0]
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_band_ren2[:,T,1] - self.extr_full_band_ren2[:,T,0]

                    x0,x1 = np.polyfit(self.extr_pressure1,self.extr_full_gap_energy1[:,T],1)
                    print('trivial side : {} GPa'.format(-x1/x0))
                    self.pc1[T] = -x1/x0

                    x0,x1 = np.polyfit(self.extr_pressure2,self.extr_full_gap_energy2[:,T],1)
                    print('topol side : {} GPa'.format(-x1/x0))
                    self.pc2[T] = -x1/x0

                else:
                    #trivial side
                    x0,x1 = np.polyfit(self.pressure[crit_index-1:crit_index+1], self.full_gap_ren[crit_index-1:crit_index+1,T], 1)
                    print('For T={} K:'.format(self.temp[T]))
                    self.extr_full_gap_ren1[:,T] = x1 + x0*self.extr_pressure1
                    self.extr_full_gap_energy1[:,T] = self.extr_full_energy01 + self.extr_full_gap_ren1[:,T]
                    x0,x1 = np.polyfit(self.extr_pressure1,self.extr_full_gap_energy1[:,T],1)
                    print('trivial side : {} GPa'.format(-x1/x0))
                    self.pc1[T] = -x1/x0
                    # topol side
                    x0,x1 =  np.polyfit(self.pressure[crit_index+1:crit_index+3], self.full_gap_ren[crit_index+1:crit_index+3,T], 1)
                    self.extr_full_gap_ren2[:,T] = x1 + x0*self.extr_pressure2
                    self.extr_full_gap_energy2[:,T] = self.extr_full_energy02 + self.extr_full_gap_ren2[:,T]
                    x0,x1 = np.polyfit(self.extr_pressure2,self.extr_full_gap_energy2[:,T],1)
                    print('topol side : {} GPa'.format(-x1/x0))
                    self.pc2[T] = -x1/x0

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
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='Static T=0')
                            if self.crit_index2 is not None:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:crit_index2+2], self.explicit_gap_energies[:crit_index2+2],marker='d', color='black', label='Static T=0')

                            else:
                                if only:
                                    _arr[0][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', color='black', label='Static T=0')

                        else:
                            if only:
                                _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='Static T=0')
                            else:
                                _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black', label='Static T=0')
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
                                    _arr[0][0].plot(self.explicit_pressures[crit_index2+2,:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[crit_index2+2:], self.explicit_gap_energies[crit_index2+2:],marker='d', color='black', label='Static T=0')

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
                                _arr2[0][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                            else:
                                arr2[1][0].plot(self.explicit_pressures[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                        else:
                            if only:
                                _arr2[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
                            else:
                                _arr2[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,1],marker='d', color='black', label='Static T=0')
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

#########################################
       #-------- crit_index but no split or split2
                else:
                    if T==self.temp[0]:
                        ###### FIX ME : keep explicit_pressures and explicit_gap_energies, just get a new crit index for this array! 
                        if self.gap_fname is not None:
#                            _arr[1][0].plot(self.explicit_pressures[0:s], self.explicit_gap_energies[0:s],marker='d', color='black', label='Static T=0')
                            if crit_index2 is not None:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[:crit_index2+1], self.explicit_gap_energies[:crit_index2+1],marker='d',markersize=8, color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:crit_index2+1], self.explicit_gap_energies[:crit_index2+1],marker='d',markersize=8, color='black', label='Static T=0')
                            else:
                                if only:
                                     _arr[0][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d', markersize=8,color='black', label='Static T=0')
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[:], self.explicit_gap_energies[:],marker='d',markersize=8, color='black', label='Static T=0')

                        else:
                            if only:
                                _arr[0][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black', label='Static T=0')
                            else:
                                _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='d', color='black', label='Static T=0')
#                    else:
#                        _arr[1][0].plot(self.pressure[0:s], self.full_energy0[0:s,0],marker='o', color='black')
    
                    if only:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='d',markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
                    else:
                        _arr[0][0].plot(self.pressure[0:s], self.full_gap_ren[0:s,T], marker='d',markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
        #                self.set_legend_gap(_arr[0][0])
                        _arr[1][0].plot(self.pressure[0:s], self.full_gap_energy[0:s,T], marker='d',markersize=8, linewidth=2.0, color=self.color[T], label=str(self.ref_temp[T])+' K')
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
                                     _arr[0][0].plot(self.explicit_pressures[crit_index2+1:], self.explicit_gap_energies[crit_index2+1:],marker='o',markersize=8, color='black', label=None)
                                else:
                                    _arr[1][0].plot(self.explicit_pressures[crit_index2+1:], self.explicit_gap_energies[crit_index2+1:],marker='o',markersize=8, color='black', label=None)

                    if only:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o',markersize=8, linewidth=2.0, color=self.color[T])
                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])

                    else:
                        _arr[0][0].plot(self.pressure[s:], self.full_gap_ren[s:,T], marker='o',markersize=8, linewidth=2.0, color=self.color[T])
                        _arr[1][0].plot(self.pressure[s:], self.full_gap_energy[s:,T], marker='o',markersize=8, linewidth=2.0, color=self.color[T])

                        _arr[0,0].plot(self.extr_pressure2, self.extr_full_gap_ren2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])
                        _arr[1,0].plot(self.extr_pressure2, self.extr_full_gap_energy2[:,T], linestyle='--', linewidth=2.0, label = None, color = self.color[T])


                    if self.gap_fname is not None:
                        pass
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
                    raise Exception('Only keyword not implemented yet without critical pressure index')
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
        if self.gap_ylims is not None:
            limms[0] = self.gap_ylims
        if self.egap_ylims is not None:
            limms[1] = self.egap_ylims

        if only:
            ylims = limms[1]    
            self.set_vrefs(_arr[0][0], self.pressure, 0.)
#            if i==0:
            if self.crit_pressure:
                _arr[0][0].set_rasterization_zorder(0)
                pc1 = self.crit_pressure
                self.set_hrefs(ylims, _arr[0][0], pc1,'gray')

                if self.crit_pressure2:
                    pc2 = self.crit_pressure2
                    self.set_hrefs(ylims, _arr[0][0], pc2,'gray')
                    _arr[0][0].fill([pc1,pc2,pc2,pc1],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2,zorder=-1, rasterized=True)
                else:
                    print('No value was given for PC2')
            else:
                print('No value was given for PC1')
#            self.set_hrefs(ylims, _arr[0,0], 2.18,'orange')
#            self.set_hrefs(ylims, _arr[0,0], 2.50,'orange')
#            _arr[0,0].fill([2.18,2.50,2.50,2.18],[ylims[0],ylims[0],ylims[1],ylims[1]],'orange',alpha=0.2)
#            self.set_hrefs(ylims, _arr[0,0], 2.68,'magenta')


        else:
            for i in range(2):
    #            lims = _arr[i,0].get_ylim() 
                _arr[i][0].set_rasterization_order(0)
    #            ylims = lims
                ylims = limms[i]    
                self.set_vrefs(_arr[i][0], self.pressure, 0.)
    #            if i==0:
                pc1 = self.crit_pressure
                pc2 = self.crit_pressure2

                self.set_hrefs(ylims, _arr[i][0], pc1,'black')
                self.set_hrefs(ylims, _arr[i][0], pc2,'black')
                _arr[i][0].fill([pc1,pc2,pc2,pc1],[ylims[0],ylims[0],ylims[1],ylims[1]],'gray',alpha=0.2,zorder=-1, rasterized=True)
    #            else:
                self.set_hrefs(ylims, _arr[i,0], 2.28,'black')
                self.set_hrefs(ylims, _arr[i,0], 2.65,'black')
                _arr[i,0].fill([2.28,2.65,2.65,2.28],[ylims[0],ylims[0],ylims[1],ylims[1]],'orange',alpha=0.2,zorder=-1, rasterized=True)
 

        if only:
            self.set_xaxis(_arr[0][0], self.pressure)
            self.set_yaxis(_arr[0][0], 'Gap energy ({})'.format(self.gap_units))
            _arr[0,0].set_ylim(limms[1])

        else:
            self.set_xaxis(_arr[1][0], self.pressure)
            self.set_yaxis(_arr[0][0], 'Gap renormalization ({})'.format(self.gap_units))
            self.set_yaxis(_arr[1][0], 'Gap energy ({})'.format(self.gap_units))

            _arr[1,0].set_ylim(limms[1])


#        if only:
#            _arr[0,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=24) 
#            _arr[0,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=24) 
#
##            _arr[0,0].text(self.explicit_pressures[5],290, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=26)
##            _arr[0,0].text(self.explicit_pressures[38],290, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=26)
#
#
#            _arr[0,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=26, weight='bold', color='gray')
#            _arr[0,0].text(self.explicit_pressures[25],250, r'WSM 300K', fontsize=26, weight='bold', color='orange')
#
#        else:
#            _arr[1,0].text(self.explicit_pressures[12],180, r'$\mathbf{Z_2=0}$', fontsize=24) 
#            _arr[1,0].text(self.explicit_pressures[30],180, r'$\mathbf{Z_2=1}$', fontsize=24) 
#
##            _arr[1,0].text(self.explicit_pressures[5],290, r'$\frac{k_z c}{2\pi} = 0.5$',fontsize=26)
##            _arr[1,0].text(self.explicit_pressures[38],290, r'$\frac{k_z c}{2\pi} = 0.5349$',fontsize=26)
#
#
#            _arr[1,0].text(self.explicit_pressures[12],250, r'WSM 0K', fontsize=26, weight='bold', color='gray')
#            _arr[1,0].text(self.explicit_pressures[25],250, r'WSM 300K', fontsize=26, weight='bold', color='orange')



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

        # Custom text
        if only:
            _arr[0][0].text(0.7,250, r'$\mathbb{Z}_2\!=\!0$',fontsize=32,color='k')
            _arr[0][0].text(2.7,250, r'$\mathbb{Z}_2\!=\!1$',fontsize=32,color='k')
            _arr[0][0].text(1.90,200,r'$\Rightarrow$',fontsize=26,color='#5A5A5A')
            _arr[0][0].text(1.44,210,r'Static',fontsize=26,color='#5A5A5A',weight='bold')
            _arr[0][0].text(1.44,190,r'WSM',fontsize=26,color='#5A5A5A',weight='bold')

            legend_handles = []
            legend_handles.append(Line2D([0],[0],color='k',linewidth=2.5,label=r'Static'))
            for t,T in enumerate(self.ref_temp):
                legend_handles.append(Line2D([0],[0],color=self.color[t],linewidth=2.5,label=r'{:>3.0f} K'.format(T)))

            legend1 = _arr[0][0].legend(handles=legend_handles, loc=9, bbox_to_anchor=(0.5,1.16),fontsize=24,handletextpad=0.4,handlelength=1.0,frameon=True,ncol=len(self.ref_temp)+1,columnspacing=0.7)
            _arr[0][0].add_artist(legend1)

            legend2_handles = []
            legend2_handles.append(Line2D([0],[0],marker='d',markersize=12,markerfacecolor='None',markeredgecolor='k', linestyle='None',label=r'P$_{\text{C1}}$ plane'))
            legend2_handles.append(Line2D([0],[0],marker='o',markersize=12,markerfacecolor='None',markeredgecolor='k', linestyle='None',label=r'P$_{\text{C2}}$ plane'))
            legend2 = _arr[0][0].legend(handles=legend2_handles, loc=1,bbox_to_anchor=(1.0,1.0),fontsize=24, handletextpad=0.4,handlelength=1.4,frameon=True,ncol=1,labelspacing=0.1)
            _arr[0][0].add_artist(legend2)

        else:
            print('Complete legend not implemented yet')


#        plt.show()
#        plt.subplots_adjust(hspace=0.05, top=0.95) 

        if self.split:
            fig.align_ylabels()
            fig2.align_ylabels()
            self.save_figure_split(fig,fig2)
        else:
            fig.align_ylabels()
            if self.diagram_together:
                _arr[0,0].text(-0.5,460*350/500, r'a)',fontsize=28)
                self.fig = fig
                self.arr = _arr[0,1]
            else:
                self.save_figure(fig)

        self.write_full_gap_loc()


    def get_rashba_info(self, p):

        # Get Rashba energy and Rashba parameter alpha for each P,T

        if self.pressure[p] < self.crit_pressure:
            a_idx = self.find_kpt_index([0.0, 0.0, 0.5])
        else:
            a_idx = self.find_kpt_index([0.0,0.0,0.53849])

        energy = self.full_gap_energy_band[p,:,:]

        return energy

    def write_full_gap_loc(self):

        # write full EPI+TE indirect gap info to netcdf file
        fname = 'total_gap_PT_loc.nc'

        with nc.Dataset(fname,'w') as ncdata:

            ncdata.createDimension('number_of_temperatures',self.ntemp)
            ncdata.createDimension('number_of_pressures', len(self.pressure))
            ncdata.createDimension('two',2)
            ncdata.createDimension('three', 3)

            data = ncdata.createVariable('temperature','d',('number_of_temperatures'))
            data[:] = self.ref_temp[:]
            data.units = 'K'

            data = ncdata.createVariable('pressure','d',('number_of_pressures'))
            data[:] = self.pressure[:]
            data.units = 'GPa'

#np.zeros((file_qty, self.ntemp+1,2)) , self.kpoints
            data = ncdata.createVariable('total_bare_gap_location','d',('number_of_pressures','two', 'three'))
            data[:,0,:] = self.gap_loc_coordinates[:,0,0,:]       
            data[:,1,:] = self.gap_loc_coordinates[:,0,1,:] 
                
            data = ncdata.createVariable('total_gap_location_temperature_dependent','d',('number_of_pressures','number_of_temperatures','two','three'))
            data[:,:,0,:] = self.gap_loc_coordinates[:,1:,0,:]
            data[:,:,1,:] = self.gap_loc_coordinates[:,1:,1,:]

        #    data = ncdata.createVariable('bare_bandgap_energy','d',('number_of_pressures'))
        #    data[:] = self.full_gap_energy[:,0#]
#            data.units = self.units

            data = ncdata.createVariable('total_bandgap_energy_temperature_dependent','d', ('number_of_pressures','number_of_temperatures'))
            data[:,:] = self.full_gap_energy[:,:]
            data.units = self.units

            data = ncdata.createVariable('bare_rashba_energy','d',('number_of_pressures','two'))
            data[:,:] = self.rashba_energy[:,0,:]
            data.units = self.units

            data = ncdata.createVariable('rashba_energy','d',('number_of_pressures','number_of_temperatures','two'))
            data[:,:,:] = self.rashba_energy[:,1:,:]
            data.units = self.units



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
 

    def plot_phase_diagram(self):

        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        # Plot (P,T) phase diagram, from output of plot_pgap_indirect

        if self.diagram_together:
            _fig = self.fig
            _arr = self.arr
        else:
            _fig, _arr = plt.subplots(1,1,figsize=self.figsize,squeeze=True)

        _arr.set_rasterization_zorder(0)
        _arr.plot(self.pc1,self.ref_temp,marker='o',markersize=6,color='k',linewidth=1.5)
        _arr.plot(self.pc2,self.ref_temp,marker='o',markersize=6,color='k',linewidth=1.5)
        if not self.crossover:
            _arr.fill_betweenx(self.ref_temp,self.pc1,self.pc2,color='gray',alpha=0.2,zorder=-1, rasterized=True)
    #        _arr.fill_betweenx(self.ref_temp,self.pc2,self.pressure[-1]*np.ones(len(self.pc2)),color='#363256',alpha=0.6,zorder=-1, rasterized=True) #darkbluepurple
            _arr.fill_betweenx(self.ref_temp,self.pc2,self.pressure[-1]*np.ones(len(self.pc2)),color='#6C70C8',alpha=0.6,zorder=-1, rasterized=True) #periwinkle


            legend_handles = [Line2D([0],[0],color='k',linewidth=1,linestyle='dotted',dashes=(2,4),label=r'Static')]
            legend1 = _arr.legend(handles=legend_handles, loc=3,fontsize=24, handletextpad=0.4,handlelength=1.4,frameon=True,ncol = 1)
            _arr.add_artist(legend1)

            _arr.text(0.8,400,r'$\mathbb{Z}_2=0$',fontsize=32)
            _arr.text(3.7,400,r'$\mathbb{Z}_2=1$',fontsize=32)
            _arr.text(2.3,400,r'WSM',fontsize=26,color='#5A5A5A')


        else:

            dim = 1001
            # Read delta_pc from file
            if not self.deltap_fname:
                raise Exception('A file containing delta_p must be provided for crossover PT diagram.')

#            self.deltap_fname = '/Users/Veronique/Google_Drive/doctorat/work/TI/BiTeI/NC/phonons/0gpa/11-TE/delta_p_from_kbt.nc'
            with nc.Dataset(self.deltap_fname,'r') as ncdata:
                delta_pc1 = ncdata.variables['delta_p_trivial_phase'][:]
                delta_pc2 = ncdata.variables['delta_p_topological_phase'][:]

            b1 = np.zeros(len(delta_pc1))
            b2 = np.zeros(len(delta_pc2))

            inc = True
            if inc:
                inc1 = [1.0,1.00,1.00,1.00,1.01,1.03]
                inc2 = [1.0,1.00,1.00,1.20,1.50,1.8]

                for i in range(len(delta_pc1)):
                    b1[i] = inc1[i]*delta_pc1[i]
                    b2[i] = inc2[i]*delta_pc2[i]

                self.delta_pc1 = self.pc1-b1
                self.delta_pc2 = self.pc2+b2

            else:
                self.delta_pc1 = self.pc1-delta_pc1
                self.delta_pc2 = self.pc2+delta_pc2

            self.delta_pc1_plus = self.pc1+delta_pc1
            self.delta_pc2_minus = self.pc2-delta_pc2

            #Define meshgrid
            gridx, gridy = np.meshgrid(np.linspace(self.pressure[0], self.pressure[-1], dim), np.linspace(self.temp[0], self.temp[-1], dim))
            # pc1 and pc2 are still stored
            #_arr.plot(gridx,gridy, marker='.', color='k', linestyle='none')
#            _arr.plot(self.delta_pc1, self.ref_temp, marker='s',color='k', linestyle = 'dashed',linewidth=1.5)
#            _arr.plot(self.delta_pc2, self.ref_temp, marker='s',color='k', linestyle = 'dashed',linewidth=1.5)

            # define a value for each point in the grid, depending on its location vs the phase boundaries

#            _arr.plot(self.pc1-delta_pc1, self.temp, 'r',linewidth=2)
#            _arr.plot(self.pc2+delta_pc2, self.temp, 'r',linewidth=2)
            self.phase = self.get_phase_grid(gridx,gridy,_arr)

#            cmap = mcolors.LinearSegmentedColormap.from_list("", ["cornsilk","palegreen","mediumorchid"])
            if self.cmap_cols:
                cmap = mcolors.LinearSegmentedColormap.from_list("", self.cmap_cols)
            else:
                cmap = 'Blues'

            _arr.imshow(self.phase, extent = [self.pressure[0],self.pressure[-1], self.temp[0], self.temp[-1]], origin='lower', 
                    aspect = 0.65*self.pressure[-1]/self.temp[-1],interpolation='bicubic', vmin=0, vmax=2,cmap=cmap,zorder=1, rasterized=False)            
                                                                    # spline16, bilinear, bicubic

            legend_handles = [Line2D([0],[0],color=self.cmap_cols[0],linestyle='None',marker='s',markersize=20, label=r'$\mathbb{Z}_2=0$',markeredgecolor='k'),
                                Line2D([0],[0],color=self.cmap_cols[1],linestyle='None',marker='s',markersize=20, label=r'WSM',markeredgecolor='k'),
                                Line2D([0],[0],color=self.cmap_cols[2],linestyle='None',marker='s',markersize=20, label=r'$\mathbb{Z}_2=1$',markeredgecolor='k'),
                                Line2D([0],[0],color='k',linestyle='dotted',linewidth=2.5, dashes=(2,4),label=r'Static'),
                                Line2D([0],[0],color='k', linestyle='solid',linewidth=2.5, label=r'P$_{\text{C}1}$,P$_{\text{C2}}$'),
]

            legend2 = _arr.legend(handles = legend_handles, loc=9, bbox_to_anchor = (0.5,1.16), ncol = len(legend_handles),
                            fontsize=24,handletextpad=0.4,frameon=True,handlelength=1.4,columnspacing=1.0)
            _arr.add_artist(legend2)


        #_arr.text(0.5,400,r'T crossover',fontsize=24,zorder=10)
        x = np.ones((len(self.ref_temp)))
        _arr.plot(self.crit_pressure*x, self.ref_temp,'k:',linewidth=1.5, dashes=(2,4))
        _arr.plot(self.crit_pressure2*x, self.ref_temp,'k:',linewidth=1.5,dashes=(2,4))
        

        self.set_xaxis(_arr, self.pressure)
        _arr.set_xlim(self.pressure[0],self.pressure[-1])
        self.set_temp_yaxis(_arr, r'Temperature (K)',self.ref_temp)



        if self.pgap:
            _arr.set_title('EPI',fontsize=24)
        if self.te_pgap:
            _arr.set_title('TE',fontsize=24)
        if self.total_pgap:
            if not self.crossover:
                _arr.set_title('EPI+TE', fontsize=24)
        if self.diagram_together:
            _arr.text(-0.5,460, r'b)',fontsize=28)

#        _inset = inset_axes(_arr, loc=2,width='30%', height='30%',borderpad=3)
#        _inset.plot(self.ref_temp,self.pc2-self.pc1,'b')

#        _arr.set_xlims(

        self.save_phase_diagram(_fig)

        self.write_phase_diagram()

    def get_phase_grid(self, x,y ,f):

        dim = x.shape[0]
        res = np.zeros((x.shape[0], y.shape[1]))
        # Define interpolated lines of pc1-delta_pc1, pc1, pc2, pc2+delta_pc2 as P(T) with interp1D
        # Store lines in a long array of lenght dim for each boundary
        bound1m = np.zeros((dim))
        bound1 = np.zeros((dim))
        bound1p = np.zeros((dim))
        bound2m = np.zeros((dim))
        bound2 = np.zeros((dim))
        bound2p = np.zeros((dim))


        for i in range(dim):
            t = y[i,i]
            loc = self.locate_temp_index(t)

            if loc<len(self.temp)-1:
                line = interp1d(self.temp[loc:loc+2],self.delta_pc1[loc:loc+2], kind='slinear')
                bound1m[i] = line(t)

                line = interp1d(self.temp[loc:loc+2],self.pc1[loc:loc+2], kind='slinear')
                bound1[i] = line(t)

                line = interp1d(self.temp[loc:loc+2],self.delta_pc1_plus[loc:loc+2], kind='slinear')
                bound1p[i] = line(t)

                line = interp1d(self.temp[loc:loc+2],self.delta_pc2_minus[loc:loc+2], kind='slinear')
                bound2m[i] = line(t)

                line = interp1d(self.temp[loc:loc+2],self.pc2[loc:loc+2], kind='slinear')
                bound2[i] = line(t)

                line = interp1d(self.temp[loc:loc+2],self.delta_pc2[loc:loc+2], kind='slinear')
                bound2p[i] = line(t)

            else:
                bound1m[i] = self.delta_pc1[-1]
                bound1[i] = self.pc1[-1]
                bound1p[i] = self.delta_pc1_plus[-1]
                bound2m[i] = self.delta_pc2_minus[-1]
                bound2[i] = self.pc2[-1]
                bound2p[i] = self.delta_pc2[-1]

#        f.plot(bound1m, y[:,0], 'r',linewidth=2)
#        f.plot(bound2, y[:,0], 'b:',linewidth=2)
#        f.plot(bound3, y[:,0], 'c:',linewidth=2)
#        f.plot(bound2p, y[:,0], 'm',linewidth=2)

#        bound1m = 0.9*bound1m
#        bound2p = 1.1*bound2p

        for i,j in itt.product(range(dim),range(dim)):
            '''ligne horizontale'''
#            f.plot(x[i,j], y[0,0], 'ob') 
            ''' ligne verticale'''
#            f.plot(x[11,11],y[i,j], 'or')
            '''point individuel'''
#            f.plot(x[20,20], y[40,40],'og')

            #Define individual point
            point = [x[i,j],y[i,j]]
#            print(i,j,point[0], point[1])

            #For a given Temp value point[1], check where is P value point[0] is located, using boundary arrays
            if point[0] <= bound1m[i]:
                # If trivial, assign value = 0
                res[i,j] = 0 

            if point[0] >= bound1[i] and point[0] <= bound2[i]:
                # If WSM (between pc1 and pc2) assign 0 : for single crossover on P axis
                res[i,j] = 1

                # for double crossover on P and T axis
                #linear
                res[i,j] =  1- 0.3*(point[1]-self.temp[0])/(self.temp[-1]-self.temp[0])
                #gaussian
#                res[i,j] = np.exp(-(point[1]-self.temp[0])**2/(300**2))

            # for double crossover on P axis
#            if point[0] >= bound1p[i] and point[0] <= bound2m[i]:
#                # If WSM (between pc1_plus and pc2_minus) assign 1
#                res[i,j] = 1



            if point[0] >= bound2p[i]: 
                # If topological, assign 2
                res[i,j] = 2

            # If within uncertaintoes from kbT, define value from linear interpolation between boundaries at T on the P axis
            #between pc1_m and pc1: for single crossover on P axis
            if point[0] > bound1m[i] and point[0] < bound1[i]:
                res[i,j] = (point[0]-bound1m[i])/(bound1[i]-bound1m[i])
                #res[i,j] = 0.1 + 0.9*((point[0]-bound1[i]))/(bound2[i]-bound1[i])

                # for double crossover on P and T axis
                #linear
                res[i,j] = res[i,j]*(1 - 0.3*(point[1]-self.temp[0])/(self.temp[-1]-self.temp[0]))
                #gaussian
#                res[i,j] = res[i,j]*np.exp(-(point[1]-self.temp[0])**2/(300**2))

            #between pc1_m and pc1_p : for double crossover on P axis
#            if point[0] > bound1m[i] and point[0] < bound1p[i]:
#                res[i,j] = (point[0]-bound1m[i])/(bound1p[i]-bound1m[i])
                #res[i,j] = 0.1 + 0.9*((point[0]-bound1[i]))/(bound2[i]-bound1[i])


            # between pc2 and pc2_plus : for single crossover on P axis
            if point[0] > bound2[i] and point[0] < bound2p[i]:
                res[i,j] = 1 + (point[0]-bound2[i])/(bound2p[i]-bound2[i])
#                res[i,j] = 1 + 0.85*(point[0]-bound3[i])/(bound4[i]-bound3[i])

                # for double crossover on P and T axis
                #linear renormalized (not really good)
#                res[i,j] = (point[0]-bound2[i])/(bound2p[i]-bound2[i]) +  (1 - 0.3*(point[1]-self.temp[0])/(self.temp[-1]-self.temp[0]))
                #gaussian
#                res[i,j] = res[i,j]*(1+np.exp(-(point[1]-self.temp[0])**2/(300**2)))

                #this one is the good one
                res[i,j] = res[i,j]*(1 + 0.2*(point[1]-self.temp[0])/(self.temp[-1]-self.temp[0]))


#            # between pc2_minus and pc2_plus : for double crossover on P axis
#            if point[0] > bound2m[i] and point[0] < bound2p[i]:
#                res[i,j] = 1 + (point[0]-bound2m[i])/(bound2p[i]-bound2m[i])
#                res[i,j] = 1 + 0.85*(point[0]-bound3[i])/(bound4[i]-bound3[i])




        #self.phase a une shape de (2,dim,dim), comme meshgrid)

#            if self.phase[:,i,j] == -1:
#                f.plot(x[i,j], y[i,j],'ro')
#            if self.phase[:,
        return res

    def write_phase_diagram(self):

        outfile = '{}_phase_diagram.dat'.format(self.rootname)

        create_directory(outfile)

        with open(outfile,'w') as f:

            f.write('Topological phase transition boundaries\n\n')
            f.write('{:>15}  {:>9}  {:>9}  {:>11}  {:>11}\n'.format('Temperature (K)','Pc1 (GPa)','Pc2 (GPa)','width (GPa)','width vs T=0(%)'))
            f.write('{:>15}  {:>9.5f}  {:>9.5f}  {:>11.5f}  {:>11.5f}\n'.format('Static',self.crit_pressure,self.crit_pressure2,self.crit_pressure2-self.crit_pressure,1.0))
            for t, T in enumerate(self.ref_temp):
                f.write('{:>15.2f}  {:>9.5f}  {:>9.5f}  {:>11.5f}  {:>11.5f}\n'.format(T,self.pc1[t],self.pc2[t],self.pc2[t]-self.pc1[t], (self.pc2[t]-self.pc1[t])/(self.crit_pressure2-self.crit_pressure)))

        f.close()
        return

    def plot_self_energy(self):
        
        plot_qty = len(self.bands_to_print)
        file_qty = len(self.zpr_fnames)

        if file_qty > len(self.color):
            raise Exception('color array must have at least {} entries, but it has only {}'.format(file_qty,len(self.color)))

        if self.units == 'eV':
            self.fermi = self.fermi*cst.ha_to_ev
        elif self.units == 'meV':
            self.fermi = self.fermi*1000

        fig, arr = plt.subplots( plot_qty,1, squeeze = False, sharex = True,  figsize = self.figsize)

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
    
                arr[iplot,0].plot(self.omega_se, self.self_energy[0,self.point_for_se, self.bands_to_print[iplot]-1, :, 0], color=self.color[ifile],label=self.labels[ifile], linewidth=1.5)

                self.set_vrefs(arr[iplot,0],self.omega_se, 0.)


            self.set_xaxis(arr[plot_qty-1,0],self.omega_se)
            #    # Add the x=y reference line, which will indicate the fully dynamical renormalized eigenvalue. See Antonius PRB 2016, Fig.1
            #    x = np.linspace(self.xlims[0],self.xlims[1],50)
            #    arr[iplot,0].plot(x,x,color='black',linestyle='dashed')

        if not self.ylims:

            self.ylims = arr[iplot,0].get_ylim()


        # Add the x=y reference line, which will indicate the fully dynamical renormalized eigenvalue. See Antonius PRB 2016, Fig.1
        for iplot in range(plot_qty):
            x = np.linspace(self.xlims[0],self.xlims[1],50)
            arr[iplot,0].plot(x,x,color='black',linestyle='dashed')

            if self.title is not None:
                self.set_title(arr[iplot,0],self.title[iplot])



            # Set axis properties
            self.set_yaxis(arr[iplot,0],'Re[Self-energy] ({})'.format(self.units))

            self.set_hrefs(self.ylims, arr[iplot,0],0.,'black')

        self.set_legend_se(arr[0,0])

#        if plot_qty == 2:
#            self.set_title(arr[0,0],'Valence band')
#            self.set_title(arr[0,1], 'Conduction band')

        self.set_main_title(fig)
        self.save_figure(fig)

        plt.show()

    def plot_spectral_function(self):
        print('spectral')

    def find_gap_index(self, loc):

        loclst = list(loc)
        lst = [list(x) for x in self.kpoints]

        #index = lst.index(loclst)

        #if index is not None:
        #    return index
        #else:
        for k,kpoint in enumerate(lst):
            index = np.allclose(loclst,kpoint)
            if index==True:
        #        print(index,k)
                return k


    def find_temp_index(self):

        lst = list(self.temp) 

        if self.vbcb or self.split_occupied_subspace:
            if self.temp_to_print in lst:
                return lst.index(self.temp_to_print)
            else:
                return None
        if self.renormalization:
            if self.tmax in lst:
                return lst.index(self.tmax)
            else:
                return None

    def locate_temp_index(self, t):
        lst = list(self.temp)

        if t in lst:
            return lst.index(t)
        else:
            for i in range(len(self.temp)-1):
                if self.temp[i] < t and t < self.temp[i+1]:
                    return i
            print('{} not found in temperature array'.format(t))
            return None

    def find_crit_pressure(self, crit_pressure, pressure):
        if self.pgap or self.te_pgap or self.total_pgap:
            lst = list(pressure)
            if crit_pressure in lst:
                return lst.index(crit_pressure)
            else:
                for i in range(len(pressure)-1):
                    if pressure[i]<crit_pressure and crit_pressure < pressure[i+1]:
                        return i
                return None

    def find_crit_index(self):

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

        if self.explicit_pressures is not None:
            pressure = self.explicit_pressures
        else:
            pressure = self.pressure

        lst = list(pressure) 
        if p in lst:
            return lst.index(p)
        else:
            for i in range(len(pressure)-1):
                if pressure[i]<p and p< pressure[i+1]:
#                    print("in between")
                    return i
            return None


    def get_gap_info(self, unpert,pert,split):

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
            for t in range(self.ntemp):
                tmp[:,:,:,t] = self.interpolate_eigenvalues(pert[:,:,:,t])

            pert = tmp

        if split == False:

            loc = np.zeros((self.ntemp+1),dtype=int)
            ener = np.zeros((self.ntemp+1), dtype=float)
            ren = np.zeros((self.ntemp), dtype=float)
        
            arr0 = unpert[0,:,cond] - unpert[0,:,val]
            loc[0] = np.argmin(arr0)
            ener[0] = np.amin(arr0)
              
            for i in range(self.ntemp):
           
                arr = pert[0,:,cond,i]-pert[0,:,val,i]
                loc[i+1] = np.argmin(arr)
                ener[i+1] = np.amin(arr)
                ren[i] = ener[i+1]-ener[0]

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
            ener[0,0] = np.amin(arr0)
            diff = (unpert[:,:kmid,28]-unpert[:,:kmid,27])

             
            for i in range(self.ntemp):
                arr = pert[0,:kmid,cond,i]-pert[0,:kmid,val,i]
                loc[i+1,0] = np.argmin(arr)
                ener[i+1,0] = np.amin(arr) 
                ren[i,0] = ener[i+1,0]-ener[0,0]
       
            #Treat right gap
            arr0 = unpert[0,kmid:,cond] - unpert[0,kmid:,val]
            loc[0,1] = np.argmin(arr0)+kmid
            ener[0,1] = np.amin(arr0)

            
            for i in range(self.ntemp):
                arr = pert[0,kmid:,cond,i]-pert[0,kmid:,val,i]
                loc[i+1,1] = np.argmin(arr)+kmid
                ener[i+1,1] = np.amin(arr) 
                ren[i,1] = ener[i+1,1]-ener[0,1]


            return loc, ener, ren

    def get_indirect_gap_info(self,unpert,pert,split):
        # Gets indirect gap information

        # Finds the location of the unperturbed and perturbed direct gap + computes the gap energy
        val = self.valence - 1
        cond = val + 1
        #print(val, cond)

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

            loc = np.zeros((self.ntemp+1,2),dtype=int)
            ener = np.zeros((self.ntemp+1), dtype=float)
            ener_band = np.zeros((self.ntemp+1,2), dtype=float)
            ren = np.zeros((self.ntemp), dtype=float)
            ren_band = np.zeros((self.ntemp,2), dtype=float)

            arrc = unpert[0,:,cond]
            arrv = unpert[0,:,val]
            loc[0,0] = np.argmax(arrv)
            loc[0,1] = np.argmin(arrc)
            ener_band[0,0] = arrv[loc[0,0]]
            ener_band[0,1] = arrc[loc[0,1]] 
            ener[0] = ener_band[0,1] - ener_band[0,0]
              
            for i in range(self.ntemp):
           
                arrc = pert[0,:,cond,i]
                arrv = pert[0,:,val,i]
                loc[i+1,0] = int(np.argmax(arrv))
                loc[i+1,1] = int(np.argmin(arrc))
                ener_band[i+1,0] = arrv[loc[i+1,0]] 
                ener_band[i+1,1] = arrc[loc[i+1,1]] 
                ener[i+1] = ener_band[i+1,1] - ener_band[i+1,0] 
                for j in range(2):
                    ren_band[i,j] = ener_band[i+1,j] - ener_band[0,j]
                ren[i] = ener[i+1]-ener[0]

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


            #Treat left gap
            arrc = unpert[0,:kmid,cond]
            arrv = unpert[0,:kmid,val]
            loc[0,0,0] = np.argmax(arrv)
            loc[0,0,1] = np.argmin(arrc)
            ener_band[0,0,0] = arrv[loc[0,0,0]]
            ener_band[0,0,1] = arrc[loc[0,0,1]] 
            ener[0,0] = ener_band[0,0,1] - ener_band[0,0,0]
             
            for i in range(self.ntemp):
                arrc = pert[0,:kmid,cond,i]
                arrv = pert[0,:kmid,val,i]
                loc[i+1,0,0] = np.argmax(arrv)
                loc[i+1,0,1] = np.argmin(arrc)
                ener_band[i+1,0,0] = arrv[loc[i+1,0,0]] 
                ener_band[i+1,0,1] = arrc[loc[i+1,0,1]] 
                ener[i+1,0] = ener_band[i+1,0,1] - ener_band[i+1,0,0] 
                for j in range(2):
                    ren_band[i,0,j] = ener_band[i+1,0,j] - ener_band[0,0,j]
                ren[i,0] = ener[i+1,0]-ener[0,0]


            #Treat right gap
            arrc = unpert[0,kmid:,cond]
            arrv = unpert[0,kmid:,val]
            loc[0,1,0] = np.argmax(arrv)
            loc[0,1,1] = np.argmin(arrc)
            ener_band[0,1,0] = arrv[loc[0,1,0]]
            ener_band[0,1,1] = arrc[loc[0,1,1]] 
            ener[0,1] = ener_band[0,1,1] - ener_band[0,1,0]
            loc[0,1,0] += kmid
            loc[0,1,1] += kmid
            
            for i in range(self.ntemp):
                arrc = pert[0,kmid:,cond,i]
                arrv = pert[0,kmid:,val,i]
                loc[i+1,1,0] = np.argmax(arrv)
                loc[i+1,1,1] = np.argmin(arrc)
                ener_band[i+1,1,0] = arrv[loc[i+1,1,0]] 
                ener_band[i+1,1,1] = arrc[loc[i+1,1,1]] 
                ener[i+1,1] = ener_band[i+1,1,1] - ener_band[i+1,1,0] 
                for j in range(2):
                    ren_band[i,1,j] = ener_band[i+1,1,j] - ener_band[0,1,j]
                ren[i,1] = ener[i+1,1]-ener[0,1]
                loc[i+1,1,0] += kmid
                loc[i+1,1,1] += kmid

            return loc, ener, ren, ener_band, ren_band

    def interpolate_kpoints(self):
        
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

        return

    def interpolate_eigenvalues(self,eig0):

         nkpt = self.nkpt

         x = np.arange(nkpt)
         xfine = np.linspace(x[0],x[-1],self.nkpt_fine)

#         if self.nsppol is None:
#             self.nsppol = eig0.nsppol
#         if self.max_band is None:
#             self.max_band = eig0.max_band

         eig = np.zeros((self.nsppol,self.nkpt_fine,self.max_band))

         for j in range(self.max_band):
             spl = interp1d(x,eig0[0,:,j], kind='cubic')
             eig[0,:,j] = spl(xfine)

         return eig

    def find_kpt_index(self, loc):

        # This locates a given kpoint in a list, within a given tolerance
        loclst = list(loc)
        if self.spline:
            lst = [list(x) for x in self.kpoints_fine]
#        elif self.reduce_path:
#            lst = [list(x) for x in self.reduced_kpath]
        else:
            lst = [list(x) for x in self.kpoints]

        for k,kpoint in enumerate(lst):
            index = np.allclose(loclst,kpoint)
            if index==True:
               # print(kpoint)
                return k

#    def renormalize_0K_data(self):
#        print(np.shape(self.full_band_ren), np.shape(self.full_gap_ren))
        #self.full_band_ren = self.full_band_ren - 
    #full_band_ren, full_gap_ren

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
            cut = 0.02*max(self.temp)
            lims = (-5, 1.01*max(self.temp))
            f.set_xlim(lims)
            f.set_xlabel('Temperature (K)', fontsize=32)
            plt.setp(f.get_xticklabels(), fontsize=26, weight='bold')
#            f.xaxis.set_tick_params(labelsize=16)
#            if not self.xlims:
#                self.xlims = [min(x), max(x)]
#            f.set_xlim(self.xlims)

        if self.vbcb:
            f.set_xlim((0,len(x)-1))
            if self.xticks is not None:
                plt.tick_params(axis='x',labelsize = 18)
#                plt.setp(f.get_xticklabels(), fontsize='large')
                plt.setp(f,xticks=self.xticks[0], xticklabels=self.xticks[1])
    
                if self.xticks_alignment is not None:
                    for itick, tick in enumerate(f.xaxis.get_majorticklabels()):
                        tick.set_horizontalalignment(self.xticks_alignment[itick]) 

        if self.pgap or self.split_occupied_subspace or self.te_pgap or self.total_pgap:
            lims = (min(self.pressure)-0.1, max(self.pressure)+0.1)
            f.set_xlim(lims)
            f.xaxis.set_major_formatter(FuncFormatter(self.label_formatter))

            plt.setp(f.get_xticklabels(), fontsize=26, weight='bold')
            if self.separate_bands:
                f.set_xlabel('Pressure (GPa)', fontsize=32)
            else:
                f.set_xlabel('Pressure (GPa)',fontsize=32)

    def set_temp_yaxis(self,f,lab,temp):

        f.set_ylabel(lab, fontsize = 32)

        f.set_ylim(temp[0],temp[-1])

        f.yaxis.set_major_formatter(FuncFormatter(self.label_formatter))

        plt.setp(f.get_yticklabels(), fontsize=26, weight='bold')
#        if self.yminorticks:
#            f.yaxis.set_minor_locator(AutoMinorLocator(self.yminorticks+1))



    def set_yaxis(self,f, lab):

        if self.vbcb:
            f.set_ylabel(lab, fontsize=14)
            
        else:
            f.set_ylabel(lab, fontsize = 32)


            if self.ylims is not None:
                f.set_ylim(self.ylims) 
            else:
                lims = f.get_ylim()
                f.set_ylim(lims)

        f.yaxis.set_major_formatter(FuncFormatter(self.label_formatter))

        if self.yticks is not None:
            plt.setp(f.get_yticklabels(), fontsize=26)
            f.set_yticks(np.arange(self.yticks[0], self.yticks[1]+self.yticks[2], self.yticks[2]))
        else:
            plt.setp(f.get_yticklabels(), fontsize=26, weight='bold')
        if self.yminorticks:
            f.yaxis.set_minor_locator(AutoMinorLocator(self.yminorticks+1))

    def set_yaxis_separate(self,f, lab,lims=None,ticks=None):

        f.set_ylabel(lab, fontsize = 32)
        
        if lims is not None:
            f.set_ylim(lims) 
        else:
            lims = f.get_ylim()
            f.set_ylim(lims)

        f.yaxis.set_major_formatter(FuncFormatter(self.label_formatter))

        if ticks is not None:
            plt.setp(f.get_yticklabels(), fontsize=26)
            f.set_yticks(np.arange(ticks[0], ticks[1]+ticks[2], ticks[2]))
        else:
            plt.setp(f.get_yticklabels(), fontsize=26, weight='bold')
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

    def label_formatter(self,x, pos):
        return "%i" %x

    def set_hrefs(self, lims, f, val,col,style='solid'):

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
        if style=='dashed':
            f.plot(zer,y,color=col,linestyle=style, dashes=(4,4))
        else:
            f.plot(zer,y,color=col,linestyle=style)


    def set_vrefs(self,f,x,val,style='solid'):

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
        if style == 'dashed':
            f.plot(x, zer, color='black', linestyle=style,dashes=(4,4))
        else:
            f.plot(x, zer, color='black', linestyle=style)


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

        if self.labels:
            f.legend(loc=1)

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
                t = self.main_title+', Kpt : '+str(self.kpoints[self.point_for_se])
            else:
                t = self.main_title
            g.suptitle(t, fontsize=18)

    def save_figure(self,g):

        create_directory('figures/')

        if self.savefile:
            if self.save_format == 'eps' or self.save_format == 'pdf':
                g.savefig('figures/{}.{}'.format(self.savefile,self.save_format), format=self.save_format,dpi=1200)
                os.system('open figures/{}.{}'.format(self.savefile,self.save_format))
            else:
                g.savefig('figures/{}.{}'.format(self.savefile,self.save_format))
                os.system('open figures/{}.{}'.format(self.savefile,self.save_format))

        else:
            plt.show()
    
    def save_phase_diagram(self,g):

        create_directory('figures/')
        if self.savefile:
            if self.save_format == 'eps' or self.save_format == 'pdf':
                outfile = 'figures/{}_PTphase.{}'.format(self.savefile,self.save_format)
                g.savefig(outfile,format='{}'.format(self.save_format),dpi=1200)
                os.system('open {}'.format(outfile))
            else:
                outfile = 'figures/{}_PTphase.{}'.format(self.savefile,self.save_format)
                g.savefig(outfile)
                os.system('open {}'.format(outfile))
        else:
            plt.show()


    def save_figure2(self,g):

        create_directory('figures/')

        g.savefig('figures/unperturbed_energy.png')

    def save_figure_split(self,g1,g2):

        create_directory('figures/')

        if self.savefile:
            g1.savefig('figures/{}_lgap.png'.format(self.savefile))
            g2.savefig('figures/{}_rgap.png'.format(self.savefile))

    def save_figure_subspaces(self,g1,g2,g3):

        create_directory('figures/')

        if self.savefile:
            g1.savefig('figures/{}.png'.format(self.savefile))
            g2.savefig('figures/{}_Fan-DW.png'.format(self.savefile))
            g3.savefig('figures/{}_occ-unocc.png'.format(self.savefile))

    
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
        te_fnames = None,
        bare_fnames = None,
        deltap_fname = None,
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
        crit_pressure2 = None,

        # Options
        figsize = (15,9),
        color = None,
        labels = None,
        linestyle = None,
        ylims = None,
        yticks = None,
        cond_ylims = None,
        val_ylims = None,
        gap_ylims = None,
        egap_ylims = None,

        cond_yticks = None,
        val_yticks = None,
        gap_yticks = None,

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
        te_pgap = False,
        pgap = False,
        total_pgap = False,
        vbcb = False,
        separate_bands = False,
        phase_diagram = False,
        split = False,
        split2 = False,
        follow = False,
        indirect = False,
        unpert = False,
        split_contribution = False,
        split_occupied_subspace = False,
        modes = False,
        verbose = False,
        zero_gap_value = None,
        zero_gap_units = 'eV',
        experimental_data = None,
        valence = None,
        spline = False,
        spline_factor = 5,
        extrapolate_ren = False,
        read_correction_from_spline = False,
        cmap_cols = None,
        crossover = False,
        diagram_together = False,
        save_format = 'png',

        **kwargs):

        zpr_plot = ZPR_plotter(
            zpr_fnames = zpr_fnames,
            rootname = rootname,
            gap_fname = gap_fname,
            te_fnames = te_fnames,
            bare_fnames = bare_fnames,
            deltap_fname = deltap_fname,

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
            crit_pressure2 = crit_pressure2,


            xlims = xlims,
            xticks = xticks,
            xticks_alignment = xticks_alignment,
            ylims = ylims,
            yticks = yticks,
            cond_ylims = cond_ylims,
            val_ylims = val_ylims,
            gap_ylims = gap_ylims,
            egap_ylims = egap_ylims,

            cond_yticks = cond_yticks,
            val_yticks = val_yticks,
            gap_yticks = gap_yticks,

            yminorticks = yminorticks,
            figsize = figsize,
            cmap_cols = cmap_cols,

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
            te_pgap = te_pgap,
            total_pgap = total_pgap,
            separate_bands = separate_bands,
            phase_diagram = phase_diagram,
            split = split,
            split2 = split2,
            follow = follow,
            indirect = indirect,
            unpert = unpert,
            split_contribution = split_contribution,
            split_occupied_subspace = split_occupied_subspace,
            modes = modes,
            verbose = verbose,
            valence = valence,

            spline = spline,
            spline_factor = spline_factor,
            extrapolate_ren = extrapolate_ren,
            read_correction_from_spline = read_correction_from_spline,
            crossover = crossover,
            diagram_together = diagram_together,
            save_format = save_format,

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

                    if phase_diagram:
                        zpr_plot.plot_phase_diagram()

            else:
                if separate_bands:
                    zpr_plot.plot_pgap_separate()
                else:
                    zpr_plot.plot_pgap()

        if te_pgap:
            if indirect:
                if separate_bands:
                    zpr_plot.plot_te_pgap_separate_indirect()
                else:
                    zpr_plot.plot_te_pgap_indirect()

                    if phase_diagram:
                        zpr_plot.plot_phase_diagram()
            else:
                raise Exception('plot_te_pgap not implemented for direct gap yet')

        if total_pgap:
            if indirect:
                if separate_bands:
                    zpr_plot.plot_total_pgap_separate_indirect()
                else:
                    zpr_plot.plot_total_pgap_indirect()

                    if phase_diagram:
                        zpr_plot.plot_phase_diagram()
            else:
                raise Exception('plot_total_pgap not implemented for direct gap yet')



        if split_occupied_subspace:
            zpr_plot.plot_splitted_subspaces()

        

        return zpr_plot







