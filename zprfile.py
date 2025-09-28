from cdffile import CDFfile
import netCDF4 as nc
import numpy as np


class ZPRfile(CDFfile):

    def __init__(self,  *args,  **kwargs):

        super(ZPRfile,  self).__init__(*args,  **kwargs)

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
        self.fan_g2 = None
        self.zpr_qpath = None

    # Open the ZPR.nc file and read it
    def read_nc(self,  fname=None):

        fname = fname if fname else self.fname
        super(ZPRfile,  self).read_nc(fname)

        with nc.Dataset(fname,  'r') as ncdata:

            self.eig0 = ncdata.variables['reduced_bare_eigenvalues'][:, :, :]
            self.kpoints = ncdata.variables['reduced_coordinates_of_reduced_kpath'][:, :]
            self.qpoints = ncdata.variables['reduced_coordinates_of_qpoints'][:, :]

            self.temp = ncdata.variables['temperatures'][:]
            self.eigcorr = ncdata.variables['reduced_corrected_eigenvalues'][:, :, :, :]
            self.correction = ncdata.variables['reduced_eigenvalue_corrections'][:, :, :, :]
            self.omega_se = ncdata.variables['omega_se'][:]
            self.self_energy = ncdata.variables['reduced_self_energy'][:, :, :, :, :]
            self.spectral_function = ncdata.variables['spectral_function'][:, :, :, :]
            self.broadening = ncdata.variables['broadening'][:, :, :]
            self.smearing = ncdata.variables['smearing'][:]

#            self.gridsize = ncdata.gridsize
            self.unperturbed_gap_location = ncdata.variables['unperturbed_gap_location'][:]
            self.unperturbed_gap_energy = ncdata.variables['unperturbed_gap_energy'][:]
            self.gap_location = ncdata.variables['gap_location'][:, :]
            self.gap_energy = ncdata.variables['gap_energy'][:]
            self.gap_energy_units = ncdata.variables['gap_energy'].getncattr('units')
            self.gap_renormalization = ncdata.variables['gap_renormalization'][:]
#            self.band_energy = ncdata.variables['band_energy'][:, :]
#            self.band_renormalization = ncdata.variables['band_renormalization'][:, :]

            self.eigcorr_modes = ncdata.variables['eigenvalue_corrections_modes'][:, :, :, :]

            self.reduced_zpr_qpt = ncdata.variables['reduced_zpr_qpt_contribution'][:, :, :, :]  # [s,k,n,q]
            self.reduced_zpr_qpt_mode = ncdata.variables['reduced_qpt_contribution_modes'][:, :, :, :, :]  # [s,k,n,q,v]

            self.unperturbed_gap_location_split = ncdata.variables['unperturbed_gap_location_split'][:, :]
#            self.unperturbed_gap_energy_split = ncdata.variables['unperturbed_gap_energy_split'][:, :]
            self.gap_location_split = ncdata.variables['gap_location_split'][:, :, :]
            self.gap_energy_split = ncdata.variables['gap_energy_split'][:, :]
            self.gap_energy_units_split = ncdata.variables['gap_energy_split'].getncattr('units')
            self.gap_renormalization_split = ncdata.variables['gap_renormalization_split'][:, :]
#            self.band_energy_split = ncdata.variables['band_energy_split'][:, :, :]
#            self.band_renormalization_split = ncdata.variables['band_renormalization_split'][:, :, :]

            self.unperturbed_indirect_gap_location = ncdata.variables['unperturbed_indirect_gap_location'][:, :]
            self.indirect_gap_location = ncdata.variables['indirect_gap_location'][:, :, :]
            self.unperturbed_indirect_gap_energy = ncdata.variables['unperturbed_indirect_gap_energy'][:]
            self.indirect_gap_energy = ncdata.variables['indirect_gap_energy'][:]
            self.indirect_gap_ren = ncdata.variables['indirect_gap_renormalization'][:]
            self.indirect_gap_energy_band = ncdata.variables['indirect_gap_energy_band'][:, :]
            self.indirect_gap_ren_band = ncdata.variables['indirect_gap_renormalization_band'][:, :]

            self.unperturbed_indirect_gap_location_split = ncdata.variables['unperturbed_indirect_gap_location_split'][:, :, :]
            self.indirect_gap_location_split = ncdata.variables['indirect_gap_location_split'][:, :, :, :]
            self.unperturbed_indirect_gap_energy_split = ncdata.variables['unperturbed_indirect_gap_energy_split'][:, :]
            self.indirect_gap_energy_split = ncdata.variables['indirect_gap_energy_split'][:, :]
            self.indirect_gap_ren_split = ncdata.variables['indirect_gap_renormalization_split'][:, :]
            self.indirect_gap_energy_band_split = ncdata.variables['indirect_gap_energy_band_split'][:, :, :]
            self.indirect_gap_ren_band_split = ncdata.variables['indirect_gap_renormalization_band_split'][:, :, :]

#            self.fan_occ = ncdata.variables['reduced_fan_occ'][:, :, :, :]
#            self.fan_unocc = ncdata.variables['reduced_fan_unocc'][:, :, :, :]
#            self.ddw_occ = ncdata.variables['reduced_ddw_occ'][:, :, :, :]
#            self.ddw_unocc = ncdata.variables['reduced_ddw_unocc'][:, :, :, :]

            status = ncdata.get_variables_by_attributes(name='corrected_eigenvalue_spline_interpolation')
            if status != []:
                self.eigcorr_spline = ncdata.variables['corrected_eigenvalue_spline_interpolation'][:, :, :, :]
            status = ncdata.get_variables_by_attributes(name='eigenvalue_corrections_spline_interpolation')
            if status != []:
                self.td_ren_spline = ncdata.variables['eigenvalue_corrections_spline_interpolation'][:, :, :, :]
            status = ncdata.get_variables_by_attributes(name='reduced_coordinated_of_kpoints_spline')
            if status != []:
                self.kpoints_spline = ncdata.variables['reduced_coordinated_of_kpoints_spline'][:, :]

            status = ncdata.get_variables_by_attributes(name='qpt_contribution_temperature_dependent')
            if status != []:
                self.reduced_tdr_qpt = ncdata.variables['qpt_contribution_temperature_dependent'][:,:,:,:,:]
            status = ncdata.get_variables_by_attributes(name='qpt_contribution_temperature_dependent_modes')
            if status != []:
                self.reduced_tdr_qpt_mode = ncdata.variables['qpt_contribution_temperature_dependent_modes'][:,:,:,:,:,:]

            # Only for split contribution,  VB and CB / modes separate ###

#            self.fan_g2 = ncdata.variables['reduced_fan_g2'][:, :, :, :, :, :, :] # spin kpt 2 2 mode qpt cplex
#            self.ddw_g2 = ncdata.variables['reduced_ddw_g2'][:, :, :, :, :, :, :] # same
#            self.deltaE_ddw = ncdata.variables['reduced_deltaE_ddw'][:, :, :, :, :, :] # spin kpt 2 2 qpt cplex
#            self.fan_occterm = ncdata.variables['reduced_fan_occterm'][:, :, :, :, :, :, :, :] # spin kpt 2 2 mode temp qpt cplex
#            self.qpt_weight = ncdata.variables['qpt_weight'][:]
#            self.ddw_tdep = ncdata.variables['ddw_tdep'][:, :, :] # mode temp qpt

#            self.zpr_qpath = ncdata.variables['reduced_zpr_qpath'][:,:,:,:,:]
#            self.zpr_qpath_fan = ncdata.variables['reduced_zpr_qpath_fan'][:,:,:,:,:]
#            self.zpr_qpath_ddw = ncdata.variables['reduced_zpr_qpath_ddw'][:,:,:,:,:]
#            self.gkk2_qpath_fan = ncdata.variables['reduced_gkk2_qpath_fan'][:,:,:,:,:]
#            self.gkk2_qpath_ddw = ncdata.variables['reduced_gkk2_qpath_ddw'][:,:,:,:,:]


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
            return self.eigcorr_modes.shape[-1]
        elif self.zpr_qpath is not None:
            return self.zpr_qpath.shape[-1]
        else:
            return None

    @property
    def nqpt(self):
        if self.fan_g2 is not None:
            return self.fan_g2.shape[-2]
        else:
            return self.qpoints.shape[0]
