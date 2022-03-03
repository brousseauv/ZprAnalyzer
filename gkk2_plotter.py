#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable 
from epcfile import EpcGkk2File

plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')

class Gkk2Plotter(object):

    def __init__(self,

            gkk2_flist=None,
            qpt_weight=None,
            kpt_index=1,
            valence=None,
            band_range=None,
            qpt_index=None,
            mode=None,

            single=True,
            diff=False,

            title=None,
            verbose=False,
            figsize=None,
            savefig_root = None,
            diff_labels = None,

            *args, **kwargs):

        # Transfer to Python indexing

        self.kpt_index = kpt_index - 1
        self.valence = valence - 1
        self.band_range = band_range-np.ones_like(band_range)  # transfer to python indexing
        self.qpt_index = qpt_index - 1 
        self.mode = mode - 1

        self.title = title
        self.verbose = verbose
        self.figsize = figsize
        self.savefig_root = savefig_root
        self.diff_labels = diff_labels

        self.single = single
        self.diff = diff

        if self.single and self.diff:
            raise Exception('Must choose only one plotting option, "single" or "diff"')

        if not self.figsize:
            if self.single:
                self.figsize = (8,4)
            elif self.diff:
                self.figsize = (8,8)

        self.gkk2_flist = gkk2_flist
        if not self.gkk2_flist:
            raise Exception('Must provide "*EP.nc" file list, computed with "epi_matrix_elements=True"')

        nfile = len(self.gkk2_flist)
        if self.diff and nfile != 2:
            raise Exception('For "diff" option, gkk2_flist must contain 2 file names, but has only {}'.format(nfile))

        if len(self.band_range) != 2:
            raise Exception('Band_range must be of length 2, but has length {}'.format(len(self.band_range)))

        # Call plotting functions
        if self.single:

            for ifile in range(nfile):

                gkk2 = EpcGkk2File(fname=self.gkk2_flist[ifile])
                gkk2.read_nc()

                if self.mode+1 > gkk2.nmode:
                    raise Exception("""Input asked for mode {} but file {}
                                    has only {} modes""".format(self.mode+1, self.gkk2_flist[ifile], gkk2.nmode))
                if self.band_range[-1]+1 > gkk2.nband:
                    raise Exception("""Input asked for band_max {} but file {}
                                    has only {} bands""".format(self.band_range[-1]+1, self.gkk2_flist[ifile], gkk2.nband))

                fig, ax = plt.subplots(1,2, squeeze=True, figsize=self.figsize)
#                self.imshow_gkk2_single(gkk2)
                self.imshow_gkk2_single(gkk2, fig, ax)

                self.set_titles(ax, gkk2.qpoints[self.qpt_index])
                self.save_and_display(fig)

        elif self.diff:

            gkk2 = EpcGkk2File(fname=self.gkk2_flist[0])
            gkk2.read_nc()

            gkk2b = EpcGkk2File(fname=self.gkk2_flist[1])
            gkk2b.read_nc()

            if self.mode+1 > gkk2.nmode:
                raise Exception("""Input asked for mode {} but file {}
                                has only {} modes""".format(self.mode+1, self.gkk2_flist[ifile], gkk2.nmode))
            if self.band_range[-1]+1 > gkk2.nband:
                raise Exception("""Input asked for band_max {} but file {}
                                has only {} bands""".format(self.band_range[-1]+1, self.gkk2_flist[ifile], gkk2.nband))

            fig, ax = plt.subplots(2,2, squeeze=False, figsize=self.figsize)
            self.imshow_gkk2_single(gkk2, fig, ax[0,:])
            self.imshow_gkk2_difference(gkk2, gkk2b, fig, ax[1,:])
            self.set_titles(ax[0], gkk2.qpoints[self.qpt_index])
            self.save_and_display(fig)

#        plt.show()


    def imshow_gkk2_single(self, data, f, ax):

        # Store the relevant matrix elements in 2D array
        fan, ddw = self.extract_fan_ddw(data)
        print(fan[1,1], fan[1,-1])

        limits = np.array(self.band_range)
        limits = limits + np.array([-0.5, -0.5])

        # Use imshow
        #f, ax = plt.subplots(1,2, squeeze=True, figsize=self.figsize)
        extent=[self.band_range[0], self.band_range[1], self.band_range[1], self.band_range[0]]
        ax0 = ax[0].imshow(fan, origin='upper', extent=extent)
#        ax0.set_clim(gkk_min, gkk_max)
        ax1 = ax[1].imshow(ddw, origin="upper", extent=extent)
#        ax1.set_clim(gkk_min, gkk_max)

        # Add colorbars
        self.set_colorbar(f, ax[0], ax0, title=None)
        self.set_colorbar(f, ax[1], ax1, title='single')

        # Make the figure cute
        #self.set_limits(ax, limits, origin='upper')
        self.show_valence(ax, origin='upper')
        self.set_axes_and_labels(ax)

#        plt.tight_layout(h_pad=1)
#        plt.subplots_adjust(left=0.05, wspace=0.2)
#        if self.savefig_root:
#            fname = 'gkk2_{}_qpt{}_mode{}.png'.format(self.savefig_root, self.qpt_index+1, self.mode+1)
#            plt.savefig(fname)

#        plt.show()

    def imshow_gkk2_difference(self, data1, data2, f, ax):

        # Store the relevant matrix elements in 2D array
        fan, ddw = self.extract_fan_ddw_difference(data1, data2)
        print(fan[1,1], fan[1,-1])

        limits = np.array(self.band_range)
        limits = limits + np.array([-0.5, -0.5])

        # Use imshow
        extent=[self.band_range[0], self.band_range[1], self.band_range[1], self.band_range[0]]

        ax0 = ax[0].imshow(fan, origin="upper", extent=extent)
#        ax0.set_clim(gkk_min, gkk_max)
        ax1 = ax[1].imshow(ddw, origin="upper", extent=extent)
#        ax1.set_clim(gkk_min, gkk_max)

        # Add colorbars
        self.set_colorbar(f, ax[0], ax0, title=None)
        self.set_colorbar(f, ax[1], ax1, title='diff')

        # Make the figure cute
        #self.set_limits(ax, limits, origin='upper')
        self.show_valence(ax, origin='upper')
        self.set_axes_and_labels(ax)

    def extract_fan_ddw_difference(self, data1, data2):
        
        fan1, ddw1 = self.extract_fan_ddw(data1)
        fan2, ddw2 = self.extract_fan_ddw(data2)

#        return np.abs(fan1-fan2), np.abs(ddw1-ddw2)
        return fan1-fan2, ddw1-ddw2


    def extract_fan_ddw(self, data):
        
        v = self.mode
        k = self.kpt_index
        bmin = self.band_range[0]
        bmax = self.band_range[1]+1
        q = self.qpt_index

        fan = data.gkk2_fan[v,k,bmin:bmax,bmin:bmax,q]
        ddw = data.gkk2_ddw[v,k,bmin:bmax,bmin:bmax,q]

#        print('Fan shape', np.shape(fan))
        return fan, ddw

    def save_and_display(self, f):

        f.tight_layout(h_pad=1)
#        plt.subplots_adjust(left=0.05, wspace=0.2)
        if self.savefig_root:
            if self.diff:
                fname = 'diffgkk2_{}_qpt{}_mode{}.png'.format(self.savefig_root, self.qpt_index+1, self.mode+1)
            elif self.single:
                fname = 'gkk2_{}_qpt{}_mode{}.png'.format(self.savefig_root, self.qpt_index+1, self.mode+1)

            f.savefig(fname, bbox_inches='tight')
        plt.show()

    def set_colorbar(self, g, axes, img, title=None):

#        cbar = g.colorbar(myarr0, cax=cbax, orientation='vertical') #, format=FormatStrFormatter("%.1f"))
#        cbax = g.add_axes([0.885, 0.125, 0.015, 0.825])
#        cbax0 = g.add_axes([0.465, 0.195, 0.020, 0.600])
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = g.colorbar(img, cax=cax, orientation='vertical', shrink=0.5, aspect=20, format=FormatStrFormatter("%.1E"))

        if title == 'single':
            cbar.set_label(r"$|g_{knn'}(\mathbf{q}\nu)|^2$", fontsize=16)
        if title == 'diff':
#            cbar.set_label(r"$|\Delta(|g_{knn'}(\mathbf{q}\nu)|^2)|$ (%s"%self.diff_labels[0]+"-%s"%self.diff_labels[1]+")", fontsize=16)
            cbar.set_label(r"$\Delta(|g_{knn'}(\mathbf{q}\nu)|^2)$ (%s"%self.diff_labels[0]+"-%s"%self.diff_labels[1]+")", fontsize=16)


    def show_valence(self, f, origin='upper'):

        # Add vertical lines to show separation between VBM and CBM band indices
        if self.band_range[0] < self.valence < self.band_range[1]:

            # This one is without the 'extent' parameter
            #i = self.valence - self.band_range[0] + 1
            #n = self.band_range[1] + 1 - self.band_range[0]
            #x = [i, i]
            #y = [-1, n]

            x = [self.valence, self.valence]
            x2 = [self.valence-4, self.valence-4] # those are the relevant VBM degenerate bands
            x3 = [self.valence-6, self.valence-6] # those are the splitoff

            if origin == 'lower':
                y = [self.band_range[0], self.band_range[1]]

            elif origin == 'upper':
                y = [self.band_range[1], self.band_range[0]]


            for j in range(2):
                f[j].plot(x, y, 'w')
                f[j].plot(y, x, 'w')
                # show lower VBM degen
                f[j].plot(x2, y, 'c:')
                f[j].plot(y, x2, 'c:')
                # Show splitoff
                f[j].plot(x3, y, 'm:')
                f[j].plot(y, x3, 'm:')

        return

    def set_titles(self, f, qpt):

#        if self.single:
#            f[0].set_title(r"$|g_{knn'}^{\rm{Fan}}|^2$", fontsize=16, pad=10)
#            f[1].set_title(r"$|g_{knn'}^{\rm{DW}}|^2$", fontsize=16, pad=10)
        f[0].set_title(r"Fan", fontsize=16, pad=10)
        f[1].set_title(r"DW", fontsize=16, pad=10)


#        if self.diff:
#            f[0].set_title(r"$|\Delta(|g_{knn'}^{\rm{Fan}}|^2)|$ (%s"%self.diff_labels[0]+"-%s"%self.diff_labels[1]+")", fontsize=16, pad=10)
#            f[1].set_title(r"$|\Delta(|g_{knn'}^{\rm{DW}}|^2)|$ (%s"%self.diff_labels[0]+"-%s"%self.diff_labels[1]+")", fontsize=16, pad=10)

        plt.suptitle(r'$\mathbf{q}$ = (%.3f'%qpt[0]+r', %.3f'%qpt[1]+r', %.3f'%qpt[2]+r'), $\nu$ = %i'%(self.mode+1), fontsize=16)

    def set_axes_and_labels(self, f):

        for j in range(2):
#            f[j].set_xlabel(r"n $\rightarrow$", fontsize=16)
#            f[j].set_ylabel(r"$\leftarrow$ n'", fontsize=16)
            f[j].set_xlabel(r"n", fontsize=16)
            f[j].set_ylabel(r"n'", fontsize=16)

            for side in ['top', 'bottom', 'left', 'right']:
                f[j].spines[side].set_visible(False)

#            f[j].axes.xaxis.set_ticks([])
#            f[j].axes.yaxis.set_ticks([])

    def set_limits(self, f, lims, origin='upper'):

        for j in range(2):

            f[j].set_xlim(lims)
            if origin == 'upper':
                f[j].set_ylim(np.flip(lims))
            elif origin == 'lower':
                f[j].set_ylim(lims)

##############################################################################

# Main function
def imshow_gkk2(

        gkk2_flist=None,
        qpt_weight=None,
        kpt_index=1, #Given in Fortran indexing
        valence=None,
        band_range=None,

        qpt_index=None,
        mode=None,

        # Plotting options
        single = True, # Plot only one gkk^2 matrix
        diff = False,  # Plot the difference and relative difference between 2 gkk2_fname files

        title=None,
        verbose=False,
        figsize=None,
        savefig_root=None,
        diff_labels=None,
        
        *args, **kwargs):


    gkk2_ = gkk2_flist

    qpt_weight_fname = qpt_weight
    kpt_index = kpt_index
    valence = valence 

    if not band_range:
        raise Exception('Must provide the band_range for imshow_gkk2 (in Fortran indexing)')
    band_range = band_range

    qpt_index = qpt_index
    if not qpt_index: # Eventually, modify to find the kpt and qpt coordinates in the data
        raise Exception("Must provide one qpt_index (in Fortran indexing)")

    mode = mode
    if not mode:
        raise Exception("Must provide phonon mode index (in Fortran indexing)")

    figsize = figsize
    title = title
    diff_labels = diff_labels
    
    gkk2 = Gkk2Plotter(

        gkk2_flist=gkk2_flist,

        kpt_index=kpt_index,
        qpt_index=qpt_index,
        band_range=band_range,
        mode=mode,
        valence=valence,

        figsize=figsize,
        title=title,
        verbose=verbose,
        savefig_root=savefig_root,
        diff_labels=diff_labels,

        single=single,
        diff=diff,

        **kwargs)

    return

