#! usr/bin/env python

import numpy as np
from zprfile import ZPRfile
import matplotlib.pyplot as plt
from constants import ha_to_ev
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import os
from colorpalettes import bright, vibrant, highcontrast 
import seaborn as sns
from scipy.interpolate import UnivariateSpline
rc('text', usetex=True)
rc('font', family='sans-serif')

class QptContribution(object):

    def __init__(self,

                 qptcont_flist=None,
                 qptmode_flist=None,
                 qptcont_farr=None,
                 qptmode_farr=None,

                 qpt_weight_fname=None,
                 mode=False,
                 cumulative=True,
                 ratio=False,
                 band_index=None,
                 bz_weight=True,
                 kpt_index=0,

                 figsize=(8, 8),
                 title=None,
                 verbose=False,
                 color = [highcontrast['blue'], highcontrast['yellow'], bright['green'], highcontrast['red'],
                          bright['cyan'], bright['gray'], bright['purple'], vibrant['orange'], vibrant['teal']],
                 linestyle = ['solid', 'dashed', 'dotted', 'dashdot'],
#                 linestyle = ['solid', 'dashed', (0, (7, 5)), 'dashdot', (0, (3, 1, 1, 1, 1, 1, 1, 1)), 
#                              (0, (1,1)), (0, (3,5,1,5)), (0, (2,2)), (0, (1,10))],
                 savefig=None,
                 flist_labels=None,

                 **kwargs):

        self.mode = mode
        self.cumulative = cumulative
        self.ratio = ratio
        self.qpt_weight_fname = qpt_weight_fname

        self.figsize = figsize
        self.title = title
        self.verbose = verbose
        self.bz_weight = bz_weight

        self.color = color
        self.linestyle = linestyle
        self.savefig = savefig
        self.flist_labels = flist_labels

        if qptcont_flist and qptcont_farr:
            raise Exception('Cannot simultaneously define qptcont_flist (for individual plots) and qptcont_farr (for ratios).')
        if qptcont_flist:
            self.nfile = len(qptcont_flist)
            self.band_index = band_index

            if len(self.band_index) != self.nfile:
                raise Exception('Band_index must have {} entries, but has {}.'.format(self.nfile, len(self.band_index)))

        elif qptcont_farr:
            self.nratio, self.nfile = np.shape(qptcont_farr)
            self.band_index = band_index

            if np.shape(self.band_index) != (self.nratio, self.nfile):
                raise Exception('Band_index must have shape ({}, {}), but has shape {}.'.format(self.nratio, self.nfile, np.shape(self.band_index)))

            soc_order = np.sign(self.band_index[0][0]-self.band_index[0][1])
            if soc_order<0:
                raise Exception('First fname in each flist of qptcont_farr should include SOC and have more bands that second fname')
            for n in range(self.nratio-1):
                if np.sign(self.band_index[n+1][0]-self.band_index[n+1][1]) != soc_order:
                    raise Exception('Check fname ordering for flist {}; it differs from flist 0'.format(n+1))


        if not qpt_weight_fname:
            raise Exception('Must provide a file for qpt weights.')
        else:
            self.set_wtq(qpt_weight_fname)

        if self.mode and len(qptmode_flist) != self.nfile:
            raise Exception('''For mode decomposition, qptcont_flist and qptmode_flist should have
                            have the same lenght, but have {} and {}.'''.format(self.nfile, len(self.qptmode_flist)))

        if self.mode and len(self.color) < 4:
            raise Exception('For Mode decomposition, custom color list must contain 4 entries.')

        if len(self.linestyle) < self.nfile:
            raise Exception('For {} files in flist, linestyle need {} entries, but has only {}.'.format(
                            self.nfile, self.nfile, len(self.linestyle)))


        if self.ratio and self.nfile != 2:
            raise Exception('Ratio between histograms requires 2 files in flist, but got {}'.format(self.nfile))

        if self.ratio:

            if self.mode:
                raise Exception('Mode not yet implemented for histogram ratio')

            else:
                for n in range(self.nratio):
                    print('\nComputing the ratio of the histograms of qpt contribution of\n{}\n    vs\n{}'.format(qptcont_farr[n][0], qptcont_farr[n][1]))
                    qptcont1 = ZPRfile(fname=qptcont_farr[n][0], read=False)
                    qptcont1.read_nc()

                    qptcont2 = ZPRfile(fname=qptcont_farr[n][1], read=False)
                    qptcont2.read_nc()

                    # get data
                    self.zpr_qpt = qptcont1.reduced_zpr_qpt[0, kpt_index, self.band_index[n][0], :]*ha_to_ev*1E3
                    self.qpoints = qptcont1.qpoints
                    self.nqpt = qptcont1.nqpt
                    if self.nqpt is None:
                        raise Exception('nqpt is None, check your input file')
                    self.zpr_qpt2 = qptcont2.reduced_zpr_qpt[0, kpt_index, self.band_index[n][1], :]*ha_to_ev*1E3

                    if (self.qpoints != qptcont2.qpoints).any():
                        raise Exception('Q-point lists differ between both files.')

                    # This is quite unnecessary...
                    if self.nqpt != qptcont2.nqpt:
                        raise Exception('Q-point number differ between both files.')

                    self.get_qpt_norm()

                    if self.cumulative:
                        raise  Exception('Cumulative not implemented for ratio')
                        #fig, ax = plt.subplots(3, 1, squeeze=True, sharex=True, figsize=self.figsize)
                    else:
                        if n == 0:
                            fig = plt.figure(figsize=self.figsize)
                            spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig)
                            if self.bz_weight:
                                ax = [fig.add_subplot(spec[0:2, 0]),  fig.add_subplot(spec[2, 0])]
                            else:
                                # This one is to have only the histograms, no BZ weight
                                ax = [fig.add_subplot(spec[0:2, 0])]

                    self.plot_qpt_contribution_ratio(n, ax)

        else:
            for n in range(self.nfile):
                print('\nTreating {}...'.format(os.path.basename(qptcont_flist[n])))
                qptcont = ZPRfile(fname=qptcont_flist[n], read=False)
                qptcont.read_nc()

                if self.mode:
                    qptmode = ZPRfile(qptmode_flist[n], read=False)
                    qptmode.read_nc()
                    self.nmode = qptmode.nmodes
                    self.natom = np.int(self.nmode/3)

                # get data
                self.zpr_qpt = qptcont.reduced_zpr_qpt[0, kpt_index, self.band_index[n], :]*ha_to_ev*1E3
                self.qpoints = qptcont.qpoints
                self.nqpt = qptcont.nqpt
                if self.nqpt is None:
                    raise Exception('nqpt is None, check your input file')

                self.get_qpt_norm()

                if self.mode:
                    self.zpr_mode_qpt = qptmode.reduced_zpr_qpt_mode[0, kpt_index, self.band_index[n], :, :]*ha_to_ev*1E3
                    self.split_modes()

                if self.mode:
                    if n == 0:
                        if self.cumulative:
                            # Add the if not self.bz_weight
                            fig, ax = plt.subplots(3, 1, squeeze=True, sharex=True, figsize=self.figsize)
                        else:
                            fig = plt.figure(figsize=self.figsize)
                            spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig)
                            if not self.bz_weight:
                                ax = [fig.add_subplot(spec[0:2, 0])]
                            else:
                                ax = [fig.add_subplot(spec[0:2, 0]),  fig.add_subplot(spec[2, 0])]
    #                        fig, ax = plt.subplots(2, 1, squeeze=True, sharex=True, figsize=self.figsize)
                        self.ref_zpr = np.sum(self.zpr_qpt)
                    self.plot_qpt_mode_contribution(n, ax)
                else:
                    if n == 0:
                        if self.cumulative:
                            fig, ax = plt.subplots(3, 1, squeeze=True, sharex=True, figsize=self.figsize)
                        else:
                            fig = plt.figure(figsize=self.figsize)
                            spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig)
                            if self.bz_weight:
                                ax = [fig.add_subplot(spec[0:2, 0]),  fig.add_subplot(spec[2, 0])]
                            else:
                                # This one is to have only the histograms, no BZ weight
                                ax = [fig.add_subplot(spec[0:2, 0])]

                    self.plot_qpt_contribution(n, ax)

    def set_wtq(self, fname):
        # Read qpt weights from file
        wtq = []
        with open(fname, 'r') as f:
            lines = f.readlines()
            nline = len(lines)
            for i, line in enumerate(lines):
                if i == nline-1:
                    w = (line.split('\n')[0]).split(',')
                else:
                    w = (line.split('\n')[0]).split(',')[:-1]
                for iw in w:
                    wtq.append(np.float(iw.split('[')[-1]))

        self.wtq = np.array(wtq)/np.sum(wtq)

    def get_qpt_norm(self):

        self.qpt_norm = np.zeros((self.nqpt))

        for i, qpt in enumerate(self.qpoints):

            self.qpt_norm[i] = np.linalg.norm(qpt)*np.sqrt(2)

    def split_modes(self):

        self.acoustic = self.zpr_mode_qpt[:, :3]
        self.to = self.zpr_mode_qpt[:, 3:(2*(self.natom-1)+3)]
        self.lo = self.zpr_mode_qpt[:,-(self.natom-1)]

        self.acoustic = np.einsum('qv->q', self.acoustic)
        self.to = np.einsum('qv->q', self.to)

        if self.natom >2:
            self.lo = np.einsum('qv->q', self.lo)

        self.check_mode_sum()

    def check_mode_sum(self):

        total_zpr = np.sum(self.zpr_qpt)
        total_modes = np.sum(self.acoustic)+np.sum(self.to)+np.sum(self.lo)

        print('Checking mode decomposition...')
        if np.abs(total_zpr - total_modes) > 1E-4:
            print('Check your input files/calculation : mode decomposition does not sum up to ZPR.')
            print('total zpr: {}'.format(np.sum(self.zpr_qpt)))
            print('total modes: {}'.format(np.sum(self.acoustic)+np.sum(self.to)+np.sum(self.lo)))

        if self.verbose:
            print('  total zpr: {:>8.4f} meV'.format(np.sum(self.zpr_qpt)))
            print('  acoustic: {:>8.4f} meV'.format(np.sum(self.acoustic)))
            print('  TO: {:>8.4f} meV'.format(np.sum(self.to)))
            print('  LO: {:>8.4f} meV'.format(np.sum(self.lo)))
        print('... all good!')

    def plot_qpt_mode_contribution(self, n, ax):

        alpha = [0.6, 1.0]

        if self.cumulative:
            # Add if not bz_weight:
    #        ax.plot(self.qpt_norm, self.zpr_qpt, 'o')
            # ZPR contribution
            ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='step', linewidth=1.5, stacked=True,
                       color=self.color[0], linestyle = self.linestyle[n])
            ax[0].hist(self.qpt_norm, bins=50, weights=self.acoustic, histtype='step', linewidth=1.5, stacked=True,
                       color=self.color[1], linestyle = self.linestyle[n])
            ax[0].hist(self.qpt_norm, bins=50, weights=self.to, histtype='step', linewidth=1.5, stacked=True,
                       color=self.color[2], linestyle = self.linestyle[n])
            ax[0].hist(self.qpt_norm, bins=50, weights=self.lo, histtype='step', linewidth=1.5, stacked=True, 
                       color=self.color[3], linestyle = self.linestyle[n])
            histzpr, binszpr = self.get_histogram(self.zpr_qpt)
            ax[0].plot(binszpr[:-1], histzpr, 'r')
            spl_hist = UnivariateSpline(bins[15:-1], histzpr, k=3)
            xq = np.linspace(binszpr[15], binszpr[-1], 100)
            ax[0].plot(xq, spl_hist(xq), 'y', lw=3)
            # CUmulative ZPR contribution
    #        ax[1].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, cumulative=True, density=True,
    #                   histtype='step', linewidth=1.5, color=self.color[0], linestyle = self.linestyle[n])
            hist, bins = self.get_cumulative_histogram(self.zpr_qpt)
            ax[1].hist(bins[:-1], bins=bins, weights=hist, histtype='step', linewidth=1.5, color=self.color[0], linestyle = self.linestyle[n])
            hist, bins = self.get_cumulative_histogram(self.acoustic)
            ax[1].hist(bins[:-1], bins=bins, weights=hist, histtype='step', linewidth=1.5, color=self.color[1], linestyle = self.linestyle[n])
            hist, bins = self.get_cumulative_histogram(self.to)
            ax[1].hist(bins[:-1], bins=bins, weights=hist, histtype='step', linewidth=1.5, color=self.color[2], linestyle = self.linestyle[n])
            hist, bins = self.get_cumulative_histogram(self.lo)
            ax[1].hist(bins[:-1], bins=bins, weights=hist, histtype='step', linewidth=1.5, color=self.color[3], linestyle = self.linestyle[n])


            # Qpoint density
            hist, bins = self.get_bz_histogram()  # FIX ME : this does not work very well.... should be removed

    #        sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50, ax=ax[2], color='gray', element='step', fill=False, 
    #                     kde=True, kde_kws={'cut':0}, line_kws={'color':'r'})
    #        sns.histplot(x=bins[:-1], weights=hist, bins=50, ax=ax[2], color='gray', element='step', fill=False, 
    #                     kde=True, kde_kws={'cut':0}, line_kws={'color':'r'})

            histbz, bins = self.get_histogram(self.wtq)
            line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50,kde=True).get_lines()[0].get_data()
            ax[2].clear()
            ax[2].plot(line[0], line[1], 'r')

            # extend histogram data to the same lenght as the kde
            longhist = np.zeros((200))
            longhistzpr = np.zeros((200))
            for i in range(50):
                longhistzpr[4*i:4*(i+1)] = histzpr[i]
                longhist[4*i:4*(i+1)] = histbz[i]
            longcorr = line[1]/longhist
            corr = line[1][::4]/histbz
            ax[0].plot(bins[:-1], corr*histzpr, 'g')
            ax[0].plot(line[0], longcorr*longhistzpr, 'm')

            sns.kdeplot(x=self.qpt_norm, weights=self.wtq, ax=ax[2], color=self.color[0], linestyle='solid', fill=False, cut=0)
    #        ax[2].hist(self.qpt_norm, bins=50, weights=self.wtq, histtype='step', linewidth=5.0, color=self.color[0], linestyle = self.linestyle[n])
            #self.shade_lowq(ax, 0.2)
    #        print(len(sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50,kde=True).get_lines()))

    #        print(ax[1].get_ylim())
            xlim = ax[0].get_xlim()
            ax[0].plot(np.linspace(xlim[0], xlim[1], 10), np.zeros((10)), 'k', zorder=-1, linewidth=0.75)
            for i in range(3):
                ax[i].set_xlim(xlim[0], xlim[1])

    #        ax[0].set_ylim(-0.1, 1.505)
            ax[1].set_ylim(0.0, 1.05)

        else:
            # No culumative histogram
    #        ax.plot(self.qpt_norm, self.zpr_qpt, 'o')
            # ZPR contribution
            ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='stepfilled', linewidth=1.5, stacked=True,
                       color=self.color[0], linestyle = self.linestyle[n],alpha=alpha[n], zorder=n)
            ax[0].hist(self.qpt_norm, bins=50, weights=self.lo, histtype='stepfilled', linewidth=1.5, stacked=True, 
                       color=self.color[3], linestyle = self.linestyle[n],alpha=alpha[n], zorder=2+n)
            ax[0].hist(self.qpt_norm, bins=50, weights=self.acoustic, histtype='stepfilled', linewidth=1.5, stacked=True,
                       color=self.color[1], linestyle = self.linestyle[n],alpha=alpha[n], zorder=6+n)
            ax[0].hist(self.qpt_norm, bins=50, weights=self.to, histtype='stepfilled', linewidth=1.5, stacked=True,
                       color=self.color[2], linestyle = self.linestyle[n],alpha=alpha[n], zorder=4+n)
            # This is just to outline the LO for Ge...
            #ax[0].hist(self.qpt_norm, bins=50, weights=self.lo, histtype='step', linewidth=1.5, stacked=True, 
            #           color=self.color[3], linestyle = self.linestyle[n],alpha=alpha[n], zorder=8+n)

            histzpr, binszpr = self.get_histogram(self.zpr_qpt)
      #      ax[0].plot(binszpr[:-1], histzpr, 'r')
#            spl_hist = UnivariateSpline(binszpr[15:-1], histzpr[15:], k=4)
#            xq = np.linspace(binszpr[15], binszpr[-1], 100)
#            ax[0].plot(xq, spl_hist(xq), 'y', lw=3)
            # CUmulative ZPR contribution

            if self.bz_weight:
                # Qpoint density
                hist, bins = self.get_bz_histogram()  # FIX ME : this does not work very well.... should be removed

                line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50, ax=ax[1], color='gray', element='step', fill=True, 
                              kde=True, kde_kws={'cut':0}).get_lines()[0].get_data()
        #        sns.histplot(x=bins[:-1], weights=hist, bins=50, ax=ax[2], color='gray', element='step', fill=False, 
        #                     kde=True, kde_kws={'cut':0}, line_kws={'color':'r'})

                histbz, bins = self.get_histogram(self.wtq)
                if self.verbose:
                    print('histbz sum (should be very close to 1):', np.sum(histbz))
    #            line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50,kde=True).get_lines()[0].get_data()
          #      ax[1].clear()
                ax[1].plot(line[0], line[1], 'k', linewidth=1.5)
        #        ax[1].plot(bins[15], histbz[15],'ok')
                # extend histogram data to the same lenght as the kde
    #            longhist = np.zeros((200))
    #            longhistzpr = np.zeros((200))
    #            for i in range(50):
    #                longhistzpr[4*i:4*(i+1)] = histzpr[i]
    #                longhist[4*i:4*(i+1)] = histbz[i]
    #            longcorr = line[1]/longhist
    #            corr = line[1][::4]/histbz
        #        ax[0].plot(binszpr[:-1], corr*histzpr, 'g')
    #            ax[0].plot(line[0], longcorr*longhistzpr, 'm')

    #            sns.kdeplot(x=self.qpt_norm, weights=self.wtq, ax=ax[1], color=self.color[0], linestyle='solid', fill=False, cut=0)
        #        ax[2].hist(self.qpt_norm, bins=50, weights=self.wtq, histtype='step', linewidth=5.0, color=self.color[0], linestyle = self.linestyle[n])
                #self.shade_lowq(ax, 0.2)
        #        print(len(sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50,kde=True).get_lines()))

        #        print(ax[1].get_ylim())
                #xlim = ax[0].get_xlim()
            xlim = [0,1]
            ax[0].plot(np.linspace(xlim[0], xlim[1], 10), np.zeros((10)), 'k', zorder=-1, linewidth=0.75)
            if self.bz_weight:
                for i in range(2):
                    ax[i].set_xlim(xlim[0], xlim[1])
                    ax[i].tick_params(axis='both', labelsize=12)
                    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f')) 
                ax[0].xaxis.set_visible(False)

            else:
                ax[0].set_xlim(xlim[0], xlim[1])
                ax[0].tick_params(axis='both', labelsize=12)
                ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f')) 

#            ax[0].set_ylim(-0.1, 1.505)

        if n+1  == self.nfile:

            # For CdTe
#            ax[0].annotate('',xy=(0.15,0.70), xycoords='data', xytext=(0.20, 0.90),arrowprops=dict(color='k',headlength=12,headwidth=8,width=0.2))
#            ax[0].annotate('',xy=(0.70,0.30), xycoords='data', xytext=(0.70, 0.50),arrowprops=dict(color='k',headlength=12,headwidth=8,width=0.2))
#            ax[0].text(0.10, 0.92, r'change of $m^\ast$', fontsize=20)
#            ax[0].text(0.52, 0.51, r'global $\varepsilon_{\mathbf{k}n}$ reduction', fontsize=20)
            # For CdS
            ax[0].annotate('',xy=(0.15,2.20), xycoords='data', xytext=(0.20, 2.50),arrowprops=dict(color='k',headlength=12,headwidth=8,width=0.2))
            ax[0].annotate('',xy=(0.70,0.63), xycoords='data', xytext=(0.70, 1.20),arrowprops=dict(color='k',headlength=12,headwidth=8,width=0.2))
            ax[0].text(0.10, 2.58, r'change of $m^\ast$', fontsize=20)
            ax[0].text(0.52, 1.23, r'global $\varepsilon_{\mathbf{k}n}$ reduction', fontsize=20)

            self.set_hist_labels(ax)

            self.set_legend(ax[0])
            
            if self.title:
                plt.suptitle(self.title, fontsize=20)

            if self.bz_weight:
                plt.subplots_adjust(hspace=0.07, top=0.93)
            else:
                plt.subplots(top=0.3)


            if self.savefig:
                plt.savefig('{}_qpt_mode.pdf'.format(self.savefig), dpi=1200, bbox_inches='tight')
                os.system('open {}_qpt_mode.pdf'.format(self.savefig))
            else:
                plt.show()

    def plot_qpt_contribution(self, n, ax):

#        ax.plot(self.qpt_norm, self.zpr_qpt, 'o')

        if self.cumulative:
            ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='step', linestyle = self.linestyle[n], color=self.color[n], lw=2.0)
            ax[1].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, cumulative=True, density=True, histtype='step', linestyle = self.linestyle[n], color=self.color[n], lw=2.0)
            ax[2].hist(self.qpt_norm, bins=50, weights=self.wtq, histtype='step', linestyle = self.linestyle[n], color=self.color[n], lw=2.0)

            if n+1 == self.nfile:
                self.set_hist_labels(ax)
                self.set_legend(ax[0])

                if self.title:
                    plt.suptitle(self.title, fontsize=14)
                if self.savefig:
                    plt.savefig('{}_qpt.pdf'.format(self.savefig), dpi=1200)
                    os.system('open {}_qpt.pdf'.format(self.savefig))
                else:
                    plt.show()
        
        else:
            ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='step', linestyle = self.linestyle[n], color=self.color[n], lw=3.0)

            if self.bz_weight:
                # Qpoint density
                hist, bins = self.get_bz_histogram()  # FIX ME : this does not work very well.... should be removed

                line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50, ax=ax[1], color='gray', element='step', fill=True, 
                              kde=True, kde_kws={'cut':0}).get_lines()[0].get_data()

                histbz, bins = self.get_histogram(self.wtq)
                if self.verbose:
                    print('histbz sum (should be very close to 1):', np.sum(histbz))
                ax[1].plot(line[0], line[1], 'k', linewidth=1.5)

            if n+1 == self.nfile:
                self.set_hist_labels(ax)
                self.set_legend(ax[0])

                xlim = [0,1]
                if self.bz_weight:
                    for i in range(2):
                        ax[i].set_xlim(xlim[0], xlim[1])
                        ax[i].tick_params(axis='both', labelsize=14)
                    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f')) 
                    ax[0].xaxis.set_visible(False)

                else:
                    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                    ax[0].plot(np.linspace(xlim[0], xlim[1], 10), np.zeros((10)), 'k', zorder=-1, linewidth=0.75)
                    ax[0].set_xlim(xlim[0], xlim[1])
                    ax[0].tick_params(axis='both', labelsize=14)

                if self.title:
                    plt.suptitle(self.title, fontsize=20)

                if self.bz_weight:
                    plt.subplots_adjust(hspace=0.07, top=0.93)
                else:
                    plt.subplots_adjust(top=0.93)
                if self.savefig:
    #                plt.savefig('{}_qpt.pdf'.format(self.savefig), dpi=1200)
                    plt.savefig('{}_qpt.pdf'.format(self.savefig), dpi=1200, bbox_inches='tight')

                    os.system('open {}_qpt.pdf'.format(self.savefig))
                else:
                    plt.show()

    def plot_qpt_contribution_ratio(self, n, ax):

        hist, bins = self.get_histogram(self.zpr_qpt)
        hist2, bins2 = self.get_histogram(self.zpr_qpt2)

        if (bins != bins2).any():
            raise Exception('Bins are different')

        ratio = hist/hist2
#        print(ratio)
#        ratio = self.zpr_qpt/self.zpr_qpt2
#        ax[0].stairs(ratio[1:-1], bins[1:-1], color=self.color[0], fill=False)
        ax[0].stairs(ratio[:], bins[:], color=self.color[n], fill=False, lw=3.0)


         # Since the drop is almost constant in each bin, I guess that taking the ratio before taking the histogram distribution comes up to 
         # summing BZ weight
#        ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt/self.zpr_qpt2, histtype='step', color=self.color[1], lw=3.0)
        if self.bz_weight:
            # Qpoint density
            hist, bins = self.get_bz_histogram()  # FIX ME : this does not work very well.... should be removed

            line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50, ax=ax[1], color='gray', element='step', fill=True, 
                          kde=True, kde_kws={'cut':0}).get_lines()[0].get_data()

            histbz, bins = self.get_histogram(self.wtq)
            if self.verbose:
                print('histbz sum (should be very close to 1):', np.sum(histbz))
            ax[1].plot(line[0], line[1], 'k', linewidth=1.5)

        if n+1 == self.nfile:
            self.set_hist_labels(ax)
            self.set_legend(ax[0])

            xlim = [0,1]
            if self.bz_weight:
                for i in range(2):
                    ax[i].set_xlim(xlim[0], xlim[1])
                    ax[i].tick_params(axis='both', labelsize=14)
                ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f')) 
                ax[0].xaxis.set_visible(False)

            else:
                ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f')) 
                ax[0].plot(np.linspace(xlim[0], xlim[1], 10), np.zeros((10)), 'k', zorder=-1, linewidth=0.75)
                ax[0].plot(np.linspace(xlim[0], xlim[1], 10), np.ones((10)), 'k', linestyle='dashed', zorder=-1, linewidth=0.75)
                ax[0].set_xlim(xlim[0], xlim[1])
                ax[0].tick_params(axis='both', labelsize=14)

            if self.title:
                plt.suptitle(self.title, fontsize=20)

            if self.bz_weight:
                plt.subplots_adjust(hspace=0.07, top=0.93)
            else:
                plt.subplots_adjust(top=0.93)
            if self.savefig:
    #                plt.savefig('{}_qpt.pdf'.format(self.savefig), dpi=1200)
                plt.savefig('{}_ratio_qpt.pdf'.format(self.savefig), dpi=1200, bbox_inches='tight')

                os.system('open {}_ratio_qpt.pdf'.format(self.savefig))
            else:
                plt.show()


    def set_hist_labels(self, ax):

        fs = 18
        if self.cumulative:
            # Add if not self.bz_weight
            ax[2].set_xlabel(r'Q-point norm ($\frac{2\pi}{a}$)', fontsize=fs)
            if self.ratio:
                ax[0].set_ylabel(r'Ratio of contribution to ZPR SOC/noSOC', fontsize=fs)
            else:
                ax[0].set_ylabel(r'Contribution to ZPR (meV)', fontsize=fs)

            if self.flist_labels and not self.ratio:
                ax[1].set_ylabel(r'Cumulative contribution\\to ZPR / Total ZPR ({})'.format(self.flist_labels[0]), fontsize=12)
            else:
                ax[1].set_ylabel(r'Cumulative contribution\\to ZPR / Total ZPR', fontsize=fs)

    #        ax[2].set_ylabel(r'Number of q-points', fontsize=fs)
            ax[2].set_ylabel(r'BZ weight', fontsize=fs)

    #        ax[0].set_title(r'Contribution to ZPR vs $\vert q\vert$', fontsize=12)
    #        ax[1].set_title(r'Cumulative contribution vs $\vert q\vert$', fontsize=12)
    #        ax[2].set_title(r'Q-point distribution', fontsize=12)
        else:
            if not self.bz_weight:
                ax[0].set_xlabel(r'$\mathbf{q}$-point norm $q$ ($\frac{2\pi}{a}$)', fontsize=fs)
            else:
                ax[1].set_xlabel(r'$\mathbf{q}$-point norm $q$ ($\frac{2\pi}{a}$)', fontsize=fs)
                ax[1].set_ylabel(r'BZ weight', fontsize=fs)
            if self.ratio:
                ax[0].set_ylabel(r'Ratio of contribution to ZPR SOC/noSOC', fontsize=fs)
            else:
                ax[0].set_ylabel(r'Contribution to ZPR (meV)', fontsize=fs)

    #        ax[2].set_ylabel(r'Number of q-points', fontsize=12)

    #        ax[0].set_title(r'Contribution to ZPR vs $\vert q\vert$', fontsize=12)
    #        ax[1].set_title(r'Cumulative contribution vs $\vert q\vert$', fontsize=12)
    #        ax[2].set_title(r'Q-point distribution', fontsize=12)




    def set_legend(self, ax):

        handles = []
        if self.mode:
            labels = ['Total', 'Acoustic', 'TO', 'LO']
            for i in range(4):
                handles.append(Line2D([0], [0], linewidth=4.0, color=self.color[i], label=labels[i]))

            if self.cumulative:
                artist = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.62, 0.75), ncol=4, fontsize=12)
                ax.add_artist(artist)
            else:
    #            artist = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.50, 1.01), ncol=4, fontsize=16)
                artist = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.60, 0.85), ncol=4, fontsize=16,
                                   columnspacing=1.1, handlelength=1.5)

                ax.add_artist(artist)

        if self.ratio:
            if self.nfile>1 and self.flist_labels:
                handles = []
                for i in range(self.nfile):
                    handles.append(Line2D([0], [0], linewidth=3.0, color=self.color[i], label=self.flist_labels[i], linestyle=self.linestyle[i]))
#                artist2 = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.885, 0.665), ncol=1, fontsize=16, handlelength=1.0,
#                                    handletextpad=0.5)
                artist2 = ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=self.nfile, fontsize=16, handlelength=1.5,
                                    handletextpad=0.5)

                ax.add_artist(artist2)

        else:
            if self.nfile > 1 and self.flist_labels:
                handles = []
                alpha = [0.5, 1.0]
                for i in range(self.nfile):
                    if self.mode:
                        #handles.append(Line2D([0], [0], linewidth=3.0, color='k', label=self.flist_labels[i], linestyle='solid', alpha=alpha[i]))
                        handles.append(Patch(edgecolor='k', facecolor='k', alpha=alpha[i], label=self.flist_labels[i], linewidth=4.0))
                        #[Patch(edgecolor='k', facecolor='k', hatch='///', label=r'Ce-$f$')]
                    else:
                        handles.append(Line2D([0], [0], linewidth=3.0, color=self.color[i], label=self.flist_labels[i], linestyle=self.linestyle[i]))


                if self.cumulative:
                    artist2 = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.81, 0.57), ncol=4, fontsize=12)
                else:
                    # Main paper, no modes, 5 histograms
                    #artist2 = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.835, 0.53), ncol=1, handlelength=2.5, fontsize=16)
                    # Main paper, CdTe and CdS
                    #artist2 = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.885, 0.665), ncol=1, fontsize=16, handlelength=1.0,
                    #                    handletextpad=0.5)
                    # SM, no modes, 9 histograms
                    artist2 = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.635, 0.565),
                                        ncol=2, handlelength=2.5, fontsize=14, columnspacing=1.5)



    def get_cumulative_histogram(self, y):

        hist, bins = np.histogram(self.qpt_norm, bins=50, weights=y)
        cumul_hist = np.zeros_like(hist)

        for i in range(len(hist)):
            cumul_hist[i] = np.sum(hist[:i])

        cumul_hist = cumul_hist/cumul_hist[-1]*np.sum(y)/np.sum(self.ref_zpr)

        return cumul_hist, bins

    def get_bz_histogram(self):
        
        hist, bins = np.histogram(self.qpt_norm, bins=50, weights=self.wtq)
        trpz = np.trapz(y=hist, x=bins[:-1], axis=0)

        return hist/trpz, bins

    def get_histogram(self, y):

        hist, bins = np.histogram(self.qpt_norm, bins=50, weights=y)

        return hist, bins


class QptPathContribution(object):

    def __init__(self,

                 qptpath_flist=None,
                 qpt_weight_fname=None,
                 band_index=None,
                 kpt_index=0,
                 cumulative=False,
                 with_gkk=True,

                 figsize=(8, 8),
                 title=None,
                 verbose=False,
                 color = [highcontrast['blue'], highcontrast['yellow'], bright['green'], highcontrast['red'], bright['cyan'], bright['gray']],
                 linestyle = ['solid', 'dashed', 'dotted', 'dashdot', (0, (3, 1, 1, 1, 1, 1))],
                 savefig=None,
                 flist_labels=None,

                 **kwargs):

        self.band_index = band_index
        self.qpt_weight_fname = qpt_weight_fname

        self.figsize = figsize
        self.title = title
        self.verbose = verbose

        self.color = color
        self.linestyle = linestyle
        self.savefig = savefig
        self.flist_labels = flist_labels
        self.cumulative = cumulative
        self.with_gkk = with_gkk

        self.nfile = len(qptpath_flist)

        if not qpt_weight_fname:
            raise Exception('Must provide a file for qpt weights.')
        else:
            self.set_wtq(qpt_weight_fname)

        #if self.mode and len(qptmode_flist) != self.nfile:
        #    raise Exception('''For mode decomposition, qptcont_flist and qptmode_flist should have
        #                    have the same lenght, but have {} and {}.'''.format(nfile, len(self.qptmode_flist)))

        #if self.mode and len(self.color) < 4:
        #    raise Exception('For Mode decomposition, custom color list must contain 4 entries.')

        if len(self.linestyle) < self.nfile:
            raise Exception('For {} files in flist, linestyle need {} entries, but has only {}.'.format(
                            self.nfile, self.nfile, len(self.linestyle)))

        if len(self.band_index) != self.nfile:
            raise Exception('Band_index must have {} entries, but has {}.'.format(self.nfile, len(self.band_index)))

        for n in range(self.nfile):
            qptpath = ZPRfile(fname=qptpath_flist[n], read=False)
            qptpath.read_nc()

            self.nmode = qptpath.nmodes
            self.natom = np.int(self.nmode/3)

            # get data
            zpr = qptpath.zpr_qpath[0, kpt_index, self.band_index[n], :, :]*ha_to_ev*1E3
            # sum on modes, as the qpath retains the mode dependency
            self.zpr_qpt = np.einsum('qv->q', zpr)

            if self.with_gkk:
                '''FFIX ME : units?!?'''
                gkk2_fan = qptpath.gkk2_qpath_fan[0, kpt_index, self.band_index[n], :, :]*ha_to_ev*1E3
                gkk2_ddw = qptpath.gkk2_qpath_ddw[0, kpt_index, self.band_index[n], :, :]*ha_to_ev*1E3
                gkk2 = gkk2_fan - gkk2_ddw

                self.gkk2 = np.einsum('qv->q', gkk2)
            
            self.qpoints = qptpath.qpoints
            self.nqpt = qptpath.nqpt
            if self.nqpt is None:
                raise Exception('nqpt is None, check your input file')

            self.get_qpt_norm()

            #if self.mode:
            #    self.zpr_mode_qpt = qptmode.reduced_zpr_qpt_mode[0, kpt_index, self.band_index[n], :, :]*ha_to_ev*1E3
            #    self.split_modes()

#            if mode:
#                if n == 0:
#                    if self.cumulative:
#                        fig, ax = plt.subplots(3, 1, squeeze=True, sharex=True, figsize=self.figsize)
#                    else:
#                        fig = plt.figure(figsize=self.figsize)
#                        spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig)
#                        ax = [fig.add_subplot(spec[0:2, 0]),  fig.add_subplot(spec[2, 0])]
##                        fig, ax = plt.subplots(2, 1, squeeze=True, sharex=True, figsize=self.figsize)
#                    self.ref_zpr = np.sum(self.zpr_qpt)
#                self.plot_qpt_mode_contribution(n, ax)
#            else:
            if n == 0:
                if self.with_gkk:
                    fig, ax = plt.subplots(3, 1, squeeze=True, sharex=True, figsize=self.figsize)
                else:
                    fig, ax = plt.subplots(2, 1, squeeze=True, sharex=True, figsize=self.figsize)
            self.plot_qpt_contribution(n, ax)

    def set_wtq(self, fname):
        # Read qpt weights from file
        wtq = []
        with open(fname, 'r') as f:
            lines = f.readlines()
            nline = len(lines)
            for i, line in enumerate(lines):
                if i == nline-1:
                    w = (line.split('\n')[0]).split(',')
                else:
                    w = (line.split('\n')[0]).split(',')[:-1]
                for iw in w:
                    wtq.append(np.float(iw.split('[')[-1]))

        self.wtq = np.array(wtq)/np.sum(wtq)

    def get_qpt_norm(self):

        self.qpt_norm = np.zeros((self.nqpt))

        for i, qpt in enumerate(self.qpoints):

            self.qpt_norm[i] = np.linalg.norm(qpt)*np.sqrt(2)

    def split_modes(self):

        self.acoustic = self.zpr_mode_qpt[:, :3]
        self.to = self.zpr_mode_qpt[:, 3:(2*(self.natom-1)+3)]
        self.lo = self.zpr_mode_qpt[:,-(self.natom-1)]

        self.acoustic = np.einsum('qv->q', self.acoustic)
        self.to = np.einsum('qv->q', self.to)

        if self.natom >2:
            self.lo = np.einsum('qv->q', self.lo)

        self.check_mode_sum()

    def check_mode_sum(self):

        total_zpr = np.sum(self.zpr_qpt)
        total_modes = np.sum(self.acoustic)+np.sum(self.to)+np.sum(self.lo)

        print('Checking mode decomposition...')
        if np.abs(total_zpr - total_modes) > 1E-4:
            print('Check your input files/calculation : mode decomposition does not sum up to ZPR.')
            print('total zpr: {}'.format(np.sum(self.zpr_qpt)))
            print('total modes: {}'.format(np.sum(self.acoustic)+np.sum(self.to)+np.sum(self.lo)))

        if self.verbose:
            print('  total zpr: {:>8.4f} meV'.format(np.sum(self.zpr_qpt)))
            print('  acoustic: {:>8.4f} meV'.format(np.sum(self.acoustic)))
            print('  TO: {:>8.4f} meV'.format(np.sum(self.to)))
            print('  LO: {:>8.4f} meV'.format(np.sum(self.lo)))
        print('... all good!')

    def plot_qpt_mode_contribution(self, n, ax):

        alpha = [0.6, 1.0]

        if self.cumulative:
    #        ax.plot(self.qpt_norm, self.zpr_qpt, 'o')
            # ZPR contribution
            ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='step', linewidth=1.5, stacked=True,
                       color=self.color[0], linestyle = self.linestyle[n])
            ax[0].hist(self.qpt_norm, bins=50, weights=self.acoustic, histtype='step', linewidth=1.5, stacked=True,
                       color=self.color[1], linestyle = self.linestyle[n])
            ax[0].hist(self.qpt_norm, bins=50, weights=self.to, histtype='step', linewidth=1.5, stacked=True,
                       color=self.color[2], linestyle = self.linestyle[n])
            ax[0].hist(self.qpt_norm, bins=50, weights=self.lo, histtype='step', linewidth=1.5, stacked=True, 
                       color=self.color[3], linestyle = self.linestyle[n])
            histzpr, binszpr = self.get_histogram(self.zpr_qpt)
            ax[0].plot(binszpr[:-1], histzpr, 'r')
            spl_hist = UnivariateSpline(bins[15:-1], histzpr, k=3)
            xq = np.linspace(binszpr[15], binszpr[-1], 100)
            ax[0].plot(xq, spl_hist(xq), 'y', lw=3)
            # CUmulative ZPR contribution
    #        ax[1].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, cumulative=True, density=True,
    #                   histtype='step', linewidth=1.5, color=self.color[0], linestyle = self.linestyle[n])
            hist, bins = self.get_cumulative_histogram(self.zpr_qpt)
            ax[1].hist(bins[:-1], bins=bins, weights=hist, histtype='step', linewidth=1.5, color=self.color[0], linestyle = self.linestyle[n])
            hist, bins = self.get_cumulative_histogram(self.acoustic)
            ax[1].hist(bins[:-1], bins=bins, weights=hist, histtype='step', linewidth=1.5, color=self.color[1], linestyle = self.linestyle[n])
            hist, bins = self.get_cumulative_histogram(self.to)
            ax[1].hist(bins[:-1], bins=bins, weights=hist, histtype='step', linewidth=1.5, color=self.color[2], linestyle = self.linestyle[n])
            hist, bins = self.get_cumulative_histogram(self.lo)
            ax[1].hist(bins[:-1], bins=bins, weights=hist, histtype='step', linewidth=1.5, color=self.color[3], linestyle = self.linestyle[n])


            # Qpoint density
            hist, bins = self.get_bz_histogram()  # FIX ME : this does not work very well.... should be removed

    #        sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50, ax=ax[2], color='gray', element='step', fill=False, 
    #                     kde=True, kde_kws={'cut':0}, line_kws={'color':'r'})
    #        sns.histplot(x=bins[:-1], weights=hist, bins=50, ax=ax[2], color='gray', element='step', fill=False, 
    #                     kde=True, kde_kws={'cut':0}, line_kws={'color':'r'})

            histbz, bins = self.get_histogram(self.wtq)
            line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50,kde=True).get_lines()[0].get_data()
            ax[2].clear()
            ax[2].plot(line[0], line[1], 'r')

            # extend histogram data to the same lenght as the kde
            longhist = np.zeros((200))
            longhistzpr = np.zeros((200))
            for i in range(50):
                longhistzpr[4*i:4*(i+1)] = histzpr[i]
                longhist[4*i:4*(i+1)] = histbz[i]
            longcorr = line[1]/longhist
            corr = line[1][::4]/histbz
            ax[0].plot(bins[:-1], corr*histzpr, 'g')
            ax[0].plot(line[0], longcorr*longhistzpr, 'm')

            sns.kdeplot(x=self.qpt_norm, weights=self.wtq, ax=ax[2], color=self.color[0], linestyle='solid', fill=False, cut=0)
    #        ax[2].hist(self.qpt_norm, bins=50, weights=self.wtq, histtype='step', linewidth=5.0, color=self.color[0], linestyle = self.linestyle[n])
            #self.shade_lowq(ax, 0.2)
    #        print(len(sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50,kde=True).get_lines()))

    #        print(ax[1].get_ylim())
            xlim = ax[0].get_xlim()
            ax[0].plot(np.linspace(xlim[0], xlim[1], 10), np.zeros((10)), 'k', zorder=-1, linewidth=0.75)
            for i in range(3):
                ax[i].set_xlim(xlim[0], xlim[1])

    #        ax[0].set_ylim(-0.1, 1.505)
            ax[1].set_ylim(0.0, 1.05)

        else:
    #        ax.plot(self.qpt_norm, self.zpr_qpt, 'o')
            # ZPR contribution
            ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='stepfilled', linewidth=1.5, stacked=True,
                       color=self.color[0], linestyle = self.linestyle[n],alpha=alpha[n], zorder=n)
#            ax[0].hist(self.qpt_norm, bins=50, weights=self.lo, histtype='stepfilled', linewidth=1.5, stacked=True, 
#                       color=self.color[3], linestyle = self.linestyle[n],alpha=alpha[n], zorder=2+n)
#            ax[0].hist(self.qpt_norm, bins=50, weights=self.acoustic, histtype='stepfilled', linewidth=1.5, stacked=True,
#                       color=self.color[1], linestyle = self.linestyle[n],alpha=alpha[n], zorder=6+n)
#            ax[0].hist(self.qpt_norm, bins=50, weights=self.to, histtype='stepfilled', linewidth=1.5, stacked=True,
#                       color=self.color[2], linestyle = self.linestyle[n],alpha=alpha[n], zorder=4+n)
            histzpr, binszpr = self.get_histogram(self.zpr_qpt)
      #      ax[0].plot(binszpr[:-1], histzpr, 'r')
#            spl_hist = UnivariateSpline(binszpr[15:-1], histzpr[15:], k=4)
#            xq = np.linspace(binszpr[15], binszpr[-1], 100)
#            ax[0].plot(xq, spl_hist(xq), 'y', lw=3)
            # CUmulative ZPR contribution

            # Qpoint density
            hist, bins = self.get_bz_histogram()  # FIX ME : this does not work very well.... should be removed

            line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50, ax=ax[1], color='gray', element='step', fill=True, 
                          kde=True, kde_kws={'cut':0}).get_lines()[0].get_data()
    #        sns.histplot(x=bins[:-1], weights=hist, bins=50, ax=ax[2], color='gray', element='step', fill=False, 
    #                     kde=True, kde_kws={'cut':0}, line_kws={'color':'r'})

            histbz, bins = self.get_histogram(self.wtq)
            print('histbz sum:', np.sum(histbz))
#            line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50,kde=True).get_lines()[0].get_data()
      #      ax[1].clear()
            ax[1].plot(line[0], line[1], 'k', linewidth=1.5)
    #        ax[1].plot(bins[15], histbz[15],'ok')
            # extend histogram data to the same lenght as the kde
#            longhist = np.zeros((200))
#            longhistzpr = np.zeros((200))
#            for i in range(50):
#                longhistzpr[4*i:4*(i+1)] = histzpr[i]
#                longhist[4*i:4*(i+1)] = histbz[i]
#            longcorr = line[1]/longhist
#            corr = line[1][::4]/histbz
    #        ax[0].plot(binszpr[:-1], corr*histzpr, 'g')
#            ax[0].plot(line[0], longcorr*longhistzpr, 'm')

#            sns.kdeplot(x=self.qpt_norm, weights=self.wtq, ax=ax[1], color=self.color[0], linestyle='solid', fill=False, cut=0)
    #        ax[2].hist(self.qpt_norm, bins=50, weights=self.wtq, histtype='step', linewidth=5.0, color=self.color[0], linestyle = self.linestyle[n])
            #self.shade_lowq(ax, 0.2)
    #        print(len(sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50,kde=True).get_lines()))

    #        print(ax[1].get_ylim())
            #xlim = ax[0].get_xlim()
            xlim = [0,1]
            ax[0].plot(np.linspace(xlim[0], xlim[1], 10), np.zeros((10)), 'k', zorder=-1, linewidth=0.75)
            for i in range(2):
                ax[i].set_xlim(xlim[0], xlim[1])
                ax[i].tick_params(axis='both', labelsize=12)
#            ax[0].set_ylim(-0.1, 1.505)

        if n+1  == self.nfile:

            # For CdTe
#            ax[0].annotate('',xy=(0.15,0.70), xycoords='data', xytext=(0.20, 0.90),arrowprops=dict(color='k',headlength=12,headwidth=8,width=0.2))
#            ax[0].annotate('',xy=(0.70,0.30), xycoords='data', xytext=(0.70, 0.50),arrowprops=dict(color='k',headlength=12,headwidth=8,width=0.2))
#            ax[0].text(0.10, 0.92, r'change of $m^\ast$', fontsize=20)
#            ax[0].text(0.52, 0.51, r'global $\varepsilon_{kn}$ reduction', fontsize=20)
            # For CdS
            ax[0].annotate('',xy=(0.15,2.20), xycoords='data', xytext=(0.20, 2.50),arrowprops=dict(color='k',headlength=12,headwidth=8,width=0.2))
            ax[0].annotate('',xy=(0.70,0.63), xycoords='data', xytext=(0.70, 1.20),arrowprops=dict(color='k',headlength=12,headwidth=8,width=0.2))
            ax[0].text(0.10, 2.58, r'change of $m^\ast$', fontsize=20)
            ax[0].text(0.52, 1.23, r'global $\varepsilon_{kn}$ reduction', fontsize=20)

            self.set_hist_labels(ax)

            self.set_legend(ax[0])
            
            if self.title:
                plt.suptitle(self.title, fontsize=12)

            if self.savefig:
                plt.savefig('{}_qpt_mode.pdf'.format(self.savefig), dpi=1200)
                os.system('open {}_qpt_mode.pdf'.format(self.savefig))
            else:
                plt.show()

    def plot_qpt_contribution(self, n, ax):

        alpha = [0.6, 1.0] 
#        ax.plot(self.qpt_norm, self.zpr_qpt, 'o')

        ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='stepfilled', alpha=alpha[n], color=self.color[0],  linestyle = self.linestyle[n], zorder=n)
#        ax[1].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, cumulative=True, density=True, histtype='step', linestyle = self.linestyle[n])
        #ax[1].hist(self.qpt_norm, bins=50, weights=self.wtq, histtype='stepfilled', linestyle = self.linestyle[n])

        if self.with_gkk:
            ax[1].hist(self.qpt_norm, bins=50, weights=self.gkk2, histtype='stepfilled', alpha=alpha[n], color=self.color[0],  linestyle = self.linestyle[n], zorder=n)

            line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50, ax=ax[2], color='gray', element='step', fill=True, 
                          kde=True, kde_kws={'cut':0}).get_lines()[0].get_data()
            ax[2].plot(line[0], line[1], 'k', linewidth=1.5)
        else:
            line = sns.histplot(x=self.qpt_norm, weights=self.wtq, bins=50, ax=ax[1], color='gray', element='step', fill=True, 
                          kde=True, kde_kws={'cut':0}).get_lines()[0].get_data()
            ax[1].plot(line[0], line[1], 'k', linewidth=1.5)

        if n+1 == self.nfile:
            self.set_hist_labels(ax)

            self.set_legend(ax[0])

            if self.title:
                plt.suptitle(self.title, fontsize=14)
            if self.savefig:
                plt.savefig('{}_qpt.pdf'.format(self.savefig), dpi=1200)
                os.system('open {}_qpt.pdf'.format(self.savefig))
            else:
                plt.show()

    def set_hist_labels(self, ax):

        fs = 18
        if self.cumulative:
            ax[2].set_xlabel(r'Q-point norm ($\frac{2\pi}{a}$)', fontsize=fs)
            ax[0].set_ylabel(r'Contribution to ZPR (meV)', fontsize=fs)
            if self.flist_labels:
                ax[1].set_ylabel(r'Cumulative contribution\\to ZPR / Total ZPR ({})'.format(self.flist_labels[0]), fontsize=12)
            else:
                ax[1].set_ylabel(r'Cumulative contribution\\to ZPR / Total ZPR', fontsize=fs)

    #        ax[2].set_ylabel(r'Number of q-points', fontsize=fs)
            ax[2].set_ylabel(r'BZ weight', fontsize=fs)

    #        ax[0].set_title(r'Contribution to ZPR vs $\vert q\vert$', fontsize=12)
    #        ax[1].set_title(r'Cumulative contribution vs $\vert q\vert$', fontsize=12)
    #        ax[2].set_title(r'Q-point distribution', fontsize=12)
        else:
            if self.with_gkk:
                ax[2].set_xlabel(r'Q-point norm ($\frac{2\pi}{a}$)', fontsize=fs)
                ax[0].set_ylabel(r'Contribution to ZPR (meV)', fontsize=fs)

                ax[1].set_ylabel(r'$|gkk|^2$ (arb. units)', fontsize=12)
                ax[2].set_ylabel(r'BZ weight', fontsize=fs)

        #        ax[0].set_title(r'Contribution to ZPR vs $\vert q\vert$', fontsize=12)
        #        ax[1].set_title(r'Cumulative contribution vs $\vert q\vert$', fontsize=12)
        #        ax[2].set_title(r'Q-point distribution', fontsize=12)
            else:
                ax[1].set_xlabel(r'Q-point norm ($\frac{2\pi}{a}$)', fontsize=fs)
                ax[0].set_ylabel(r'Contribution to ZPR (meV)', fontsize=fs)

        #        ax[2].set_ylabel(r'Number of q-points', fontsize=12)
                ax[1].set_ylabel(r'BZ weight', fontsize=fs)

        #        ax[0].set_title(r'Contribution to ZPR vs $\vert q\vert$', fontsize=12)
        #        ax[1].set_title(r'Cumulative contribution vs $\vert q\vert$', fontsize=12)
        #        ax[2].set_title(r'Q-point distribution', fontsize=12)


    def set_legend(self, ax):

        handles = []
        labels = ['Total', 'Acoustic', 'TO', 'LO']
        for i in range(4):
            handles.append(Line2D([0], [0], linewidth=4.0, color=self.color[i], label=labels[i]))

        if self.cumulative:
            artist = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.62, 0.75), ncol=4, fontsize=12)
            ax.add_artist(artist)
        else:
#            artist = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.50, 1.01), ncol=4, fontsize=16)
            artist = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.60, 0.85), ncol=4, fontsize=16,
                               columnspacing=1.1, handlelength=1.5)
            ax.add_artist(artist)

        if self.nfile > 1 and self.flist_labels:
            handles = []
            for i in range(self.nfile):
                handles.append(Line2D([0], [0], linewidth=3.0, color='k', label=self.flist_labels[i], linestyle=self.linestyle[i]))

            if self.cumulative:
                artist2 = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.81, 0.57), ncol=4, fontsize=12)
#            else:
#                artist2 = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.75, 0.70), ncol=4, fontsize=16)


    def get_cumulative_histogram(self, y):

        hist, bins = np.histogram(self.qpt_norm, bins=50, weights=y)
        cumul_hist = np.zeros_like(hist)

        for i in range(len(hist)):
            cumul_hist[i] = np.sum(hist[:i])

        cumul_hist = cumul_hist/cumul_hist[-1]*np.sum(y)/np.sum(self.ref_zpr)

        return cumul_hist, bins

    def get_bz_histogram(self):
        
        hist, bins = np.histogram(self.qpt_norm, bins=50, weights=self.wtq)
        trpz = np.trapz(y=hist, x=bins[:-1], axis=0)

        return hist/trpz, bins

    def get_histogram(self, y):

        hist, bins = np.histogram(self.qpt_norm, bins=50, weights=y)

        return hist, bins


def compute_contribution(

        qptcont_flist=None,
        qptmode_flist=None,
        qptpath_flist=None,
        qptcont_farr=None,
        qptmode_farr=None,
        qpt_weight=None,
        kpt_index=1,
        band_index=None,

        grid=True,
        path=False,
        mode=False,
        with_gkk=True,
        cumulative=True,
        ratio=True,
        title=None,
        verbose=False,
        bz_weight=True,
        figsize=(8, 8),
        
        *args, **kwargs):

    grid = grid
    path = path

    qptcont_flist = qptcont_flist
    qptmode_flist = qptmode_flist
    qptpath_flist = qptpath_flist
    qptcont_farr = qptcont_farr
    qptmode_farr = qptmode_farr

    qpt_weight_fname = qpt_weight
    mode = mode
    cumulative = cumulative
    ratio = ratio
    kpt_index = kpt_index-1
    bz_weight = bz_weight

    if not band_index:
        raise Exception('Must provide the band_index band index, starting at 1')

    band_index = band_index-np.ones_like(band_index)  # transfer to python indexing
    figsize = figsize
    title = title
    
    if grid:
        qptdata = QptContribution(

            qptcont_flist=qptcont_flist,
            qptmode_flist=qptmode_flist,
            qptcont_farr=qptcont_farr,
            qptmode_farr=qptmode_farr,

            qpt_weight_fname=qpt_weight_fname,
            aode=mode,
            cumulative=cumulative,
            ratio=ratio,

            kpt_index=kpt_index,
            band_index=band_index,
            bz_weight=bz_weight,

            figsize=figsize,
            title=title,
            verbose=verbose,

            **kwargs)

    elif path:
        qptdata = QptPathContribution(

            qptpath_flist=qptpath_flist,
            qpt_weight_fname=qpt_weight_fname,
            with_gkk = with_gkk,

            kpt_index=kpt_index,
            band_index=band_index,

            figsize=figsize,
            title=title,
            verbose=verbose,

            **kwargs)


#    if mode:
#        qptdata.plot_qpt_mode_contribution()
#    else:
#        qptdata.plot_qpt_contribution()
