#! usr/bin/env python

import numpy as np
from zprfile import ZPRfile
import matplotlib.pyplot as plt
from constants import ha_to_ev
from matplotlib.lines import Line2D
from matplotlib import rc
import matplotlib.gridspec as gridspec
import os
from colorpalettes import bright, vibrant 
import seaborn as sns
from scipy.interpolate import UnivariateSpline
rc('text', usetex=True)
rc('font', family='sans-serif')

class QptContribution(object):

    def __init__(self,

                 qptcont_flist=None,
                 qptmode_flist=None,
                 qpt_weight_fname=None,
                 mode=False,
                 cumulative=True,
                 band_index=None,
                 kpt_index=0,

                 figsize=(8, 8),
                 title=None,
                 verbose=False,
                 color = [bright['blue'], bright['green'], bright['yellow'], bright['red']],
                 linestyle = ['solid', 'dashed', 'dotted', 'dashdot'],
                 savefig=None,
                 flist_labels=None,

                 **kwargs):

        self.mode = mode
        self.cumulative = cumulative
        self.band_index = band_index
        self.qpt_weight_fname = qpt_weight_fname

        self.figsize = figsize
        self.title = title
        self.verbose = verbose

        self.color = color
        self.linestyle = linestyle
        self.savefig = savefig
        self.flist_labels = flist_labels

        self.nfile = len(qptcont_flist)

        if not qpt_weight_fname:
            raise Exception('Must provide a file for qpt weights.')
        else:
            self.set_wtq(qpt_weight_fname)

        if self.mode and len(qptmode_flist) != self.nfile:
            raise Exception('''For mode decomposition, qptcont_flist and qptmode_flist should have
                            have the same lenght, but have {} and {}.'''.format(nfile, len(self.qptmode_flist)))

        if self.mode and len(self.color) < 4:
            raise Exception('For Mode decomposition, custom color list must contain 4 entries.')

        if len(self.linestyle) < self.nfile:
            raise Exception('For {} files in flist, linestyle need {} entries, but has only {}.'.format(
                            self.nfile, self.nfile, len(self.linestyle)))

        if len(self.band_index) != self.nfile:
            raise Exception('Band_index must have {} entries, but has {}.'.format(self.nfile, len(self.band_index)))

        for n in range(self.nfile):
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

            if mode:
                if n == 0:
                    if self.cumulative:
                        fig, ax = plt.subplots(3, 1, squeeze=True, sharex=True, figsize=self.figsize)
                    else:
                        fig = plt.figure(figsize=self.figsize)
                        spec = gridspec.GridSpec(ncols=1, nrows=3, figure=fig)
                        ax = [fig.add_subplot(spec[0:2,0]), fig.add_subplot(spec[2,0])]
#                        fig, ax = plt.subplots(2, 1, squeeze=True, sharex=True, figsize=self.figsize)
                    self.ref_zpr = np.sum(self.zpr_qpt)
                self.plot_qpt_mode_contribution(n, ax)
            else:
                if n == 0:
                    fig, ax = plt.subplots(3, 1, squeeze=True, sharex=True, figsize=self.figsize)
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

        alpha = [0.5, 1.0]

        if self.cumulative:
    #        ax.plot(self.qpt_norm, self.zpr_qpt, 'o')
            # ZPR contribution
            ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='step', linewidth=1.5, stacked=True,
                       color=self.color[0], linestyle = self.linestyle[n])
    #        ax[0].hist(self.qpt_norm, bins=50, weights=self.acoustic, histtype='step', linewidth=1.5, stacked=True,
    #                   color=self.color[1], linestyle = self.linestyle[n])
    #        ax[0].hist(self.qpt_norm, bins=50, weights=self.to, histtype='step', linewidth=1.5, stacked=True,
    #                   color=self.color[2], linestyle = self.linestyle[n])
    #        ax[0].hist(self.qpt_norm, bins=50, weights=self.lo, histtype='step', linewidth=1.5, stacked=True, 
    #                   color=self.color[3], linestyle = self.linestyle[n])
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
            ax[0].hist(self.qpt_norm, bins=50, weights=self.lo, histtype='stepfilled', linewidth=1.5, stacked=True, 
                       color=self.color[3], linestyle = self.linestyle[n],alpha=alpha[n], zorder=2+n)
            ax[0].hist(self.qpt_norm, bins=50, weights=self.acoustic, histtype='stepfilled', linewidth=1.5, stacked=True,
                       color=self.color[1], linestyle = self.linestyle[n],alpha=alpha[n], zorder=6+n)
            ax[0].hist(self.qpt_norm, bins=50, weights=self.to, histtype='stepfilled', linewidth=1.5, stacked=True,
                       color=self.color[2], linestyle = self.linestyle[n],alpha=alpha[n], zorder=4+n)
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
    #        ax[0].set_ylim(-0.1, 1.505)

        if n+1  == self.nfile:
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

#        ax.plot(self.qpt_norm, self.zpr_qpt, 'o')

        ax[0].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, histtype='step', linestyle = self.linestyle[n])
        ax[1].hist(self.qpt_norm, bins=50, weights=self.zpr_qpt, cumulative=True, density=True, histtype='step', linestyle = self.linestyle[n])
        ax[2].hist(self.qpt_norm, bins=50, weights=self.wtq, histtype='step', linestyle = self.linestyle[n])

        if n+1 == self.nfile:
            self.set_hist_labels(ax)

            if self.title:
                plt.suptitle(self.title, fontsize=14)
            if self.savefig:
                plt.savefig('{}_qpt.pdf'.format(self.savefig), dpi=1200)
                os.system('open {}_qpt.pdf'.format(self.savefig))
            else:
                plt.show()

    def set_hist_labels(self, ax):

        fs = 16
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
            handles.append(Line2D([0], [0], linewidth=3.0, color=self.color[i], label=labels[i]))

        if self.cumulative:
            artist = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.62, 0.75), ncol=4, fontsize=12)
            ax.add_artist(artist)
        else:
            artist = ax.legend(handles=handles, loc=8, bbox_to_anchor=(0.50, 1.03), ncol=4, fontsize=16)
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
        qpt_weight=None,
        kpt_index=1,
        band_index=None,

        mode=False,
        cumulative=True,
        title=None,
        verbose=False,
        figsize=(8, 8),
        
        *args, **kwargs):

    qptcont_flist = qptcont_flist
    qptmode_flist = qptmode_flist
    qpt_weight_fname = qpt_weight
    mode = mode
    cumulative = cumulative
    kpt_index = kpt_index-1

    if not band_index:
        raise Exception('Must provide the band_index band index, starting at 1')

    band_index = band_index-np.ones_like(band_index)  # transfer to python indexing
    figsize = figsize
    title = title
    qptdata = QptContribution(

        qptcont_flist=qptcont_flist,
        qptmode_flist=qptmode_flist,
        qpt_weight_fname=qpt_weight_fname,
        mode=mode,
        cumulative=cumulative,

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
