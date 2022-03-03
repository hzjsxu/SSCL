#!/bin/python
# -*- coding:utf-8 -*-
###
# File: plot_HiC.py
# Author: jsxu
# Contact: <hzaujsxu@163.com>
# TDDO: plot heatmap and bw tracks of given regions.
###

import argparse
from scipy import ndimage
from math import floor, ceil
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import log2, max
from matplotlib.colors import LogNorm, LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from typing import get_args
import numpy as np
import pandas as pd
import pyBigWig as pb
import strawC
import cooler
import logging
import sys
import math
import matplotlib
matplotlib.use('Agg')

LOGLEVEL = logging.WARNING


def get_logger(name, file_=sys.stderr, level=LOGLEVEL):
    FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
    formatter = logging.Formatter(fmt=FORMAT)
    if isinstance(file_, str):
        handler = logging.FileHandler(file_)
    else:
        handler = logging.StreamHandler(file_)
    handler.setFormatter(formatter)
    log = logging.getLogger(name)
    log.addHandler(handler)
    log.setLevel(level)
    return log


log = get_logger(__name__)


class GenomeRange(object):
    """
    get the range region on the genome.
    """

    def __init__(self, *args):
        if len(args) == 1:
            if isinstance(args[0], GenomeRange):
                chrom, start, end = tuple(args[0])
            else:
                chrom, start, end = GenomeRange.parse_region_string(args[0])
        elif len(args) == 3:
            chrom, start, end = args
        else:
            raise ValueError("Inappropriate initial arguments."
                             "Correct example: range1 = GenomeRanges(\"chr1:1000-2000\") or "
                             "range1 = GenomeRange(\"chr1\", 1000, 2000)")

        if end < start:
            raise ValueError(f"Please check that the region end is larger than the region start. "
                             f"Values given: start: {start}, end: {end}")

        self.chrom = chrom
        self.start = start
        self.end = end

    @staticmethod
    def parse_region_string(region_string):
        """
        split a region into a tuple (chrom, start, end).
        """
        try:
            chrom, position = region_string.strip().split(':')
            for char in ",.;|!{}()":
                position = position.replace(char, "")
            position_list = position.split('-')
            region_start = int(position_list[0])
            region_end = int(position_list[1])
            return chrom, region_start, region_end

        except Exception:
            raise ValueError(f"Failed to parse region {region_string}, please check the region format "
                             f"should be like \"chr:start-end\".")

    @property
    def length(self):
        """
        return the length of region.
        """
        return self.end - self.start


class STRAW(object):
    """
    A wrap for straw Python API, for read .hic file.
    Parameters
    ----------
    path : str
        Path to the '.hic' file
    normalization : {False, no, 'VC', 'VC_SQRT', 'KR'}
        Method for the matrix normalization.
        default 'KR'
    binsize : int
        resolution of the data. for example 5000.
    >> test = STRAW(path=hic_file)
    >> test_mat = test.fetch('chr11:109413000-109513000')
    """

    def __init__(self, path, normalization='KR', binsize=25000):
        self.hic_file = path
        if normalization == "no" or normalization is False:
            normalization = 'NONE'
        elif normalization:
            normalization = normalization
        self.binsize = binsize
        self.normalization = normalization

    def fetch(self, genome_region1, genome_region2=None):
        """
        Return matrix np.ndarray.
        """
        if genome_region2 is None:
            genome_region2 = genome_region1
        #genome_region1 = to_gr(genome_region1)
        #genome_region2 = to_gr(genome_region2)
        straw_iter = self.__fetch_straw_iter(genome_region1, genome_region2)
        return self.__straw_to_matrix(straw_iter, genome_region1, genome_region2)

    def __fetch_straw_iter(self, genome_region1, genome_region2):
        slist = self.__fetch_straw_list_strawc(genome_region1, genome_region2)
        siter = [(r.binX, r.binY, r.counts) for r in slist]
        return siter

    def __fetch_straw_list_strawc(self, genome_region1, genome_region2):
        chr1loc = str(genome_region1).replace('-', ':')
        chr2loc = str(genome_region2).replace('-', ':')
        slist = []
        try:
            slist = strawC.strawC(
                self.normalization, self.hic_file, chr1loc, chr2loc, 'BP', self.binsize)
        except Exception as e:
            log.warning("Error occurred when reading the file with straw:")
            log.warning(str(e))
        return slist

    def __straw_to_matrix(self, straw_iter, genome_region1, genome_region2):
        is_cis = (genome_region1 == genome_region2)
        genome_region1 = to_gr(genome_region1)
        genome_region2 = to_gr(genome_region2)
        binlen1 = genome_region1.length // self.binsize + 1
        binlen2 = genome_region2.length // self.binsize + 1
        mat = np.zeros((binlen1, binlen2), dtype=np.float64)
        for rec in straw_iter:
            loc1, loc2, c = rec[0], rec[1], rec[2]
            bin1id = (loc1 - genome_region1.start) // self.binsize
            bin2id = (loc2 - genome_region2.start) // self.binsize
            if is_cis:
                mat[bin1id, bin2id] = c
                mat[bin2id, bin1id] = c
            else:
                if (0 <= bin1id < mat.shape[0]) and (0 <= bin2id < mat.shape[1]):
                    mat[bin1id, bin2id] = c
        return mat


def to_gr(obj):
    """
    convert obj to GenomeRange object.
    """
    if isinstance(obj, GenomeRange):
        return obj
    else:
        return GenomeRange(obj)


def read_bw(filename, region, nBins=700, type='mean'):
    """
    Todo:
        reading bigwig files.
    ------
    Parameters:
        filename: bw file.
        region: the format could be 'chr:start-end'
        nBins: the number of bins that region needs to be divided.
    ------
    Return:
        x_values, scores_per_bin
    """
    bw = pb.open(filename)
    region = to_gr(region)
    scores_per_bin = bw.stats(
        region.chrom, region.start, region.end, nBins=nBins, type=type)
    x_values = np.linspace(region.start, region.end, num=nBins)
    x_values = np.linspace(1, len(x_values), num=len(x_values))
    scores_per_bin = np.array(scores_per_bin)
    bw.close()
    return x_values, scores_per_bin


def plot_coverage(ax, x_values, scores, plot_type, size=1, color='red', negative_color='blue', alpha=1, grid=0):
    if grid:
        ax.grid(axis='y', zorder=0)
    if plot_type == 'line':
        if color == negative_color:
            ax.plot(x_values, scores, '-', linewidth=size,
                    color=color, alpha=alpha)
        else:
            log.warning(
                'Line plots with a different negative color might not look pretty.\n')
            pos_x_values = x_values.copy()
            pos_x_values[scores < 0] = np.nan
            ax.plot(pos_x_values, scores, '-',
                    linewidth=size, color=color, alpha=alpha)
            neg_x_values = x_values.copy()
            neg_x_values[scores >= 0] = np.nan
            ax.plot(neg_x_values, scores, '-', linewidth=size,
                    color=negative_color, alpha=alpha)

    elif plot_type == 'point':
        if color == negative_color:
            ax.plot(x_values, scores, '.', markersize=size,
                    color=color, alpha=alpha)
        else:
            pos_x_values = x_values.copy()
            pos_x_values[scores < 0] = np.nan
            ax.plot(pos_x_values, scores, '.',
                    markersize=size, color=color, alpha=alpha)
            neg_x_values = x_values.copy()
            neg_x_values[scores >= 0] = np.nan
            ax.plot(neg_x_values, scores, '.', markersize=size,
                    color=negative_color, alpha=alpha)

    else:
        if plot_type != 'fill':
            log.warning(
                'The plot type was not part of known types: {fill, line, points}.')
        if color == negative_color:
            ax.fill_between(x_values, scores, linewidh=0.1, color=color,
                            facecolor=color, alpha=alpha, interpolate=True)
        else:
            ax.fill_between(x_values, scores, where=scores >= 0, color=color,
                            facecolor=color, alpha=alpha, interpolate=True)
            ax.fill_between(x_values, scores, where=scores < 0, color=negative_color,
                            facecolor=negative_color, alpha=alpha, interpolate=True)
    return ax


def getopt():
    parser = argparse.ArgumentParser(
        usage='%(prog)s -f test.hic -o.test.pdf -r chr12:74358000-82302000 -s square --vMin 0.0 --vmax 99 -cm classic -bs 25000 \
                                    -bw test_1.bw test_2.bw test_3.bw --track_labels test_1 test_2 test_3 --plot_types line fill fill \
                                    --colors gray red yellow --negative_colors gray blue blue --vMinBigwig -1 -150 0 --vMaxBigwig 1 150 30',
        description='Plot genomic tracks on specified regions.',
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # parameter Group 1
    group_required = parser.add_argument_group(
        title='Required parameters',
        description='\n'.join([
            '-h, --help         Show this help message and exit',
            '-f, --hicFile      Input HiC contact map file, must be .hic',
            '-r, --region       The regions to plot, the format is "chr:start-end".',
            '-o, --outfile      The Plotted figure name'
        ])
    )
    group_required.add_argument(
        '-f', '--hicFile', type=str, help=argparse.SUPPRESS, required=True)
    group_required.add_argument(
        '-r', '--region', type=str, help=argparse.SUPPRESS, required=True)
    group_required.add_argument(
        '-o', '--outfile', type=str, help=argparse.SUPPRESS, required=True)

    # parameter Group 2
    group_optional = parser.add_argument_group(
        title='Optional parameters',
        description='\n'.join([
            '-s, --shape             HiC heatmap could be either "square" or "triangular". default: "square"',
            '--vMin                  Minimum percentile score value. Range from 0.0 to 100.0',
            '--vMax                  Maximum percentile score value. Range from 0.0 to 100.0',
            '--annoRegion            Region needs to annotate using rectangle.',
            '-cm, --colormap         The colorMap used in HiC heatmap.',
            '-bs, --binsize          The resolution (bp).',
            '-bw, --bwfiles          The (multiple) bigwig files used to plot tracks.',
            '--track_labels          Name of bigwig tracks.',
            '--plot_types            Type of bigwig tracks, choices: ["line" , "point", "fill"].',
            '--colors                Colors of positive values for bigwig tracks.',
            '--negative_colors       Colors of negative values for bigwig tracks.',
            '--vMinBigwig            Minimum score value for bigwig tracks.',
            '--vMaxBigwig            Maximum score value for bigwig tracks.'
        ])
    )
    group_optional.add_argument(
        '-h', '--help', action='help', help=argparse.SUPPRESS)
    group_optional.add_argument('-s', '--shape', nargs='?', choices=(
        'square', 'tri'), default='square', help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--vMin', type=float, default=0.0, help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--vMax', type=float, default=99.0, help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--annoRegion', type=str, help=argparse.SUPPRESS)
    group_optional.add_argument(
        '-cm', '--colormap', type=str, default='classic', help=argparse.SUPPRESS)
    group_optional.add_argument(
        '-bs', '--binsize', type=int, default=25000, help=argparse.SUPPRESS)
    group_optional.add_argument(
        '-bw', '--bwfiles', type=str, nargs='+', help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--track_labels', type=str, nargs='+', help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--plot_types', type=str, nargs='+', help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--colors', type=str, nargs='+', help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--negative_colors', type=str, nargs='+', help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--vMinBigwig', type=float, nargs='+', help=argparse.SUPPRESS)
    group_optional.add_argument(
        '--vMaxBigwig', type=float, nargs='+', help=argparse.SUPPRESS)
    return parser.parse_args()


if __name__ == "__main__":
    args = getopt()

    hic_tmp = STRAW(path=args.hicFile, binsize=args.binsize)
    hic_mat = hic_tmp.fetch(args.region)
    hic_mat = np.nan_to_num(hic_mat)
    length = hic_mat.shape[0]

    # you can define colormap yourself. If so, you need to modify paramater: '-cm', '--colormap'
    if args.colormap == "classic":
        cmap = LinearSegmentedColormap.from_list(
            'normal', [(1, 1, 1), (1, 0, 0)], N=1000)  # white -> red
    else:
        cmap = plt.get_cmap(args.colormap)

    vmin, vmax = np.percentile(
        hic_mat, args.vMin), np.percentile(hic_mat, args.vMax)

    left, bottom, width, height = 0.25, 0.5, 0.45, 0.45
    if args.shape == 'square':
        size_heatmap = [left, bottom, width, height]
        size_colorbar = [left + width + 0.02, bottom, width/20, height]
    elif args.shape == 'tri':
        size_heatmap = [left, bottom, width, height/2]
        size_colorbar = [left + width + 0.02, bottom, width/20, height/2]
        hic_mat = ndimage.rotate(hic_mat, 45, order=0,
                                 reshape=True, prefilter=False, cval=np.nan)
        length = hic_mat.shape[0]
    else:
        log.warning(
            f'Please ensure the HiC matrix shape is either "square" or "tri".')

    fig = plt.figure(figsize=(100, 100), facecolor='w', edgecolor='w')
    fig.subplots_adjust(hspace=0, wspace=0.02)

    # 1. add HiC heatmap.
    ax1 = fig.add_axes(size_heatmap)
    img = ax1.imshow(hic_mat, cmap=cmap, origin='upper', interpolation='nearest',
                     extent=(0, length, 0, length),
                     aspect='auto', vmax=vmax, vmin=vmin)  # extent: (left, right, bottom, top)
    if args.annoRegion:
        anno_region = to_gr(args.annoRegion)
        anno_start, anno_end = int(anno_region.start / args.binsize), int(anno_region.end / args.binsize)
        anno_width, anno_height = anno_end - anno_start, anno_end - anno_start
        ax1.add_patch(
            patches.Rectangle((anno_start, length-anno_end), anno_width, anno_height, edgecolor='black', fill=False)
            )
    ax1.set_ylabel(args.region, fontsize=12)
    ax1.set_xlim([0, length])
    ax1.set_ylim([0, length])
    if args.shape == 'tri':
        ax1.set_ylim([length/2, length])
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklines(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    plt.setp(ax1.get_yticklines(), visible=False)

    # 2. add colorbar for HiC heatmap.
    colorAx = fig.add_axes(size_colorbar)
    cbar = plt.colorbar(img, cax=colorAx, orientation='vertical')
    cbar.set_label('Normalized interaction counts ' +
                   str(round(vmax, 2)), fontsize=12)

    # 3. add bigwig tracks: for loop.
    if args.bwfiles:
        len_bws = len(args.bwfiles)
        step_h = (1 - bottom)/len_bws/3
        bottom_h = bottom - step_h
        for i in range(len(args.bwfiles)):
            bw_file = args.bwfiles[i]
            plot_type = args.plot_types[i]
            color = args.colors[i]
            negative_color = args.negative_colors[i]
            vminBigwig, vmaxBigwig = args.vMinBigwig[i], args.vMaxBigwig[i]
            track_label = args.track_labels[i]

            ax = fig.add_axes([left, bottom_h, width, step_h-step_h/5])
            x_values, y_scores = read_bw(bw_file, args.region, nBins=length)
            ax = plot_coverage(ax, x_values, y_scores, plot_type=plot_type,
                               color=color, negative_color=negative_color)
            ax.hlines(0, 0, len(x_values), colors="gray", linestyles="dashed")
            ax.set_ylim([vminBigwig, vmaxBigwig])
            ax.set_yticks([vminBigwig, vmaxBigwig])
            ax.set_xlim([0, length])
            ax.set_xticks([])
            ax.set_ylabel(track_label)
            ax.get_yaxis().set_label_coords(-0.125, 0.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            bottom_h -= step_h
    fig.savefig(args.outfile)
