#!/store/jsxu/biosoft/minicondas3/envs/biosoft/bin/python
# -*- coding:utf-8 -*-
###
# File: calc_HiC_feature.py
# Created Date: 2021-06-21
# Author: jsxu
# Contact: <hzaujsxu@163.com>
# Last Modified: Wednesday June 23rd 2021 3:18:28 pm
#
###

import numpy as np
import time
import sys
from scipy import sparse
from multiprocessing import Pool
import cooler
from cooler.core import CSRReader
from cooler.util import open_hdf5
import strawC
import argparse

def getopt():
    parser = argparse.ArgumentParser(
        usage           = '%(prog)s -i test.hic -o.test.bdg -r 25000 -c Chr1 -g test.chrom.size -f insulation_score',
        description     = 'Calculate the directionality index or insulation score from HiC contact map file.',
        add_help        = False,
        formatter_class = argparse.RawDescriptionHelpFormatter
        )

    opt = parser.add_argument_group(
        title = 'optional arguments',
        description = '\n'.join([
            '-h, --help            show this help message and exit',
            '-i, --input           Input HiC contact map file, must be .hic or .cool file',
            '-o, --outfile         output bdg-format di/insulation score file',
            '-r, --resolution      Resolution (bp) (default: 25000)',
            '-w, --window_size     Length of upstream array and downstream array (default: 20)',
            '    --ignore_diags    The number of diagonals to ignore (default: 1)',
            '-c, --chromosome      Chromosome name to calculate (default: all)',
            '    --norm            Normalization method used only for .hic files, could be one of NONE, KR, VC-SQRT (default: KR)',
            '-g, --genome          The chromosome size file (tab-seperated)',
            '-f, --feature         could be one of standard_DI, adaptive_DI or insulation_score',
            '-t, --threads         Number of threads (default: 1)'
        ])
    )

    opt.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)
    opt.add_argument('-i', '--input', type=str, help=argparse.SUPPRESS, required=True)
    opt.add_argument('-o', '--outfile', type=str, help=argparse.SUPPRESS, required=True)
    opt.add_argument('-r', '--resolution', type=int, help=argparse.SUPPRESS, required=True, default=25000)
    opt.add_argument('-w', '--window_size', type=int, help=argparse.SUPPRESS, default=20)
    opt.add_argument('--ignore_diags', type=int, help=argparse.SUPPRESS, default=1)
    opt.add_argument('--norm', type=str, help=argparse.SUPPRESS, default='KR')
    opt.add_argument('-c', '--chromosome', type=str, help=argparse.SUPPRESS, default='all')
    opt.add_argument('-f', '--feature', type=str, help=argparse.SUPPRESS, required=True)
    opt.add_argument('-g', '--genome', type=str, help=argparse.SUPPRESS, required=True)
    opt.add_argument('-t', '--threads', type=int, help=argparse.SUPPRESS, default=1)

    return parser.parse_args()

## extract matrix from .cool/.hic file
def extract_matrix(contact_map_file, chrom, res, norm='KR', symmetric=False):
    """[summary]
    Args:
        NOTE: the paramater res and norm is specified for .hic files; as the single resolution and normalized cool file is needed.
        contact_map_file ([str]): path to contact map file, could be .cool/.hic
        chrom ([str]): chromosome name.
        res ([int]): resolution (bp), for .hic files.
        norm ([str]): normalization method, default: KR; for .hic files.
        symmetric ([bool]): if True, return symmetric matrix; otherwise return a triangular matrix.
    """
    if contact_map_file.endswith('cool'):
        co = cooler.Cooler(contact_map_file)
        with open_hdf5(co.store) as h5:
            root = h5[co.root]
            chrom_offset = root['indexes']['chrom_offset']
            cid = co._chromids[chrom]
            row_st, row_ed = chrom_offset[cid], chrom_offset[cid + 1]
            col_st, col_ed = row_st, row_ed
            reader = CSRReader(
                h5=root,
                field='count',
                max_chunk=500000000000000000,
            )
            x, y, counts = reader.query(row_st, row_ed, col_st, col_ed)
            mat = sparse.csr_matrix(
                (counts, (x - row_st, y - col_st)),
                shape=(row_ed - row_st, col_ed - col_st)
            )

    elif contact_map_file.endswith('hic'):
        tmp = strawC.strawC(norm, contact_map_file, chrom, chrom, 'BP', res)
        trip_mat = np.array([[x.binX, x.binY, x.counts] for x in tmp])
        row = (trip_mat[:, 0] / res).astype(int)
        col = (trip_mat[:, 1] / res).astype(int)
        row_nmax = np.max(row)
        col_nmax = np.max(col)
        if row_nmax != col_nmax:
            nmax = max(row_nmax, col_nmax)
            trip_mat = np.vstack((trip_mat, np.array([nmax*res, nmax*res, 0])))
            row = (trip_mat[:, 0] / res).astype(int)
            col = (trip_mat[:, 1] / res).astype(int)
        counts = trip_mat[:, 2]
        mat = sparse.csr_matrix((counts, (row, col)))

    else:
        raise ValueError(
            'Contact map file should be .cool or .hic !'
        )

    mat = mat.toarray()
    mat[np.isnan(mat)] = 0

    if symmetric:
        mat += mat.T - np.diag(mat.diagonal())

    return mat.astype(np.float32)


def mask_array(mask, *args):
    """[Mask all ndarray in args with a given Boolean array.]
    Args:
        mask ([np.ndarray]): [Boolean array where desired values are marked with True.]
        args ([tuple]): [tuple of np.ndarray. Masking will be applied to each ndarray.]
    """
    for mat in args:
        if isinstance(mat, (tuple, list)):
            yield tuple(mask_array(mask, *mat))
        else:
            if len(mat.shape) == 1:
                yield mat[mask]
            else:
                yield mat[:, mask]

def di_score(matrix, window_size=20, ignore_diags=1, method="standard"):
    """Compute directionality index of a given 2d ndarray.\n
    For each bin in the main digonal, directionality index is calculated based on two arrays with length of windowsize: \n
    The upwards(upstream) vertical array start from this bin and the eastwards(downstream) horizontal array start from this bin.\n
    See function listed in 'DI_METHOD_MAP' for detailed description.
    :param matrix: np.ndarray/sparse.csr_matrix. The input matrix to calculate di score.
    :param window_size: int. length of upstream array and downstream array.
    :param ignore_diags: int. The number of diagonals to ignore.
    :param method: str. Method for computing directionality index. 'standard' and 'adaptive' are supported by now.
    :return: np.ndarray. Rerturn directionality index array.
    """
    def standard_di(up, down):
        """Compute directionality index described in:
        Jesse R.Dixon 2012. Topological domains in mammalian genomes identified by analysis of chromatin interactions.
        up: the number of reads that map from a given 40kb bin to the upstream 2M
        down: the number of reads that map from the same 40kb bin to the downstream 2M
        """
        np.seterr(divide='ignore', invalid='ignore')
        up = up.sum(axis=1)
        down = down.sum(axis=1)
        expected = (up + down) / 2.0
        standard_DI = (np.sign(down - up) * ((up - expected) ** 2 + (down - expected) ** 2)/expected)
        return standard_DI

    def adap_di(up, down):
        """Compute directionality index described in:\n
        Xiao-Tao Wang 2017.HiTAD: Detecting the structural and functional hierarchies of topologically associating domains from chromatin interactions.
        """
        np.seterr(divide='ignore', invalid='ignore')
        window_size = up.shape[1]
        mean_up = up.mean(axis=1)
        mean_down = down.mean(axis=1)
        var_up = np.square(up - mean_up[:, None]).sum(axis=1) ## mean_up[:, None] 增加一个维度
        var_down = np.square(down - mean_down[:, None]).sum(axis=1)
        adap_DI = (mean_up - mean_down) / np.sqrt((var_up + var_down) / (window_size * (window_size - 1)))
        return adap_DI

    method_map = {'standard': standard_di,
                  'adaptive': adap_di
                  }
    chrom_len = matrix.shape[0]
    x, y = matrix.nonzero()
    dis = y - x

    if isinstance(window_size, int):
        max_len = ignore_diags + window_size
        available = ((dis >= ignore_diags) & (dis < max_len) & (x >= max_len - 1) & (y <= chrom_len - max_len))
        x, y = mask_array(available, x, y)
        values = np.array(matrix[x, y]).ravel()
        x, y, values = mask_array(~np.isnan(values), x, y, values)
        dis = y - x
        contacts_up = np.zeros((chrom_len, window_size), dtype=matrix.dtype)
        contacts_down = np.zeros((chrom_len, window_size), dtype=matrix.dtype)
        for shift in range(ignore_diags, max_len):
            window_pos = shift - ignore_diags
            mask = dis == shift
            tmp_x, tmp_values = mask_array(mask, x, values)
            contacts_down[tmp_x, window_pos] = tmp_values
            contacts_up[tmp_x + shift, window_pos] = tmp_values
        contacts_up[:max_len, :] = 0
        contacts_down[:max_len, :] = 0

    else:
        raise ValueError('window_size should be an interger.')

    DI = method_map[method](contacts_up, contacts_down)

    return DI


def insulation_score(matrix, window_size=20, ignore_diags=1, normalize=True):
    """
    Calculate insulation score of a given 2d ndarray (chromosome) described in:
        Emily Crane 2015. Condensin-driven remodelling of X chromosome topology during dosage compensation.
    :param matrix: np.ndarray/scipy.sparse.csr_matrix. Interaction matrix representing a hic contacts.
    :param window_size: int. Diameter of square in which contacts are summed along the diagonal.
    :param ignore_diags: int. Number of diagonal to ignore, This values should be >= 1 which means ignore main diagonal.
    :param normalize: bool. If normalize the insulation score with log2 ratio of insu_score and mean insu_score.
    :return: np.ndarray. Return (chrom_len,) array.
    """
    chrom_len = matrix.shape[0]
    if isinstance(matrix, sparse.csr_matrix):
        x, y = matrix.nonzero()
        dis = y - x
        available = (dis >= ignore_diags) & (dis <= 2 * window_size)
        x, y = mask_array(available, x, y)
        values = np.array(matrix[x, y]).ravel()
        x, y, values = mask_array(~np.isnan(values), x, y, values)
        dis = y - x
        insu = np.zeros(chrom_len, dtype=matrix.dtype)
        insu[:window_size] = np.nan
        insu[chrom_len - window_size:] = np.nan
        counts = np.zeros(chrom_len, dtype=np.int)
        diag_dict = {}
        for i in range(window_size):
            for j in range(window_size):
                _dis = (abs(j - i) + (min(j, i) + 1) * 2)
                if diag_dict.get(_dis) is None:
                    mask = dis == (abs(j - i) + (min(j, i) + 1) * 2)
                    tmp_x, tmp_value = mask_array(mask, x, values)
                    diag_dict[_dis] = (tmp_x, tmp_value)
                else:
                    tmp_x, tmp_value = diag_dict[_dis]
                x_index = tmp_x + j + 1
                insu[x_index] += tmp_value
                counts[x_index] += 1
        counts[:window_size] = 0
        counts[chrom_len - window_size:] = 0
        insu /= counts

    elif isinstance(matrix, np.ndarray):
        insu = np.full(chrom_len, np.nan, dtype=matrix.dtype)
        counts = np.zeros(chrom_len, dtype=matrix.dtype)
        diamond_mask = np.full((window_size, window_size), True)
        diamond_mask = np.triu(diamond_mask, ignore_diags - window_size)
        for row in range(window_size, chrom_len - window_size):
            sub_mat = matrix[row - window_size:row, row + 1: row + window_size + 1][diamond_mask]
            insu[row] = np.nanmean(sub_mat)
            counts[row] = np.sum(~np.isnan(sub_mat))

    else:
        raise ValueError('matrix should be sparse.csr_matrix or np.ndarray.')

    if normalize:
        insu = np.log2(insu / np.nanmean(insu))

    return insu

def get_chrom_score(score, chrom, res):
    """
    param score: N*1 array, could be di_score or insulation_score array.
    param chroms: chromosome.
    param res: resolution.
    note: remove the last bin to prevent mistake from converting bedGraph to bw files using bedGraph.
            i.e. End coordinate 217475000 bigger than Chr1 size of 217471166 line 8699 of Brain.adaptive_DI.bdg
    """
    size = score.shape[0] - 1
    chroms = np.transpose([np.array([chrom] * size, dtype=str)])
    start = np.transpose([np.arange(0, size) * res])
    end = np.transpose([np.arange(1, size + 1) * res])
    scores = np.transpose([score[0:size]])
    return np.hstack((chroms, start, end, scores))

def load_genome(genomesize_file):
    tmp = np.loadtxt(genomesize_file, dtype=str)
    chrom_names = [x[0] for x in tmp]
    return chrom_names

def main():
    args = getopt()

    if args.chromosome == 'all':
        chroms = load_genome(args.genome)
    else:
        chroms = args.chromosome.split(' ')

    scores = []
    for chrom in chroms:
        chrom_mat = extract_matrix(args.input, chrom, args.resolution,
                                   norm=args.norm, symmetric=True)
        print(time.strftime('%Y-%m-%d %a %H:%M:%S', time.localtime()),
                        f': ****** Performing calculation on {chrom} ******')

        if args.feature == 'standard_DI':
            res_score = di_score(chrom_mat, ignore_diags=args.ignore_diags,
                             window_size=args.window_size, method='standard')
        elif args.feature == 'adaptive_DI':
            res_score = di_score(chrom_mat, ignore_diags=args.ignore_diags,
                             window_size=args.window_size, method='adaptive')
        elif args.feature == 'insulation_score':
            res_score = insulation_score(chrom_mat, window_size=args.window_size,
                                    ignore_diags=args.ignore_diags, normalize=True)
        else:
            print(f'Please select one feature from standard_DI, adaptive_DI and insulation_score.')
            sys.exit(0)
        score = get_chrom_score(res_score, chrom, args.resolution)
        scores.append(score)
    final_score = np.concatenate(scores, axis=0)
    np.savetxt(args.outfile, final_score, fmt='%s', delimiter='\t')
    return

if __name__ == "__main__":
    main()
