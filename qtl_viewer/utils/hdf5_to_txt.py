import h5py
import rpy2.robjects as ro
import numpy as np

from qtl_viewer import h5_check


import argparse
import csv
import os
import sys

#
# datatypes
#

vstr = h5py.special_dtype(vlen=bytes)
dt_features = np.dtype([('feature_id', vstr), ('group_id', vstr), ('chrom', 'S2'), ('location', '<f8'), ('name', vstr), ('description', vstr)])
dt_markers = np.dtype([('marker_id', vstr), ('chrom', 'S2'), ('location', '<f8'), ('name', vstr), ('description', vstr)])
dt_strains = np.dtype([('strain_id', vstr), ('name', vstr), ('description', vstr)])
dt_samples = np.dtype([('sample_id', vstr), ('name', vstr), ('description', vstr)])
dt_factors = np.dtype([('factor_id', vstr), ('name', vstr), ('description', vstr)])



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert hdf5 to txt')

    parser.add_argument('-i', '--hdf', help='Input HDF file', required=True)
    parser.add_argument("-c", "--csv", help="CSV instead of TSV", dest="csv", action='store_true')

    args = vars(parser.parse_args())

    print str(args)

    # check the arguments
    h5_check.main(args['hdf'], True)






