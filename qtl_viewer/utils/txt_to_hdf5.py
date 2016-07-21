import h5py
import rpy2.robjects as ro
import numpy as np

import argparse
import csv
from collections import OrderedDict
import logging
import os
import socket
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


#
#  utilities
#
def ffilename(filename):
    return os.path.abspath(filename)



#
# features
#

def create_features(h5_root, filename, tsv=True):
    logging.debug('{}/features: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = '\t' if tsv else ','

    features = []
    # rU is universal newline mode
    with open(filename, 'rU') as input:
        freader = csv.reader(input, delimiter=delim)
        for row in freader:
            features.append((row[0], row[1], row[2], row[3] if len(row[3]) > 0 else None, row[4] if len(row) > 4 else None, row[5] if len(row) > 5 else None))

    np_features = np.array(features, dtype=dt_features)
    np_features = np.array(features, dtype=dt_features)
    h5_root.create_dataset('features', data=np_features)

    logging.debug('{}/features: Created: {:,}'.format(h5_root.name, len(features)))

    return len(features)


#
# markers
#

def create_markers(h5_root, filename, tsv=True):
    logging.debug('{}/markers: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = '\t' if tsv else ','

    markers = []
    with open(filename, 'rb') as input:
        freader = csv.reader(input, delimiter=delim)
        for row in freader:
            markers.append((row[0], row[1], row[2], row[3] if len(row) > 3 else None, row[4] if len(row) > 4 else None))

    np_markers = np.array(markers, dtype=dt_markers)
    h5_root.create_dataset('markers', data=np_markers)

    logging.debug('{}/markers: Created: {:,}'.format(h5_root.name, len(markers)))
    return len(markers)


#
# lod
#
def create_lod(h5_root, filename, tsv=True):
    logging.debug('{}/lod/lod: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = r"\s+" if tsv else ','

    h5_root.create_group('lod')
    if tsv:
        np_lod = np.loadtxt(filename)
    else:
        np_lod = np.loadtxt(filename, delimiter=',')

    h5_root.create_dataset('lod/lod', data=np_lod)

    logging.debug('{}/lod/lod: Created: ({:,} x {:,})'.format(h5_root.name, np_lod.shape[0], np_lod.shape[1]))

    return np_lod.shape


#
# strains
#
def create_strains(h5_root, filename, tsv=True):
    logging.debug('{}/coef/strains: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = '\t' if tsv else ','

    strains = []
    with open(filename, 'rb') as input:
        freader = csv.reader(input, delimiter=delim)
        for row in freader:
            strains.append((row[0], row[1], row[2]))

    np_strains = np.array(strains, dtype=dt_strains)
    h5_root.create_dataset('coef/strains', data=np_strains)

    logging.debug('{}/coef/strains: Created: {:,}'.format(h5_root.name, len(strains)))
    return len(strains)


#
# coef
#
def create_coef(h5_root, filenames, tsv=True):
    '''
    coef file must end with the sample_id
    '''
    logging.debug('{}/coef/coef'.format(h5_root.name))

    delim = '\t' if tsv else ','

    strains_processed = OrderedDict()
    files_wo_strains = list(filenames)

    #print files_wo_strains

    strains = h5_root['coef/strains']

    for strain in strains:
        #print strain
        # find the file_name
        for f in filenames:
            if f.endswith('_{}.txt'.format(strain[0])):
                strains_processed[strain[0]] = f
                try:
                    files_wo_strains.remove(f)
                except:
                    pass
                break

    strains_wo_files = list(set(strains_processed.keys()).difference(set(s[0] for s in strains)))

    if len(strains_wo_files) != 0:
        logging.debug('There are strain ids without files: {}'.format(strains_wo_files))
        sys.exit()

    if len(files_wo_strains) > 0:
        logging.debug('There are files listed with --coef that contain no strain ids: {}'.format(files_wo_strains))
        sys.exit()

    ### TODO: SLOW

    num_features = h5_root['features'].shape[0]
    num_markers = h5_root['markers'].shape[0]
    num_strains = strains.shape[0]


    np_array = np.ndarray(shape=(num_features, num_strains, num_markers))

    strain_num = 0
    for strain, filename in strains_processed.iteritems():
        logging.debug('{}/coef/coef: Strain ID: {} Using {}'.format(h5_root.name, strain, ffilename(filename)))

        with open(filename, 'rb') as input:
            freader = csv.reader(input, delimiter=delim)
            rnum = 0
            for row in freader:
                cnum = 0
                for val in row:
                    np_array[rnum][strain_num][cnum] = val
                    cnum += 1
                rnum += 1

        strain_num += 1

        logging.debug('{}/coef/coef: Generated for Strain ID: {}'.format(h5_root.name, strain))

    h5_root.create_dataset('coef/coef', data=np_array)

    logging.debug('{}/coef/coef: Created: ({:,} x {:,} x {:,})'.format(h5_root.name, num_features, num_strains, num_markers))

    return np_array.shape


#
# samples
#
def create_samples(h5_root, filename, tsv=True):
    logging.debug('{}/samples: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = '\t' if tsv else ','

    samples = []
    with open(filename, 'rb') as input:
        freader = csv.reader(input, delimiter=delim)
        for row in freader:
            samples.append((row[0], row[1], row[2]))

    np_samples = np.array(samples, dtype=dt_samples)
    h5_root.create_dataset('samples', data=np_samples)

    logging.debug('{}/samples: Created: {:,}'.format(h5_root.name, len(samples)))
    return len(samples)


#
# factors
#
def create_factors(h5_root, filename, tsv=True):
    logging.debug('{}/phenotypes/factors: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = '\t' if tsv else ','

    factors = []
    with open(filename, 'rb') as input:
        freader = csv.reader(input, delimiter=delim)
        for row in freader:
            factors.append((row[0], row[1], row[2]))

    np_factors = np.array(factors, dtype=dt_factors)
    h5_root.create_dataset('phenotypes/factors', data=np_factors)

    logging.debug('{}/phenotypes/factors: Created: {:,}'.format(h5_root.name, len(factors)))
    return len(factors)


#
# phenotypes
#
def create_phenotypes(h5_root, filename, tsv=True):
    logging.debug('{}/phenotypes/phenotypes: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = '\t' if tsv else ','

    phenotypes = []
    max_columns = 0
    with open(filename, 'rb') as input:
        freader = csv.reader(input, delimiter=delim)
        for row in freader:
            phenotypes.append(tuple(row))
            max_columns = max(max_columns, len(row))

    np_phenotypes = np.array(phenotypes)
    h5_root.create_dataset('phenotypes/phenotypes', data=np_phenotypes)

    logging.debug('{}/phenotypes/phenotypes: Created: ({:,} x {:,})'.format(h5_root.name, len(phenotypes), max_columns))
    return (len(phenotypes), max_columns)


#
# genotypes
#
def create_genotypes(h5_root, filename, tsv=True):
    logging.debug('{}/genotypes/genotypes: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = '\t' if tsv else ','

    genotypes = []
    max_columns = 0
    with open(filename, 'rb') as input:
        freader = csv.reader(input, delimiter=delim)
        for row in freader:
            genotypes.append(tuple(row))
            max_columns = max(max_columns, len(row))

    np_genotypes = np.array(genotypes)
    h5_root.create_dataset('genotypes/genotypes', data=np_genotypes)

    logging.debug('{}/genotypes/genotypes: Created: ({:,} x {:,})'.format(h5_root.name, len(genotypes), max_columns))
    return (len(genotypes), max_columns)


#
# expression
#
def create_expression(h5_root, filename, tsv=True):
    logging.debug('{}/expression/expression: Using {}'.format(h5_root.name, ffilename(filename)))

    delim = '\t' if tsv else ','

    expression = []
    max_columns = 0
    with open(filename, 'rb') as input:
        freader = csv.reader(input, delimiter=delim)
        for row in freader:
            expression.append(tuple(row))
            max_columns = max(max_columns, len(row))

    np_expression = np.array(expression)
    h5_root.create_dataset('expression/expression', data=np_expression)

    logging.debug('{}/expression/expression: Created: ({:,} x {:,})'.format(h5_root.name, len(expression), max_columns))
    return (len(expression), max_columns)


if __name__ == '__main__':
    hostname = socket.gethostname()
    logging.basicConfig(format='{}: %(asctime)s: %(message)s'.format(hostname), datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)

    parser = argparse.ArgumentParser(description='Convert txt files to hdf5')

    parser.add_argument('-o', '--out', help='Output HDF file', required=True)
    parser.add_argument('-d', '--dataset', help='ID of the dataset', required=True)
    parser.add_argument('-f', '--features', help='features file, think genes', required=True)
    parser.add_argument('-l', '--lod', help='lod score file, features (rows) by markers (columns)', required=True)
    parser.add_argument('-m', '--markers', help='markers file, think snps', required=True)

    parser.add_argument('--dfA', help='Degress of freedom autosomes', type=int, required=True)
    parser.add_argument('--dfX', help='Degrees of freedom X', type=int, required=True)

    parser.add_argument("-c", "--csv", help="CSV instead of TSV", dest="csv", default=False, action='store_true')
    parser.add_argument('--datasetname', help='Name of the dataset, defaults to value of -d', required=False)

    # effect plot
    parser.add_argument('--coef', help='coef files, 1 per strain, features (rows) by markers (columns)', nargs='*', required=False)
    parser.add_argument('--strains', help='strains file', required=False)

    # factorial viewer stuff
    parser.add_argument('--samples', help='samples file', required=False)
    parser.add_argument('--factors', help='factors file, think diet, tissue, etc', required=False)
    parser.add_argument('--phenotypes', help='phenotypes file, samples (rows) by factors (columns)', required=False)
    parser.add_argument('--genotypes', help='genotypes file, markers (rows) by samples (columns)', required=False)
    parser.add_argument('--expression', help='expression data file, feature (rows) by samples (columns)', required=False)

    #print str(parser.parse_args())

    args = vars(parser.parse_args())

    #print str(args)

    # check the arguments

    h5_file = h5py.File(args['out'], "a")

    try:
        del h5_file[args['dataset']]
    except:
        pass

    tsv = not args['csv']

    new_root = h5_file.create_group(args['dataset'])

    if args['datasetname']:
        new_root.attrs.create('name', args['datasetname'])

    new_root.attrs.create('dfA', args['dfA'])
    new_root.attrs.create('dfX', args['dfX'])

    num_features = create_features(new_root, args['features'], tsv)
    num_markers = create_markers(new_root, args['markers'], tsv)
    shape_lod = create_lod(new_root, args['lod'], tsv)

    if num_features != shape_lod[0]:
        logging.debug('Error, rows ({}) in {} != rows ({}) in {}'.format(shape_lod[0], args['lod'], num_features, args['features']))
        sys.exit(-1)

    if num_markers != shape_lod[1]:
        logging.debug('Error, columns ({}) in {} != rows ({}) in {}'.format(shape_lod[1], args['lod'], num_markers, args['markers']))
        sys.exit(-1)

    if args['strains'] and args['coef']:
        num_strains = create_strains(new_root, args['strains'], tsv)
        shape_coef = create_coef(new_root, args['coef'], tsv)
    elif args['strains'] and not args['coef']:
        logging.debug('--strains specified, must have --coef')
        sys.exit()
    elif not args['strains'] and args['coef']:
        logging.debug('--coef specified, must have --strains')
        sys.exit()

    if args['samples']:
        num_samples = create_samples(new_root, args['samples'], tsv)

    if args['factors'] and args['phenotypes']:
        num_factors = create_factors(new_root, args['factors'], tsv)
        shape_phenotypes = create_phenotypes(new_root, args['phenotypes'], tsv)
    elif args['factors'] and not args['phenotypes']:
        logging.debug('--factors specified, must have --phenotypes')
        sys.exit()
    elif not args['factors'] and args['phenotypes']:
        logging.debug('--phenotypes specified, must have --factors')
        sys.exit()

    if args['genotypes'] and args['samples']:
        num_genotypes = create_genotypes(new_root, args['genotypes'], tsv)
    elif args['genotypes'] and not args['samples']:
        logging.debug('--genotypes specified, must have --samples')
        sys.exit()

    if args['expression'] and args['samples']:
        shape_expression = create_expression(new_root, args['expression'], tsv)
    elif args['expression'] and not args['samples']:
        logging.debug('--expression specified, must have --samples')
        sys.exit()

    h5_file.flush()
    h5_file.close()

    logging.debug("DONE: {} created".format(ffilename(args['out'])))

