# -*- coding: utf-8 -*-

"""
Copyright (c) 2015 Matthew Vincent, The Churchill Lab

This software was developed by Gary Churchill's Lab.
(see http://research.jax.org/faculty/churchill)

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this software. If not, see <http://www.gnu.org/licenses/>.
"""

from collections import OrderedDict
import logging
import math
import os
import sys

import h5py
import numpy as np
import scipy.stats as stats

log = logging.getLogger()
log.addHandler(logging.StreamHandler())

HDF5_FILENAME = None
HDF5 = None

FEATURES = {}
MARKERS = {}

LOG10 = math.log(10)

vstr = h5py.special_dtype(vlen=bytes)


class DatasetError(Exception):
    """Error accessing the dataset in the HDF5 file
    Args:
        msg (str): Human readable string describing the exception.
        code (Optional[int]): Error code.

    Attributes:
        msg (str): Human readable string describing the exception.
        code (int): Exception error code.

    """
    def __init__(self, msg, code=None):
        self.msg = msg
        self.code = code


def open_hdf5(filename=None):
    f = filename if filename else HDF5_FILENAME
    return h5py.File(f, 'r')


def get_hdf5(filename=None):
    global HDF5
    if HDF5 is None:
        f = filename if filename else HDF5_FILENAME
        HDF5 = open_hdf5(f)
    return HDF5


def get_dataset(dataset):
    """Get the dataset in the HDF5 file.

    Args:
        dataset: the dataset id

    Returns:
        the dataset

    Raises:
        DatasetError: If the dataset doesn't exist.

    """
    h5 = get_hdf5()

    if dataset not in h5:
        raise DatasetError("/{} doesn't exist".format(dataset))

    return h5[dataset]



def init():
    """Build a dictionary to look up the feature index by feature_id
    and marker index by marker_id.

    Returns:
        None

    Raises:
        DatasetError: If features or markers are not in the file

    """
    h5 = get_hdf5()

    try:
        for dataset in h5.keys():
            feature_map = {}
            features = np.array(h5[dataset]['features'])

            for idx, feature in enumerate(features):
                feature_map[feature['feature_id']] = idx

            FEATURES[dataset] = feature_map
    except KeyError, ke:
        message = "Error in /{}/features".format(dataset)
        logging.error(message)
        raise DatasetError(message)

    try:
        for dataset in h5.keys():
            marker_map = {}
            markers = np.array(h5[dataset]['markers'])

            for idx, marker in enumerate(markers):
                marker_map[marker['marker_id']] = idx

            MARKERS[dataset] = marker_map
    except KeyError, ke:
        message = "Error in /{}/markers".format(dataset)
        logging.error(message)
        raise DatasetError(message)


def get_datasets():
    """Return the datasets in the HDF5 file.

    Returns:
        OrderedDict: with keys being the id of the dataset, and values being
        the name (or the id if 'name' is not listed as an attribute)

    """
    h5 = get_hdf5()
    datasets = OrderedDict()

    for k in h5.keys():
        dataset = h5[k]

        if 'name' in dataset.attrs:
            datasets[str(k)] = str(dataset.attrs['name'])
        else:
            datasets[str(k)] = str(k)

    return datasets


def generate_matrix_file(dataset, output_file=None, delimiter='\t'):
    """Generate the data file used to draw the matrix.

    Args:
        dataset (str): dataset id
        output_file (str): name of the output file
        delimiter (str): file delimiter. Defaults to '\t'

    Returns:
        None

    """
    try:
        h5 = get_hdf5()

        if output_file is None:
            out = sys.stdout
        else:
            out = open(output_file, 'w')

        root = h5[dataset]
        lods = root['lod']['lod']
        features = root['features']
        markers = root['markers']
        has_pvalues = has_pval(dataset)
        if has_pvalues:
            pvals = root['lod']['pval']

        header = delimiter.join(["feature_id", "feature_group_id",
                                 "feature_chrom", "feature_location",
                                 "feature_name", "feature_description",
                                 "marker_id", "marker_chrom",
                                 "marker_location", "marker_name",
                                 "marker_description", "lod_score", "p_value"])

        out.write(header)
        out.write('\n')

        for i, rec in enumerate(lods):
            marker_idx = rec.argmax()
            m = markers[marker_idx]
            f = features[i]
            p = ''
            if has_pvalues:
                p = pvals[i][marker_idx]

            line = [f['feature_id'], f['group_id'],
                    f['chrom'], f['location'],
                    f['name'], f['description'],
                    m['marker_id'], m['chrom'],
                    m['location'], m['name'],
                    m['description'], rec[marker_idx], p]
            out.write(delimiter.join(map(str, line)))
            out.write('\n')

    except KeyError, ke:
        # dataset doesn't exist
        logging.error(ke.message)
        raise DatasetError("Path under /{} doesn't exist".format(dataset))


def generate_matrix_files():
    """Generate a matrix file for each dataset.

    The file will be in a sub directory named data and the file name will
    be <dataset>.tsv
    """
    h5 = get_hdf5()
    for dataset in h5.keys():
        print dataset
        filename = "data/{}.tsv".format(dataset)
        try:
            os.remove(filename)
        except:
            pass

        generate_matrix_file(dataset, filename, '\t')


def get_strains(dataset):
    """Get the strains in the dataset.

    Args:
        dataset (str): dataset id

    Returns:
        numpy.ndarray: All strains

    Raises:
        DatasetError: If the dataset doesn't exist or unable to get strains.

    """
    root = get_dataset(dataset)

    try:
        return root['coef/strains'][()]
    except KeyError, ke:
        logging.error(ke.message)
        raise DatasetError("/{}/coef/strains doesn't exist".format(dataset))


def get_samples(dataset):
    """Get the samples in the dataset.

    Args:
        dataset (str): dataset id

    Returns:
        numpy.ndarray: All samples

    Raises:
        DatasetError: If the dataset doesn't exist or unable to get factors.

    """
    root = get_dataset(dataset)

    try:
        return root['samples'][()]
    except KeyError, ke:
        logging.error(ke.message)
        raise DatasetError("/{}/samples doesn't exist".format(dataset))


def get_factors(dataset):
    """Get the factors in the dataset.

    Args:
        dataset (str): dataset id

    Returns:
        numpy.ndarray: All factors

    Raises:
        DatasetError: If the dataset doesn't exist or unable to get factors.

    """
    root = get_dataset(dataset)

    try:
        return root['phenotypes/factors'][()]
    except KeyError, ke:
        logging.error(ke.message)
        raise DatasetError("/{}/phenotypes/factors doesn't exist".format(dataset))


def get_factors_web(dataset):
    """Get the factors in the dataset.

    Args:
        dataset (str): dataset id

    Returns:
        numpy.ndarray: All factors

    Raises:
        DatasetError: If the dataset doesn't exist or unable to get factors.

    """
    root = get_dataset(dataset)

    try:

        factors = OrderedDict()

        for f in root['phenotypes/factors']:
            factors[f['factor_id']] = {'factor_id': f['factor_id'],
                                       'name': f['name'],
                                       'description': f['description'],
                                       'levels': set()}
    except KeyError, ke:
        logging.error(ke.message)
        raise DatasetError("/{}/phenotypes/factors doesn't exist".format(dataset))

    try:

        for p in root['phenotypes/phenotypes']:
            for i, f in enumerate(factors):
                factors[f]['levels'].add(p[i])

    except KeyError, ke:
        logging.error(ke.message)
        raise DatasetError("/{}/phenotypes/phenotypes doesn't exist".format(dataset))

    # check to see if we have genotypes

    if 'genotypes/genotypes' in root:
        try:
            genotypes = np.array([])

            for g in root['genotypes/genotypes']:
                genotypes = np.union1d(genotypes, g)

            factors['genotype'] = {'factor_id': 'genotype',
                                   'name': 'Genotype',
                                   'description': None,
                                   'levels': set(genotypes)}

        except KeyError, ke:
            logging.error(ke.message)
            raise DatasetError("/{}/genotypes/genotypes doesn't exist".format(dataset))

    return factors


def get_factors_phenotypes(dataset, marker_id=None):
    """

    Args:
        dataset:
        marker_id:

    Returns:

    """
    factors = get_factors_web(dataset)
    root = get_dataset(dataset)

    try:
        for k, v in factors.iteritems():
            factors[k]['values'] = []

        for p in root['phenotypes/phenotypes']:
            for i, k in enumerate(factors.keys()):
                if k != 'genotype':
                    factors[k]['values'].append(p[i])

    except KeyError, ke:
        logging.error(ke.message)
        raise DatasetError("/{}/phenotypes/phenotypes doesn't exist".format(dataset))

    try:
        if 'genotypes/genotypes' in root and marker_id is not None:
            idx = MARKERS[dataset][marker_id]
            factors['genotype']['values'] = list(root['genotypes/genotypes'][idx])
    except KeyError, ke:
        logging.error(ke.message)
        print ke
        raise DatasetError("/{}/genotypes/genotypes doesn't exist".format(dataset))


    return factors


def get_expression(dataset, identifier):
    h5 = get_hdf5()
    try:
        data = []
        idx = FEATURES[dataset][identifier]
        for i, l in enumerate(h5[dataset]['expression']['expression'][idx]):
            data.append(l)
        return data
    except Exception, e:
        print str(e)
    return None


def get_lod_data(dataset, identifier):
    """

    Args:
        dataset:
        identifier:

    Returns:

    """
    h5 = get_hdf5()
    data = [] # marker_id, chrom, location, lod

    # TODO: check that dataset exists

    idx = FEATURES[dataset][identifier]
    markers = h5[dataset]['markers'][()]

    for i, l in enumerate(h5[dataset]['lod']['lod'][idx]):
        data.append([markers[i]['marker_id'], markers[i]['chrom'], markers[i]['location'], l])

    return data


def has_samples(dataset):
    h5 = get_hdf5()
    return 'samples' in h5[dataset]


def has_genotypes(dataset):
    h5 = get_hdf5()
    return 'genotypes' in h5[dataset]


def has_pval(dataset):
    h5 = get_hdf5()
    return 'lod/pval' in h5[dataset]


def get_max_lod(dataset):
    #h5 = get_hdf5()
    #return 'lod/pval' in h5[dataset]
    raise Exception("Not yet implemented")


def get_genotypes(dataset, marker_id=None):
    h5 = get_hdf5()

    genotypes = np.array([])
    genotypes_data = h5[dataset]['genotypes']['genotypes'][()]

    for g in genotypes_data:
        genotypes = np.union1d(genotypes, g)

    if marker_id:
        idx = MARKERS[dataset][marker_id]
        genotypes_data = genotypes_data[idx]

    return {'genotypes':list(genotypes), 'values':genotypes_data}


def get_effect_data(dataset, identifier):
    """

    Args:
        dataset:
        identifier:

    Returns:

    """
    h5 = get_hdf5()
    #data = [] # marker_id, chrom, location, lod

    # TODO: check that dataset exists

    idx = FEATURES[dataset][identifier]
    strains = get_strains(dataset)

    top = ['marker_id']
    top.extend(zip(*list(h5[dataset]['markers'][()]))[0])

    data = []
    data.append(top)

    for i, c in enumerate(h5[dataset]['coef']['coef'][idx]):
        line = [strains[i]['strain_id']]
        line.extend(c)
        data.append(line)
        #print line
        #print '\n\n\n'

    # data is now list of list of (founder, coef, coef, coef, ..., coef * number of markers)

    lines = []

    # transpose
    for d in zip(*data):
        lines.append("\t".join(map(str, d)))

    # lines is now list of (marker, founder val, founder val, ..., founder val * number of founders) * number markers

    return '\n'.join(lines)


def pvalue2lod(pvalue, df):
    lod = stats.chi2.ppf(1-pvalue, df) / 2 / math.log(10)
    return lod


def pvalue2lodtable(dataset, x_min, x_max):
    root = get_dataset(dataset)

    vals = {}

    if 'dfA' not in root.attrs:
        raise DatasetError('dfA not found in dataset /{}'.format(dataset))

    if 'dfX' not in root.attrs:
        raise DatasetError('dfX not found in dataset /{}'.format(dataset))

    dfA = root.attrs['dfA']
    dfX = root.attrs['dfX']

    vals['dfA'] = {'df': dfA}
    vals['dfX'] = {'df': dfX}

    for x in xrange(x_min, x_max+1):
        p = 10**(-1 * x)
        vals['dfA'][x] = pvalue2lod(p, dfA)
        vals['dfX'][x] = pvalue2lod(p, dfX)

    return vals

if __name__ == '__main__':
    HDF5_FILENAME = 'data/btbr_compression.h5'
    init()
    #h5f = h5py.File('data/new_data.h5')
    #print get_datasets()

    #strains = get_strains('example')
    #for s in strains:
    #    print s['strain_id'], s['name']

    #print get_effect_data('example', 'ENSMUSG00000004446')

    #print get_genotypes('example', '1_2118117')
    #print get_factors_phenotypes('example', '1_2118117')


    #generate_matrix_files()

    #print has_pval('exampleB')
    #generate_matrix_file('exampleA')

    #print pvalue2lod(0.0001, 7), pvalue2lod(0.0001, 14)
    #print pvalue2lod(0.000001, 7), pvalue2lod(0.000001, 14)

    print str(pvalue2lodtable('adipose', 1, 16))







