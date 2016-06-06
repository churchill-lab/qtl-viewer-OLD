# -*- coding: utf-8 -*-

from collections import namedtuple, OrderedDict
import data_utils
import h5py
import numpy as np
import os
import sys

Item = namedtuple('Item', ['path', 'datatype', 'shape'])



items = OrderedDict()

def check_type(path, field, dtype, num=True):
    #print dtype[field].kind
    if num:
        t = np.typecodes['AllFloat'] + np.typecodes['AllInteger']
        if field in dtype.names and (dtype[field].kind in t):
            print 'SUCCESS\t{} found field {} as NUMERIC'.format(path, field)
        else:
            print 'ERROR\t{} found field {} as NON NUMERIC'.format(path, field)
    else:
        # TODO: O is python object, not really a string
        t = 'OaSc'
        if field in dtype.names and (dtype[field].kind in t):
            print 'SUCCESS\t{} found field {} as STRING'.format(path, field)
        else:
            print 'ERROR\t{} found field {} as NON STRING'.format(path, field)



def create_item(name, obj):
    global items
    n = '/{}'.format(name)
    items[n] = Item(n, type(obj), obj.shape if isinstance(obj, h5py.Dataset) else ())
    #for key, val in obj.attrs.iteritems():
    #    print "    %s: %s" % (key, val)



def main(filename, dump=False):
    data_utils.HDF5_FILENAME = os.path.abspath(filename)
    data_utils.init()

    root = data_utils.HDF5
    root.visititems(create_item)

    #global items
    #for i in items:
    #    print i, items[i]

    print 'Analysis of: {}'.format(root.filename)

    print 'Number of data sets under /: {}'.format(len(root.keys()))

    for dataset_name, dataset in root.iteritems():
        print 'Analyzing /{}'.format(dataset_name)

        attrs = dataset.attrs

        if 'dfA' in attrs:
            print 'SUCCESS\t/{} has attribute dfA: {}'.format(dataset_name, attrs['dfA'])
        else:
            print 'ERROR\t/{} has no dfA attribute/{}'.format(dataset_name)

        if 'dfX' in attrs:
            print 'SUCCESS\t/{} has attribute dfX: {}'.format(dataset_name, attrs['dfX'])
        else:
            print 'ERROR\t/{} has no dfX attribute/{}'.format(dataset_name)

        keys = dataset.keys()

        print 'Items: {}'.format(len(keys))

        # check features
        num_features = 0
        path = '/{}/features'.format(dataset_name)
        print '\nChecking {}\n'.format(path)
        if 'features' in dataset:
            features = dataset['features']
            print 'SUCCESS {} found'.format(path)

            if isinstance(features, h5py.Dataset):
                print 'SUCCESS {} is a Dataset'.format(path)

                # check types
                fields = set(features.dtype.names) #set(dtype for dtype in features.dtype.fields)
                req_fields = set(['feature_id', 'group_id', 'chrom', 'location', 'name', 'description'])

                diff = req_fields.difference(fields)
                if len(diff) != 0:
                    print 'ERROR\t{} are missing required fields: {}'.format(path, ','.join(diff))

                    diff = fields.difference(req_fields)
                    if len(diff) != 0:
                        print 'ERROR\t{} have extra fields: {}'.format(path, ','.join(diff))
                else:
                    print 'SUCCESS\t{} has all required fields'.format(path)

                check_type(path, 'feature_id', features.dtype, False)
                check_type(path, 'group_id', features.dtype, False)
                check_type(path, 'chrom', features.dtype, False)
                check_type(path, 'location', features.dtype, True)
                check_type(path, 'name', features.dtype, False)
                check_type(path, 'description', features.dtype, False)

                shape = features.shape
                if len(shape) != 1:
                    print 'ERROR\t{} is multidimensional: {}'.format(path, shape)
                else:
                    if shape[0] == 0:
                        print 'ERROR\t{} has no data'.format(path)
                    else:
                        print 'SUCCESS\t{} has {:,} items'.format(path, shape[0])
                        num_features = shape[0]

                        if dump:
                            np.savetxt("{}_features.txt".format(dataset_name), features.value, fmt="%s", delimiter="\t")

            else:
                print 'ERROR\t{} not a Dataset'.format(path)

        else:
            print 'ERROR\t{} not found'.format(path)

        # check markers
        num_markers = 0
        path = '/{}/markers'.format(dataset_name)
        print '\nChecking {}\n'.format(path)
        if 'markers' in dataset:
            markers = dataset['markers']
            print 'SUCCESS {} found'.format(path)

            if isinstance(markers, h5py.Dataset):
                print 'SUCCESS {} is a Dataset'.format(path)

                # check types
                fields = set(markers.dtype.names) #set(dtype for dtype in features.dtype.fields)
                req_fields = set(['marker_id', 'chrom', 'location'])

                diff = req_fields.difference(fields)
                if len(diff) != 0:
                    print 'ERROR\t{} are missing required fields: {}'.format(path, ','.join(diff))

                    diff = fields.difference(req_fields)
                    if len(diff) != 0:
                        print 'ERROR\t{} have extra fields: {}'.format(path, ','.join(diff))
                else:
                    print 'SUCCESS\t{} has all required fields'.format(path)

                check_type(path, 'marker_id', markers.dtype, False)
                check_type(path, 'chrom', markers.dtype, False)
                check_type(path, 'location', markers.dtype, True)

                shape = markers.shape
                if len(shape) != 1:
                    print 'ERROR\t{} is multidimensional: {}'.format(path, shape)
                else:
                    if shape[0] == 0:
                        print 'ERROR\t{} has no data'.format(path)
                    else:
                        print 'SUCCESS\t{} has {:,} items'.format(path, shape[0])
                        num_markers = shape[0]

                        if dump:
                            np.savetxt("{}_markers.txt".format(dataset_name), markers.value, fmt="%s", delimiter="\t")

            else:
                print 'ERROR\t{} not a Dataset'.format(path)

        else:
            print 'ERROR\t{} not found'.format(path)

        # check lod
        path = '/{}/lod'.format(dataset_name)
        print '\nChecking {}\n'.format(path)

        if 'lod' in dataset:
            print 'SUCCESS\t{} found'.format(path)
            lod = dataset['lod']

            if isinstance(lod, h5py.Group):
                print 'SUCCESS\t{} is a Group'.format(path)
            else:
                print 'ERROR\t{} is not a Group'.format(path)
        else:
            print 'ERROR\t{} not found'.format(path)

        # check lod/lod

        path = '/{}/lod/lod'.format(dataset_name)
        print '\nChecking {}\n'.format(path)
        if 'lod/lod' in dataset:
            lod = dataset['lod/lod']
            print 'SUCCESS\t{} found'.format(path)

            if isinstance(lod, h5py.Dataset):
                print 'SUCCESS {} is a Dataset'.format(path)
                shape = lod.shape
                if len(shape) != 2:
                    print 'ERROR\t{} does not have correct shape: {}'.format(path, shape)
                else:
                    if shape[0] == num_features:
                        print 'SUCCESS\t{} #rows ({:,}) equals #features ({:,})'.format(path, shape[0], num_features)
                    else:
                        print 'ERROR\t{} #rows ({:,}) does not equal #features ({:,})'.format(path, shape[0], num_features)

                    if shape[1] == num_markers:
                        print 'SUCCESS\t{} #columns ({:,}) equals #markers ({:,})'.format(path, shape[1], num_markers)
                    else:
                        print 'ERROR\t{} #columns ({:,}) does not equal #markers ({:,})'.format(path, shape[1], num_markers)

                    if dump and shape[0] == num_features and shape[1] == num_markers:
                        np.savetxt("{}_lod.txt".format(dataset_name), lod.value, fmt="%f", delimiter="\t")

            else:
                print 'ERROR\t{} not a Dataset'.format(path)

        else:
            print 'ERROR\t{} not found'.format(path)

        # coef
        num_strains = 0
        coef_found = False
        path = '/{}/coef'.format(dataset_name)
        print '\nChecking {}\n'.format(path)

        if 'coef' in dataset:
            print 'SUCCESS\t{} found'.format(path)
            coef = dataset['coef']

            if isinstance(coef, h5py.Group):
                print 'SUCCESS\t{} is a Group'.format(path)
                coef_found = True
            else:
                print 'ERROR\t{} is not a Group'.format(path)
        else:
            print 'WARNING\t{} not found'.format(path)

        if coef_found:
            path = '/{}/coef/strains'.format(dataset_name)
            print '\nChecking {}\n'.format(path)
            if 'coef/strains' in dataset:
                strains = dataset['coef/strains']
                print 'SUCCESS\t{} found'.format(path)

                if isinstance(strains, h5py.Dataset):
                    print 'SUCCESS {} is a Dataset'.format(path)

                    # check types
                    fields = set(strains.dtype.names)
                    req_fields = set(['strain_id', 'name', 'description'])

                    diff = req_fields.difference(fields)
                    if len(diff) != 0:
                        print 'ERROR\t{} are missing required fields: {}'.format(path, ','.join(diff))

                        diff = fields.difference(req_fields)
                        if len(diff) != 0:
                            print 'ERROR\t{} have extra fields: {}'.format(path, ','.join(diff))
                    else:
                        print 'SUCCESS\t{} has all required fields'.format(path)

                    check_type(path, 'strain_id', strains.dtype, False)
                    check_type(path, 'name', strains.dtype, False)
                    check_type(path, 'description', strains.dtype, False)

                    shape = strains.shape
                    if len(shape) != 1:
                        print 'ERROR\t{} is multidimensional: {}'.format(path, shape)
                    else:
                        if shape[0] == 0:
                            print 'ERROR\t{} has no data'.format(path)
                        else:
                            print 'SUCCESS\t{} has {:,} items'.format(path, shape[0])
                            num_strains = shape[0]

                            if dump:
                                np.savetxt("{}_strains.txt".format(dataset_name), strains.value, fmt=["%s","%s","%s"], delimiter="\t")

                else:
                    print 'ERROR\t{} not a Dataset'.format(path)

            else:
                print 'ERROR\t{} not found'.format(path)


            path = '/{}/coef/coef'.format(dataset_name)
            print '\nChecking {}\n'.format(path)
            if 'coef/coef' in dataset:
                coef = dataset['coef/coef']
                print 'SUCCESS\t{} found'.format(path)

                if isinstance(coef, h5py.Dataset):
                    print 'SUCCESS {} is a Dataset'.format(path)
                    shape = coef.shape
                    if len(shape) != 3:
                        print 'ERROR\t{} does not have correct shape: {}'.format(path, shape)
                    else:
                        if shape[0] == num_features:
                            print 'SUCCESS\t{} 1st dimension ({:,}) equals #features ({:,})'.format(path, shape[0], num_features)
                        else:
                            print 'ERROR\t{} 1st dimension ({:,}) does not equal #features ({:,})'.format(path, shape[0], num_features)

                        if shape[1] == num_strains:
                            print 'SUCCESS\t{} 2nd dimension ({:,}) equals #strains ({:,})'.format(path, shape[1], num_strains)
                        else:
                            print 'ERROR\t{} 2nd dimension ({:,}) does not equal #strains ({:,})'.format(path, shape[1], num_strains)

                        if shape[2] == num_markers:
                            print 'SUCCESS\t{} 3rd dimension ({:,}) equals #markers ({:,})'.format(path, shape[2], num_markers)
                        else:
                            print 'ERROR\t{} 3rd dimension ({:,}) does not equal #markers ({:,})'.format(path, shape[2], num_markers)

                    if dump and shape[0] == num_features and shape[1] == num_strains and shape[2] == num_markers:
                        data = np.rollaxis(coef.value, 1)
                        for x in xrange(num_strains):
                            strain_data = data[x,]
                            np.savetxt("{}_coef_{}.txt".format(dataset_name, strains[x][0]), strain_data, fmt="%s", delimiter="\t")

                else:
                    print 'ERROR\t{} not a Dataset'.format(path)

            else:
                print 'ERROR\t{} not found'.format(path)


        # samples
        num_samples = 0
        path = '/{}/samples'.format(dataset_name)
        print '\nChecking {}\n'.format(path)
        if 'samples' in dataset:
            samples = dataset['samples']
            print 'SUCCESS\t{} found'.format(path)

            if isinstance(samples, h5py.Dataset):
                print 'SUCCESS {} is a Dataset'.format(path)

                # check types
                fields = set(samples.dtype.names)
                req_fields = set(['sample_id', 'name', 'description'])

                diff = req_fields.difference(fields)
                if len(diff) != 0:
                    print 'ERROR\t{} are missing required fields: {}'.format(path, ','.join(diff))

                    diff = fields.difference(req_fields)
                    if len(diff) != 0:
                        print 'ERROR\t{} have extra fields: {}'.format(path, ','.join(diff))
                else:
                    print 'SUCCESS\t{} has all required fields'.format(path)

                check_type(path, 'sample_id', samples.dtype, False)
                check_type(path, 'name', samples.dtype, False)
                check_type(path, 'description', samples.dtype, False)

                shape = samples.shape
                if len(shape) != 1:
                    print 'ERROR\t{} is multidimensional: {}'.format(path, shape)
                else:
                    if shape[0] == 0:
                        print 'ERROR\t{} has no data'.format(path)
                    else:
                        print 'SUCCESS\t{} has {:,} items'.format(path, shape[0])
                        num_samples = shape[0]

                        if dump:
                            np.savetxt("{}_samples.txt".format(dataset_name), samples.value, fmt=["%s","%s","%s"], delimiter="\t")

            else:
                print 'ERROR\t{} not a Dataset'.format(path)

        else:
            print 'WARNING\t{} not found'.format(path)


        # check expression
        path = '/{}/expression'.format(dataset_name)
        print '\nChecking {}\n'.format(path)

        if 'expression' in dataset:
            print 'SUCCESS\t{} found'.format(path)
            expression = dataset['expression']

            if isinstance(expression, h5py.Group):
                print 'SUCCESS\t{} is a Group'.format(path)
            else:
                print 'ERROR\t{} is not a Group'.format(path)
        else:
            print 'WARNING\t{} not found'.format(path)

        # check expression/expression

        path = '/{}/expression/expression'.format(dataset_name)
        print '\nChecking {}\n'.format(path)
        if 'expression/expression' in dataset:
            expression = dataset['expression/expression']
            print 'SUCCESS\t{} found'.format(path)

            if isinstance(expression, h5py.Dataset):
                print 'SUCCESS {} is a Dataset'.format(path)
                shape = expression.shape
                if len(shape) != 2:
                    print 'ERROR\t{} does not have correct shape: {}'.format(path, shape)
                else:
                    if shape[0] == num_features:
                        print 'SUCCESS\t{} #rows ({:,}) equals #features ({:,})'.format(path, shape[0], num_features)
                    else:
                        print 'ERROR\t{} #rows ({:,}) does not equal #features ({:,})'.format(path, shape[0], num_features)

                    if shape[1] == num_samples:
                        print 'SUCCESS\t{} #columns ({:,}) equals #samples ({:,})'.format(path, shape[1], num_samples)
                    else:
                        print 'ERROR\t{} #columns ({:,}) does not equal #samples ({:,})'.format(path, shape[1], num_samples)

                    if dump and shape[0] == num_features and shape[1] == num_samples:
                        np.savetxt("{}_expression.txt".format(dataset_name), expression.value, fmt="%f", delimiter="\t")

            else:
                print 'ERROR\t{} not a Dataset'.format(path)

        else:
            print 'WARNING\t{} not found'.format(path)


        # check genotypes
        path = '/{}/genotypes'.format(dataset_name)
        print '\nChecking {}\n'.format(path)

        if 'genotypes' in dataset:
            print 'SUCCESS\t{} found'.format(path)
            genotypes = dataset['genotypes']

            if isinstance(genotypes, h5py.Group):
                print 'SUCCESS\t{} is a Group'.format(path)
            else:
                print 'ERROR\t{} is not a Group'.format(path)
        else:
            print 'WARNING\t{} not found'.format(path)

        # check genotypes/genotypes

        path = '/{}/genotypes/genotypes'.format(dataset_name)
        print '\nChecking {}\n'.format(path)
        if 'genotypes/genotypes' in dataset:
            genotypes = dataset['genotypes/genotypes']
            print 'SUCCESS\t{} found'.format(path)

            if isinstance(genotypes, h5py.Dataset):
                print 'SUCCESS {} is a Dataset'.format(path)
                shape = genotypes.shape
                if len(shape) != 2:
                    print 'ERROR\t{} does not have correct shape: {}'.format(path, shape)
                else:
                    if shape[0] == num_markers:
                        print 'SUCCESS\t{} #rows ({:,}) equals #markers ({:,})'.format(path, shape[0], num_markers)
                    else:
                        print 'ERROR\t{} #rows ({:,}) does not equal #markers ({:,})'.format(path, shape[0], num_markers)

                    if shape[1] == num_samples:
                        print 'SUCCESS\t{} #columns ({:,}) equals #samples ({:,})'.format(path, shape[1], num_samples)
                    else:
                        print 'ERROR\t{} #columns ({:,}) does not equal #samples ({:,})'.format(path, shape[1], num_samples)

                    if dump and shape[0] == num_markers and shape[1] == num_samples:
                        np.savetxt("{}_genotypes.txt".format(dataset_name), genotypes.value, fmt="%s", delimiter="\t")

            else:
                print 'ERROR\t{} not a Dataset'.format(path)

        else:
            print 'WARNING\t{} not found'.format(path)

        # phenotypes
        num_factors = 0
        pheno_found = False
        path = '/{}/phenotypes'.format(dataset_name)
        print '\nChecking {}\n'.format(path)

        if 'phenotypes' in dataset:
            print 'SUCCESS\t{} found'.format(path)
            phenotypes = dataset['phenotypes']

            if isinstance(phenotypes, h5py.Group):
                print 'SUCCESS\t{} is a Group'.format(path)
                pheno_found = True
            else:
                print 'ERROR\t{} is not a Group'.format(path)
        else:
            print 'WARNING\t{} not found'.format(path)

        if pheno_found:
            path = '/{}/phenotypes/factors'.format(dataset_name)
            print '\nChecking {}\n'.format(path)
            if 'phenotypes/factors' in dataset:
                factors = dataset['phenotypes/factors']
                print 'SUCCESS\t{} found'.format(path)

                if isinstance(factors, h5py.Dataset):
                    print 'SUCCESS {} is a Dataset'.format(path)

                    # check types
                    fields = set(factors.dtype.names)
                    req_fields = set(['factor_id', 'name', 'description'])

                    diff = req_fields.difference(fields)
                    if len(diff) != 0:
                        print 'ERROR\t{} are missing required fields: {}'.format(path, ','.join(diff))

                        diff = fields.difference(req_fields)
                        if len(diff) != 0:
                            print 'ERROR\t{} have extra fields: {}'.format(path, ','.join(diff))
                    else:
                        print 'SUCCESS\t{} has all required fields'.format(path)

                    check_type(path, 'factor_id', factors.dtype, False)
                    check_type(path, 'name', factors.dtype, False)
                    check_type(path, 'description', factors.dtype, False)

                    shape = factors.shape
                    if len(shape) != 1:
                        print 'ERROR\t{} is multidimensional: {}'.format(path, shape)
                    else:
                        if shape[0] == 0:
                            print 'ERROR\t{} has no data'.format(path)
                        else:
                            print 'SUCCESS\t{} has {:,} items'.format(path, shape[0])
                            num_factors = shape[0]

                        if dump:
                            np.savetxt("{}_factors.txt".format(dataset_name), factors.value, fmt=["%s","%s","%s"], delimiter="\t")
                else:
                    print 'ERROR\t{} not a Dataset'.format(path)

            else:
                print 'ERROR\t{} not found'.format(path)


            path = '/{}/phenotypes/phenotypes'.format(dataset_name)
            print '\nChecking {}\n'.format(path)
            if 'phenotypes/phenotypes' in dataset:
                phenotypes = dataset['phenotypes/phenotypes']
                print 'SUCCESS\t{} found'.format(path)

                if isinstance(phenotypes, h5py.Dataset):
                    print 'SUCCESS {} is a Dataset'.format(path)
                    shape = phenotypes.shape
                    if len(shape) != 2:
                        print 'ERROR\t{} does not have correct shape: {}'.format(path, shape)
                    else:
                        if shape[0] == num_samples:
                            print 'SUCCESS\t{} #rows ({:,}) equals #samples ({:,})'.format(path, shape[0], num_samples)
                        else:
                            print 'ERROR\t{} #rows ({:,}) does not equal #samples ({:,})'.format(path, shape[0], num_samples)

                        if shape[1] == num_factors:
                            print 'SUCCESS\t{} #columns ({:,}) equals #factors ({:,})'.format(path, shape[1], num_factors)
                        else:
                            print 'ERROR\t{} #columns ({:,}) does not equal #factors ({:,})'.format(path, shape[1], num_factors)

                    if dump and shape[0] == num_samples and shape[1] == num_factors:
                        np.savetxt("{}_phenotypes.txt".format(dataset_name), phenotypes.value, fmt="%s", delimiter="\t")
                else:
                    print 'ERROR\t{} not a Dataset'.format(path)

            else:
                print 'ERROR\t{} not found'.format(path)


if __name__ == '__main__':
    main(sys.argv[1])
