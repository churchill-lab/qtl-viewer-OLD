import h5py
import rpy2.robjects as ro
import numpy as np

import argparse
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


#
# features
#

def create_features(h5_root):
    R_ANNOTATION_FILE='data/annotations.Rdata'
    R_ANNOTATION_META=ro.r['load'](R_ANNOTATION_FILE)
    R_ANNOTATION_DATA=ro.r['annotations']

    features = []
    ensembl_ids=R_ANNOTATION_DATA[3]
    chroms=R_ANNOTATION_DATA[1]
    locations=R_ANNOTATION_DATA[2]
    names=R_ANNOTATION_DATA[4]

    for i in xrange(len(ensembl_ids)):
        features.append((ensembl_ids.levels[ensembl_ids[i]-1], ensembl_ids.levels[ensembl_ids[i]-1], chroms.levels[chroms[i]-1], locations[i], names.levels[names[i]-1], None))

    np_features = np.array(features, dtype=dt_features)
    h5_root.create_dataset('features', data=np_features)

#
# markers and lod
#

def create_markers_lods(h5_root, dataset):
    R_LOD_FILE='data/f2g_scan_{}.Rdata'.format(dataset)
    R_LOD_META=ro.r['load'](R_LOD_FILE)
    R_LOD_DATA=ro.r['f2g.scan.{}'.format(dataset)]

    marker_ids=R_LOD_DATA.rownames
    chroms=R_LOD_DATA[0]
    locations=R_LOD_DATA[1]

    markers = []
    for i in xrange(len(marker_ids)):
        markers.append((marker_ids[i], chroms.levels[chroms[i]-1], locations[i], None, None))

    np_markers = np.array(markers, dtype=dt_markers)
    h5_root.create_dataset('markers', data=np_markers)

    #lod_scores = np.transpose(R_LOD_DATA[2:])
    lod_scores = np.array(R_LOD_DATA[2:])

    lod_root = h5_root.create_group('lod')
    lod_root.create_dataset('lod', data=lod_scores)



#
# expression and samples
#

def create_expression_samples(h5_root, dataset):
    R_EXPRESSION_FILE='data/expression.{}.Rdata'.format(dataset)
    R_EXPRESSION_META=ro.r['load'](R_EXPRESSION_FILE)
    R_EXPRESSION_DATA=ro.r[dataset]
    sample_names = R_EXPRESSION_DATA.rownames

    samples = []
    for r in R_EXPRESSION_DATA.rownames:
        samples.append((r, r, None))
    np_samples = np.array(samples, dtype=dt_samples)
    h5_root.create_dataset('samples', data=np_samples)

    expression = R_EXPRESSION_DATA
    np_expression = np.array(expression)
    expression_root = h5_root.create_group('expression')
    expression_root.create_dataset('expression', data=np_expression)

#
# phenotypes and genotypes
#

def create_genotype_phenotype(h5_root, dataset):
    R_GENOTYPE_FILE='data/GenotypeData.Rdata'
    R_GENOTYPE_META=ro.r['load'](R_GENOTYPE_FILE)
    R_GENOTYPE_DATA=ro.r['Genotype']

    phenotypes = []
    genotypes = []
    for r in R_GENOTYPE_DATA[3:]:
        phenotypes.append(r[0])
        g = []
        for d in r[1:]:
            g.append(d)
        genotypes.append(g)

    np_phenotypes = np.array(phenotypes)
    phenotypes_root = h5_root.create_group('phenotypes')
    phenotypes_root.create_dataset('phenotypes', data=np_phenotypes)

    np_genotypes = np.transpose(genotypes)
    genotypes_root = h5_root.create_group('genotypes')
    genotypes_root.create_dataset('genotypes', data=np_genotypes)


#
# factors
#
def create_factors(h5_root):
    factors = []
    factors.append(('sex', 'Sex', None))
    np_factors = np.array(factors, dtype=dt_factors)

    h5_root.create_dataset('phenotypes/factors', data=np_factors)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert txt files to hdf5')

    parser.add_argument('-o', '--out', help='Output HDF file', required=True)
    parser.add_argument('-d', '--dataset', help='ID of the dataset', required=True)
    parser.add_argument('-f', '--features', help='features file, think genes', required=True)
    parser.add_argument('-l', '--lod', help='lod score file, features (rows) by markers (columns)', required=True)
    parser.add_argument('-m', '--markers', help='markers file, think snps', required=True)

    parser.add_argument('--dfA', help='Degress of freedom autosomes', required=True)
    parser.add_argument('--dfX', help='Degrees of freedom X', required=True)

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

    args = vars(parser.parse_args())

    print str(args)

    # check the arguments










    sys.exit()
    '''

    h5_file = h5py.File(h5_file_name, "w")

    for d in datasets:
        print d
        new_root = h5_file.create_group(d[0])
        new_root.attrs.create('name', d[1])

        create_features(new_root)
        create_markers_lods(new_root, d[0])
        create_expression_samples(new_root, d[0])
        create_genotype_phenotype(new_root, d[0])
        create_factors(new_root)

    h5_file.flush()
    h5_file.close()



    '''




