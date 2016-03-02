import h5py
import rpy2.robjects as ro
import numpy as np
import os


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
    h5_root.create_dataset('features', data=np_features, compression="lzf")

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
    h5_root.create_dataset('markers', data=np_markers, compression="lzf")

    #lod_scores = np.transpose(R_LOD_DATA[2:])
    lod_scores = np.array(R_LOD_DATA[2:])
    lod_root = h5_root.create_group('lod')
    lod_root.create_dataset('lod', data=lod_scores, compression="lzf")



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
    h5_root.create_dataset('samples', data=np_samples, compression="lzf")

    expression = R_EXPRESSION_DATA
    np_expression = np.array(expression)
    expression_root = h5_root.create_group('expression')
    expression_root.create_dataset('expression', data=np_expression, compression="lzf")

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
    phenotypes_root.create_dataset('phenotypes', data=np_phenotypes, compression="lzf")

    np_genotypes = np.transpose(genotypes)
    genotypes_root = h5_root.create_group('genotypes')
    genotypes_root.create_dataset('genotypes', data=np_genotypes, compression="lzf")


#
# factors
#
def create_factors(h5_root):
    factors = []
    factors.append(('sex', 'Sex', None))
    np_factors = np.array(factors, dtype=dt_factors)

    h5_root.create_dataset('phenotypes/factors', data=np_factors, compression="lzf")


if __name__ == '__main__':
    datasets = [('adipose', 'Adipose'),
                ('gastroc', 'Gastrocnemius'),
                ('hypo', 'Hypothalamus'),
                ('islet', 'Islet'),
                ('kidney', 'Kidney'),
                ('liver', 'Liver')]

    h5_file_name = 'data/btbr_compression.h5'

    try:
        os.remove(h5_file_name)
    except:
        pass

    h5_file = h5py.File(h5_file_name, "w")

    for d in datasets:
        print d
        new_root = h5_file.create_group(d[0])
        new_root.attrs.create('name', d[1])
        new_root.attrs.create('dfA', 7)
        new_root.attrs.create('dfX', 14)

        create_features(new_root)
        create_markers_lods(new_root, d[0])
        create_expression_samples(new_root, d[0])
        create_genotype_phenotype(new_root, d[0])
        create_factors(new_root)

    h5_file.flush()
    h5_file.close()

