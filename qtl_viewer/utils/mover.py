import h5py
import numpy as np
import os
import rpy2.robjects as ro
import random

h5_file_name = 'data/new_data2.h5'
old_h5_file_name = 'data/new_data.h5'
r_file_name = '/Users/mvincent/work/eqtl-viewer/data/R/expression.adipose.Rdata'

vstr = h5py.special_dtype(vlen=bytes)

print 'removing old file'

try:
    os.remove(h5_file_name)
except:
    pass

h5_file = h5py.File(h5_file_name, "w")

old_h5_file = h5py.File(old_h5_file_name, "r")

print 'copying file'


new_root = h5_file.create_group("example")
new_root.attrs.create('name', 'My Example')
old_root = old_h5_file['example']

old_root.copy('markers', new_root)
old_root.copy('features', new_root)
old_root.copy('lod', new_root)

print 'file copied, loading R data'

R_DATA = ro.r['load'](r_file_name)
EXP_DATA = ro.r[R_DATA[0]]

print 'r data loaded, creating samples'

# samples
dt = np.dtype([('sample_id', vstr), ('name', vstr), ('description', vstr)])
samples = []
for i in EXP_DATA.rownames:
    samples.append((i, i, None))
np_samples = np.array(samples, dtype=dt)
new_root = h5_file['example']
new_root.create_dataset('samples', data=np_samples)

print 'samples done, creating phenotypes'

# phenotypes
phenotypes = new_root.create_group("phenotypes")
dt = np.dtype([('factor_id', vstr), ('name', vstr), ('description', vstr)])
factors = []
factors.append(('sex', 'Sex', None))
factors.append(('diet', 'Diet', None))
np_factors = np.array(factors, dtype=dt)
phenotypes.create_dataset('factors', data=np_factors)

pheno_data = []
for x in xrange(0, len(samples)):
    print x
    row = []
    row.append('M' if random.randint(0,1) else 'F')
    row.append('LF' if random.randint(0,1) else 'HF')
    pheno_data.append(row)
np_pheno_data = np.array(pheno_data, ndmin=2)
phenotypes.create_dataset('phenotypes', data=np_pheno_data)

print np_pheno_data

# genotypes
print 'phenotypes done, creating genotypes'
genotypes = new_root.create_group("genotypes")

g = []
for x in xrange(0, len(new_root['markers'])):
    temp = []
    for y in xrange(0, len(samples)):
        gt = '{}{}'.format('B' if random.randint(0,1) else 'R', 'B' if random.randint(0,1) else 'R')
        temp.append(gt)
    g.append(temp)
genotypes.create_dataset('genotypes', data=g)


# expression
print 'genotypes done, creating expression'
expression = new_root.create_group("expression")
expression.create_dataset('expression', data=np.random.uniform(low=-6.1, high=5.8, size=(len(new_root['features']),len(samples),)))


print 'expression done, closing files'

old_h5_file.close()
h5_file.flush()
h5_file.close()

print 'done!'