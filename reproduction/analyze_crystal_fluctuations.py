# Analyze the atomic fluctuations in a protein crystal.
#
# This program reads a netCDF file produced by the script
# calculate_crystal_fluctuations.py.
#

#
# User-definable parameters
#

# The PDB code used in calculate_crystal_fluctuations.py.
pdb_code = '1IEE'

# The name of the PDB file to be read. This must be the same file
# as used for the calculation of the fluctuations!
pdb_file = '%s.pdb' % pdb_code

# The size of the crystal, given as the number of times the unit cell
# is repeated along each of the lattice vectors. The total number of
# copies of the unit cell is the product of these three numbers.
n1 = 1
n2 = 1
n3 = 1

#
# There should be no need to modify anything from here on, but feel free
# to modify and/or add the analysis code as you like.
#

from MMTK import *
from MMTK.PDB import PDBConfiguration
from Scientific.Geometry import Tensor
from Scientific import N, LA
from Scientific.IO.NetCDF import NetCDFFile
from Scientific.IO.ArrayIO import writeArray


# Read the PDB file and obtain the unit cell parameters.
conf = PDBConfiguration(pdb_file)
e1, e2, e3 = conf.basis
b1, b2, b3 = conf.reciprocal_basis

# Make a C-alpha model for each peptide chain in the asymmetric unit.
chains = Collection()
atom_map = {}
for c in conf.peptide_chains:
    chain = c.createPeptideChain(model='calpha')
    c.applyTo(chain, atom_map=atom_map)
    chains.addObject(chain)

# Calculate the average mass of a residue
av_mass = chains.mass()/chains.numberOfAtoms()

# Retrieve the B factors and ADPs from the PDB file. Rescale the B factors
# so that they are equal to the trace of the ADP tensor.
b_exp = []
adp_exp = []
for chain in chains:
    for residue in chain:
        b = atom_map[residue.peptide.C_alpha].properties['temperature_factor']
        b_exp.append(3.*b/(8.*N.pi**2))
        try:
            u = atom_map[residue.peptide.C_alpha].properties['u']
        except KeyError:
            u = None
        adp_exp.append(u)
print len([a for a in adp_exp if a is not None]), "out of", \
      chains.numberOfAtoms(), "atoms have ADPs"

if None not in adp_exp:
    # If all atoms have ADPs, don't use the B factors but the traces
    # of the ADPs instead
    f_exp = N.array([u.array for u in adp_exp])
    b_exp = f_exp[:,0,0]++f_exp[:,1,1]+f_exp[:,2,2]
else:
    b_exp = N.array(b_exp)

# Open the netCDF file containing the fluctuation data.
nc_filename = 'crystal_fluctuations_%s_%d_%d_%d.nc' % (pdb_code, n1, n2, n3)
results = NetCDFFile(nc_filename)
natoms = results.dimensions['atoms']
q = results.variables['q']
f_acoustic = results.variables['f_acoustic']
f_optical = results.variables['f_optical']
force_constant = results.variables['force_constant']
weight = results.variables['weight']

def plotDispersionRelation():
    dispersion = []
    for i in range(results.dimensions['wave']):
        qv = Vector(q[i])
        ev = force_constant[i][:3]
        if ev[0] < 0.:
            continue
        frequency = N.sqrt(ev/av_mass)/(2.*N.pi)
        for j in range(3):
            dispersion.append((qv.length()/(2.*N.pi), frequency[j]))
    dispersion.sort(lambda a, b: cmp(a[0], b[0]))
    writeArray(N.array(dispersion),
               "dispersion_%s_%d_%d_%d.plot" % (pdb_code, n1, n2 , n3))

def anisotropy(u):
    ev = LA.eigenvalues(u)
    return N.minimum.reduce(ev)/N.maximum.reduce(ev)

def adpCorrelation(u1, u2):
    if u1 is None or u2 is None:
        return None
    m1 = u1.inverse()
    m2 = u2.inverse()
    dm1 = N.multiply.reduce(m1.eigenvalues())
    dm2 = N.multiply.reduce(m2.eigenvalues())
    ds = N.multiply.reduce((m1+m2).eigenvalues())
    cc = N.sqrt(8.*N.sqrt(dm1*dm2)/ds)
    return cc

def sumFluctuations():
    f_a = 0.
    f_o = 0.
    for i in range(results.dimensions['wave']):
        qv = Vector(q[i])
        ql = qv.length()
        if ql == 0.:
            f_0 = (n1*n2*n3)*f_optical[i]
        f_a += f_acoustic[i]
        f_o += f_optical[i]
    return f_a, f_o, f_0

# Generate a plot of the dispersion relation.
plotDispersionRelation()

# Sum up the acoustic and optical contributions over all modes,
# and also retrieve the fluctuations f_0 of the single unit cell
# in periodic boundary conditions for comparison.
f_a, f_o, f_0 = sumFluctuations()
f = f_a + f_o
b_a = f_a[:,0,0]+f_a[:,1,1]+f_a[:,2,2]
b_o = f_o[:,0,0]+f_o[:,1,1]+f_o[:,2,2]
b_0 = f_0[:,0,0]+f_0[:,1,1]+f_0[:,2,2]
b = b_a + b_o

# Write plot files for the B factors.
writeArray(b_exp, 'bexp_%s.plot' % pdb_code)
writeArray(b_a[:len(b_exp)],
           'bacoustic_%s_%d_%d_%d.plot' % (pdb_code, n1, n2 ,n3))
writeArray(b_o[:len(b_exp)],
           'boptical_%s_%d_%d_%d.plot' % (pdb_code, n1, n2 ,n3))
writeArray(b[:len(b_exp)],
           'btotal_%s_%d_%d_%d.plot' % (pdb_code, n1, n2 ,n3))

# Write plot files for the anisotropies.
writeArray(N.array([anisotropy(u.array) for u in adp_exp]),
           'anisotropy_exp_%s_%d_%d_%d.plot' % (pdb_code, n1, n2 ,n3))
if n1*n2*n3 > 1:
    writeArray(N.array([anisotropy(u) for u in f_a[:len(b_exp)]]),
               'anisotropy_acoustic_%s_%d_%d_%d.plot' % (pdb_code, n1, n2 ,n3))
writeArray(N.array([anisotropy(u) for u in f_o[:len(b_exp)]]),
           'anisotropy_optical_%s_%d_%d_%d.plot' % (pdb_code, n1, n2 ,n3))
writeArray(N.array([anisotropy(u) for u in f[:len(b_exp)]]),
           'anisotropy_total_%s_%d_%d_%d.plot' % (pdb_code, n1, n2 ,n3))

if None not in adp_exp:
    # Write plot files for the ADP correlations between experimental
    # values and theoretical predictions.
    writeArray(N.array([adpCorrelation(Tensor(f[i]/b[i]),
                                       adp_exp[i]/adp_exp[i].trace())
                        for i in range(len(adp_exp))]),
               'adpc_%s_%d_%d_%d.plot' % (pdb_code, n1, n2, n3))

    writeArray(N.array([adpCorrelation(Tensor(f_0[i]/b_0[i]),
                                       adp_exp[i]/adp_exp[i].trace())
                        for i in range(len(adp_exp))]),
               'adpc_%s_1_1_1.plot' % pdb_code)


# Fit the experimental B factors to the model
#   B = s_acoustic*B_acoustis + s_optical*B_optical
b = b[:len(b_exp)]
s_a, s_o = LA.linear_least_squares(N.transpose([b_a, b_o]), b_exp)[0]
b_fitted = s_a*b_a+s_o*b_o
writeArray(b_fitted, 'bfitted_%s_%d_%d_%d.plot' % (pdb_code, n1, n2 ,n3))
