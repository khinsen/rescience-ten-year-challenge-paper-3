# Calculate the atomic fluctuations in a protein crystal.
#
# This program produces a netCDF file containing for each q vector
# the acoustic and optical mode contributions to the anisotropic
# atomic fluctuations. Use the program analyze_crystal_fluctuations.py
# for working with this file.
#

#
# User-definable parameters
#

# The PDB code is used to generate the names of the output files. It does not
# have to be a real PDB code; any text string is allowed.
pdb_code = '1IEE'

# The temperature (in Kelvin) at which the structure was obtained.
temperature = 120.

# A scale factor for the elastic network model.
scale_factor = 2.99

# The size of the crystal, given as the number of times the unit cell
# is repeated along each of the lattice vectors. The total number of
# copies of the unit cell is the product of these three numbers.
n1 = 1
n2 = 1
n3 = 1

# For convenience in batch mode, read paramters from the command line
import sys
_, pdb_code, n1, n2, n3 = sys.argv
n1 = int(n1)
n2 = int(n2)
n3 = int(n3)

pdb_file = '%s.pdb' % pdb_code

#
# There should be no need to modify anything from here on.
#


from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.ForceFields import CalphaForceField
from Scientific.IO.NetCDF import NetCDFFile
from Scientific import N, LA

# Read the PDB file and obtain the unit cell parameters.
conf = PDBConfiguration(pdb_file)
e1, e2, e3 = conf.basis
b1, b2, b3 = conf.reciprocal_basis

# Make a universe that corresponds to the unit cell
# with periodic boundary conditions.
universe = ParallelepipedicPeriodicUniverse((e1, e2, e3))

# Make a C-alpha model for each peptide chain in the asymmetric unit.
chains = Collection()
for c in conf.peptide_chains:
    chain = c.createPeptideChain(model='calpha')
    chains.addObject(chain)

# Make the required copies of the chains in the asymmetric unit
# by applying the crystallographic symmetry operations. Also
# build up atom lists for the asymmetric unit and for the unit
# cell in a well-defined order.
atoms = []
au_atoms = None
for so in conf.cs_transformations:
    image = deepcopy(chains)
    for atom in image.atomList():
        atom.setPosition(so(atom.position()))
    # The following lines move the center of mass of each copy
    # inside the unit cell. This is not necessary for the
    # calculation, but makes for nicer pictures.
    cm = image.centerOfMass()
    cm_fr = conf.to_fractional(cm)
    cm_fr = Vector(cm_fr[0] % 1., cm_fr[1] % 1., cm_fr[2] % 1.) \
            - Vector(0.5, 0.5, 0.5)
    cm = conf.from_fractional(cm_fr)
    image.translateTo(cm)
    # Add the copy of the asymmetric unit to the universe and its
    # atoms to the ordered atom list.
    universe.addObject(image)
    for chain in image:
        for residue in chain:
            atoms.append(residue.peptide.C_alpha)
    if au_atoms is None:
        au_atoms = copy(atoms)

natoms = len(atoms)
universe.foldCoordinatesIntoBox()
atom_indices = N.array([a.index for a in au_atoms])

# Calculate the Hessian matrix
cutoff = 2.5*Units.nm
universe.setForceField(CalphaForceField(cutoff, scale_factor))
e, fc = universe.energyAndForceConstants()

# Compile a list of all q vectors.
q_vectors = {}
for i1 in range(-n1/2+1, n1/2+1):
    for i2 in range(-n2/2+1, n2/2+1):
        for i3 in range(-n3/2+1, n3/2+1):
            q_vectors[(i1, i2, i3)] = 1
assert N.sum(q_vectors.values()) == n1*n2*n3

# Reduce the number of q vectors by applying inversion symmetry.
for i1 in range(-n1/2+1, n1/2+1):
    for i2 in range(-n2/2+1, n2/2+1):
        for i3 in range(-n3/2+1, n3/2+1):
            c1 = i1 if 2*i1 == n1 else -i1
            c2 = i2 if 2*i2 == n2 else -i2
            c3 = i3 if 2*i3 == n3 else -i3
            c = (c1, c2, c3)
            if c != (i1, i2, i3) and q_vectors.has_key(c):
                q_vectors[c] = 2
                del q_vectors[(i1, i2, i3)]
assert N.sum(q_vectors.values()) == n1*n2*n3

# Create the netCDF file for storing the results.
nc_filename = 'crystal_fluctuations_%s_%d_%d_%d.nc' % (pdb_code, n1, n2, n3)
results = NetCDFFile(nc_filename, 'w')
results.pdb_code = pdb_code
results.n1 = n1
results.n2 = n2
results.n3 = n3
results.createDimension('wave', len(q_vectors))
results.createDimension('mode', 3*universe.numberOfAtoms())
results.createDimension('atoms', len(au_atoms))
results.createDimension('xyz', 3)
results.createDimension('xyz1', 3)
results.createDimension('xyz2', 3)
results.createVariable('index', N.Int32, ('wave', 'xyz'))
results.createVariable('weight', N.Int32, ('wave',))
results.createVariable('q', N.Float32, ('wave', 'xyz'))
results.variables['q'].units = "nm-1"
results.createVariable('f_acoustic', N.Float32, ('wave', 'atoms',
                                                 'xyz1', 'xyz2'))
results.variables['f_acoustic'].units = "nm2"
results.createVariable('f_optical', N.Float32, ('wave', 'atoms',
                                                'xyz1', 'xyz2'))
results.variables['f_optical'].units = "nm2"
results.createVariable('force_constant', N.Float32, ('wave', 'mode'))
results.variables['force_constant'].units = "kJ mol-1 nm-2"

# Loop over all q vectors.
q_index = 0
for (i1, i2, i3), w in q_vectors.items():

    # Calculate the q vector and the phase factor for each atom
    # in the unit cell.
    q = 2.*N.pi*(i1*b1/n1 + i2*b2/n2 + i3*b3/n3)
    q_even = (2*i1) % n1 == 0 and (2*i2) % n2 == 0 and (2*i3) % n3 == 0
    p = N.exp(1j*(universe.configuration()*q).array)
    p = p[N.NewAxis, :, N.NewAxis]

    # Calculate the upper half of the complex K matrix.
    k = fc.array + 0j
    for i in range(natoms):
        for j in range(i+1, natoms):
            v = universe.distanceVector(atoms[i], atoms[j])
            d = v.length()
            if d < cutoff:
                k[atoms[i].index, :, atoms[j].index, :] *= N.exp(1j*(q*v))

    # Symmetrize the complex K matrix.
    n = k.shape[0]
    for i in range(n):
        for j in range(i+1, n):
            k[j, :, i, :] = N.transpose(N.conjugate(k[i, :, j, :]))

    # Calculate the eigenvalues and eigenvectors.
    k.shape = (3*n, 3*n)
    ev, modes = LA.Heigenvectors(k)
    ev = ev.real
    modes.shape = (3*n, n, 3)

    # Multiply the eigenvectors by the phase factors.
    N.multiply(modes, p, modes)

    # For q == 0, ignore the acoustic modes which are in fact the
    # global translations of the crystal.
    if q.length() == 0.:
        first_mode = 3
    else:
        first_mode = 0
    # Index of modes sorted by increasing eigenvalue.
    sort_index = N.argsort(ev)

    # Loop over the modes and add up the contributions to the
    # atomic fluctuations.
    f_acoustic = ParticleTensor(universe)
    f_optical = ParticleTensor(universe)
    for i in range(first_mode, len(modes)):
        if i < 3:
            f = f_acoustic
        else:
            f = f_optical
        index = sort_index[i]
        rm = modes[index].real
        im = modes[index].imag
        rms = N.sum(N.sum(rm*rm))
        ims = N.sum(N.sum(im*im))
        if q_even:
            scale = ev[index] * (n1*n2*n3) * rms
        else:
            scale = ev[index] * (n1*n2*n3) * (rms+ims)/2.
        f.array += (rm[ :, :, N.NewAxis]*rm[:, N.NewAxis, :])/scale
        if w == 2:
            if q_even:
                scale = ev[index] * (n1*n2*n3) * ims
            else:
                scale = ev[index] * (n1*n2*n3) * (rms+ims)/2.
            f.array += (im[ :, :, N.NewAxis]*im[:, N.NewAxis, :])/scale
    # Apply the Boltzmann factor
    f_acoustic.array *= temperature*Units.k_B
    f_optical.array *= temperature*Units.k_B

    # Store the results in the output file.
    results.variables['index'][q_index, :] = N.array([i1, i2, i3], N.Int32)
    results.variables['weight'][q_index] = w
    results.variables['q'][q_index, :] = N.array(q, N.Float32)
    results.variables['f_acoustic'][q_index, :, :, :] = \
                      N.take(N.array(f_acoustic.array, N.Float32), atom_indices)
    results.variables['f_optical'][q_index, :, :, :] = \
                      N.take(N.array(f_optical.array, N.Float32), atom_indices)
    results.variables['force_constant'][q_index, :] = \
                      N.take(N.array(ev, N.Float32), sort_index)
    results.sync()

    q_index += 1

# Close the output file.
results.close()
