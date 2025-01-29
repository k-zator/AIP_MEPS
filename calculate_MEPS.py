#!/bin/python
import os
import sys
import numpy as np
from pyscf import __config__
setattr(__config__, 'cubegen_box_margin', 4.0) ## set the padding, best fit for 6-31G* set
from pyscf import gto
from pyscf import dft
from pyscf.geomopt.geometric_solver import optimize
from pyscf.tools import cubegen
from pyscf.tools.cubegen import Cube

# constants
resolution = (1/12)
isosurface = [0.002, 0.0104, 0.0300]
tol = [0.00003, 0.00003, 0.0001]


def write_MEPS(self, field, fname, comment=None):
    """Short function to write out custom defined ".cube" file."""
    if comment is None:
        comment = 'Generic field? Supply the optional argument "comment" to define this line'

    mol = self.mol
    coord = mol.atom_coords()
    with open(fname, 'w') as f:
        f.write(comment+'\n')
        f.write(f'PySCF Custom Script for writing MEPS\n')
        f.write(f'{mol.natm:5d}')
        f.write('%12.6f%12.6f%12.6f\n' % tuple(self.boxorig.tolist()))
        dx = self.xs[-1] if len(self.xs) == 1 else self.xs[1]
        dy = self.ys[-1] if len(self.ys) == 1 else self.ys[1]
        dz = self.zs[-1] if len(self.zs) == 1 else self.zs[1]
        delta = (self.box.T * [dx,dy,dz]).T
        f.write(f'{self.nx:5d}{delta[0,0]:12.6f}{delta[0,1]:12.6f}{delta[0,2]:12.6f}\n')
        f.write(f'{self.ny:5d}{delta[1,0]:12.6f}{delta[1,1]:12.6f}{delta[1,2]:12.6f}\n')
        f.write(f'{self.nz:5d}{delta[2,0]:12.6f}{delta[2,1]:12.6f}{delta[2,2]:12.6f}\n')
        for ia in range(mol.natm):
            atmsymb = mol.atom_symbol(ia)
            f.write('%5d%12.6f'% (gto.charge(atmsymb), 0.))
            f.write('%12.6f%12.6f%12.6f\n' % tuple(coord[ia]))

        for i in range(len(field)):
            fmt = '  %13.5E  ' * len(field[i]) + '\n'
            f.write(fmt % tuple(field[i].tolist()))


def MEPS_calculator(mol_input, isosurface, tol):
    """Short function to calculate MEPS at three custom isosurfaces chosen for the AIP calculation"""

    ## getting the .xyz format just in case other was supplied
    file_ext = mol_input.split(".")[-1]
    name, _ = mol_input.split(f".{file_ext}")
    if file_ext != "xyz":
        os.system(f"obabel {mol_input} -O{name}.xyz; fi")
    ## build molecule and obtain its DFT energy and wavefunction
    mol = gto.M(atom=f"{name}.xyz")
    mol.basis="6-31G*"
    mol.verbose = 0
    mol.build();
    method = dft.RKS(mol)
    method.xc = 'b3lypg'
    method.scf();

    ## in principle this gives complete density and complete MEPS
    cc = Cube(mol, resolution=resolution)
    a = cubegen.density(mol, f'{name}_den.cube', method.make_rdm1(), resolution=resolution)
    b = cubegen.mep(mol, f'{name}_meps.cube', method.make_rdm1(), resolution=resolution)

    ## create an isosurface mask to select custom isosurface
    coords = cc.get_coords() ## coords of the cube voxels
    den = a.reshape(cc.nx*cc.ny*cc.nz) ## density ordered in line with coords
    meps = b.reshape(cc.nx*cc.ny*cc.nz)

    for i, iso in enumerate(isosurface):
        iso_mask = np.isclose(den, iso, atol=tol[i])

        coord_iso = coords[iso_mask]
        den_iso = den[iso_mask]

        ## and hence selected MEPS values
        meps_iso = meps[iso_mask]

        ## and write output
        field = np.hstack((coord_iso, meps_iso.reshape(len(meps_iso),1)))
        write_MEPS(cc, field, f"{name}_{iso}.cube", 'Molecular electrostatic potential in real space')


## and execute the functions
mol_input = str(sys.argv[1])
MEPS_calculator(mol_input, isosurface=isosurface, tol=tol)