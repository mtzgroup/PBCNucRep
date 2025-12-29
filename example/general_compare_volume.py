import numpy as np
import pyscf
from pyscf.dft.gen_grid import gen_atomic_grids
from pyscf.dft.radi import treutler_ahlrichs

cell = pyscf.M(
    a = np.eye(3) * 3.5668,
    atom = '''
        C     0.      0.      0.    
        C     0.8917  0.8917  0.8917
        C     1.7834  1.7834  0.    
        C     2.6751  2.6751  0.8917
        C     1.7834  0.      1.7834
        C     2.6751  0.8917  2.6751
        C     0.      1.7834  1.7834
        C     0.8917  2.6751  2.6751
    ''',
    basis = '6-31g',
    verbose = 4,
)

atom_grids_tab = gen_atomic_grids(cell, atom_grid = (10, 14), radi_method = treutler_ahlrichs, level = None, prune = None)
atm_coords = cell.atom_coords()
atomic_radii = pyscf.data.radii.BRAGG
lattice_vector = cell.lattice_vectors().flatten()

filename = "compare_volume.txt"

with open(filename, "w") as f:
    line = ""
    for i in range(9):
        line += f"{lattice_vector[i]:.16f} "
    line += "\n"
    f.write(line)

    for i_atom in range(cell.natm):
        line = f"{atm_coords[i_atom, 0]:.16f} {atm_coords[i_atom, 1]:.16f} {atm_coords[i_atom, 2]:.16f} {atomic_radii[i_atom]:.16f}\n"
        f.write(line)

    for i_atom in range(cell.natm):
        coords, vol = atom_grids_tab[cell.atom_symbol(i_atom)]
        coords = coords + atm_coords[i_atom]
        for i_grid in range(coords.shape[0]):
            line = f"{coords[i_grid, 0]:.16f} {coords[i_grid, 1]:.16f} {coords[i_grid, 2]:.16f} {vol[i_grid]:.16f} {i_atom}\n"
            f.write(line)

