#!/usr/bin/env python
from ase import io
import sys
import numpy as np 
from gpaw.utilities.ewald import Ewald 

def get_lattice_energy(str,q_dict, G=5, Ng=[9, 9, 9], Nl=[3, 3, 3]):
  charges = np.zeros(len(str))

  for i in range(len(str)):
    symbol = str[i].symbol
    charges[i] =  q_dict[symbol]

  k = 8.987e9
  e0 = 1.602e-19
  cell = str.cell
  basis = str.positions
  ewald = Ewald(cell, G=G, Ng=Ng, Nl=Nl) 
  U = []
  for i, r in enumerate(basis): 
    u = charges[i] * ewald.get_electrostatic_potential(r, basis, charges, excludefroml0=i) * 0.5e10 * e0  * k 
    print i, u
    U.append(u)
  return U


if __name__ == '__main__':
  '''
  file = sys.argv[1]
  q_dict = {'Cu':1, 'Sn':4, 'S':-2, 'Zn':2, 'In':3, 'O':-2}
  str = io.read(file, format = 'vasp_out')
  U = get_lattice_energy(str, q_dict)
  Utot  = np.sum(U)
  import pickle
  print Utot
  f = open('madelung.pckl', 'w')
  pickle.dump([Utot, U], f)
  f.close()
  '''
  from ase.lattice.spacegroup import crystal
  a = 4.6
  c = 2.95
  str =crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
  q_dict = {'Ti': 4, 'O':-2}
  U = get_lattice_energy(str, q_dict)
  print 'TiO2', U

  a = 5.64
  str  = crystal(['Na', 'Cl'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
               cellpar=[a, a, a, 90, 90, 90])
  q_dict = {'Na': 1, 'Cl':-1}
  U = get_lattice_energy(str, q_dict)
  print 'NaCl', U


