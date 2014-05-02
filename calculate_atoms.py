from ase import Atoms
from ase.calculators.vasp import Vasp, int_keys, bool_keys, list_keys
import os



elements = ['Ru', 'Pt', 'La', 'Ni', 'Pd', 'Pb', 'Nb', 'Rh', 'In', 'Te', 'O', 'Sr', 'K', 'Na', 'Sn', 'N', 'Sc', 'Fe', 'Ca', 'Rb', 'Mg', 'W', 'Mn', 'P', 'S', 'Ba', 'Bi', 'Sb', 'Mo', 'Cr', 'Se', 'Y', 'Hg', 'Zn', 'Co', 'Zr', 'Ag', 'Ir', 'Al', 'Hf', 'Cd', 'As', 'Au', 'Ga', 'Cu']


setups = {'Pb':'_d', 'Nb': '_pv', 'In': '_d', 'O': '_s', 'Sr': '_sv', 'K': '_sv', 'Na': '_pv', 'Sn': '_d', 'N': '_s', 'Sc': '_sv', 'Ca': '_pv', 'Rb': '_sv', 'W': '_pv', 'Ba': '_sv', 'Bi': '_d', 'Mo': '_pv', 'Cr': '_pv', 'Y': '_sv', 'Zr': '_sv', 'Hf': '_pv', 'Ga': '_d'}


for el in elements:
    os.system('mkdir -p %s'%(el))
    os.chdir('%s'%(el))


    atoms = Atoms(el, cell=(10.0,10.1,10.2))
    calc = Vasp(xc='PBE',
            gga='PE',
            kpts=(1,1,1),
            gamma=True,
            lplane=True,
            ismear=0,
            lorbit=10,
            addgrid=True,
            ibrion=0,
            ispin=2,
#            hund=True,
            sigma=0.0,
            encut=520,
            ediff=1e-06,
            setups=setups,
            prec='accurate',
            npar=4,
            lvhar=False,
            lvtot=False,
            algo='fast',
            lwave=False,
           )
    atoms.set_calculator(calc)
    atoms.get_potential_energy()


    os.chdir('../')
