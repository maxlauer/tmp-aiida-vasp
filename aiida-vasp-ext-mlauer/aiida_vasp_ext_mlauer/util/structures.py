import numpy as np

from aiida.plugins import DataFactory

# TODO - DOESN"T WORK WITH DISORDERED STRUCTURES YET - this could be problematic with species, if species is a disordered structure, then basis is a dictionary ..
def gen_aiida_structure(alat, lattice, species, basis, direct=True):

    Structure = DataFactory('core.structure')
    lattice = alat * lattice 
    structure = Structure(cell = lattice)

    for atm, pos in zip(species, basis):
        if direct:
            pos = np.dot(pos, lattice)
        structure.append_atom(position = pos, symbols = atm)
        
    return structure