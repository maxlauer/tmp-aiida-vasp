"""
Call script to calculate the total energies for different volumes of the silicon structure.

This particular call script set up a standard calculation for each structure and submits
each of them.
"""
# pylint: disable=too-many-arguments
import numpy as np

from aiida import load_profile
from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Str, Bool, load_code, load_group
from aiida.engine import submit 

load_profile('lauerm-test')


def get_structure(alat, lattice, basis, direct=True):

    Structure = DataFactory('core.structure')
    lattice = alat * lattice 
    structure = Structure(cell = lattice)

    for atm, pos in zip(*basis):
        if direct:
            pos = np.dot(pos, lattice)
        structure.append_atom(position = pos, symbols = atm)
        
    return structure
    

def main(code_string, incar, kmesh, structure, potential_family, potential_mapping, options):
    
    # Get AiiDA datatypes
    dict_data = DataFactory('core.dict')
    kpoints_data = DataFactory('core.array.kpoints')

    # Select Workchain
    workchain = WorkflowFactory('vasp.vasp')

    # Declare options, settings and input containers
    settings = AttributeDict()
    inputs = AttributeDict()

    inputs.code = load_code(code_string)
    inputs.structure = structure 

    kpoints = kpoints_data()
    kpoints.set_kpoints_mesh(kmesh)
    inputs.kpoints = kpoints 

    inputs.parameters = dict_data(dict=incar)

    inputs.potential_family = Str(potential_family)
    inputs.potential_mapping = dict_data(dict=potential_mapping)

    inputs.options = dict_data(dict = options)
    inputs.settings = dict_data(dict = settings)

    inputs.verbose = Bool(True)


    node = submit(workchain, **inputs)
    return node 


if __name__ == '__main__':

    GROUP = load_group('aiida-vasp-tests')
    # 
    CODE_STRING = 'vasp6@jhpc'

    INCAR = { 
        'incar': {
            "istart":  0,
            "icharg": 2,
            "encut":   250,
            "ismear":  0,
            "sigma":   0.1
        }
    }
    
    KMESH = [11, 11, 11]
    POT_FAM = "PBE.64"
    POT_MAP = {'Si': 'Si'}

    LAT = np.array([
        [1/2, 1/2, 0],
        [0, 1/2, 1/2],
        [1/2, 0, 1/2]
    ])
    BASE_POS = np.array([[0, 0, 0]])
    BASE_ATM = ["Si"]
    BASE = (BASE_ATM, BASE_POS)


    OPTIONS = AttributeDict()
    OPTIONS.queue_name = 'debug'
    OPTIONS.max_wallclock_seconds = 1800
    OPTIONS.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 10}



    ALATS = np.linspace(3.5, 4.4, 10)
    for ALAT in ALATS:
        STRUCTURE = get_structure(ALAT, LAT, BASE)
        
        NODE = main(CODE_STRING, INCAR, KMESH, STRUCTURE, POT_FAM, POT_MAP, OPTIONS)
        GROUP.add_nodes(NODE)        




