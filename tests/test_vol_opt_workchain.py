"""
# FOR TODAY - 
 1. TEST AND DEBUG VOLOPT WORKCHAIN 
 2. PERFORM CALCULATIONS FOR AlN, GaN, InN and ScN
    1. CONVERGENCE FOR ENERGY
    2. VOLOPT WORKCHAIN
    3. CONVERGENCE FOR BANDGAP (IF POSSIBLE, IF NOT - IGNORE FOR NOW) - END OF THIS ALSO GIVES BANDSTRUCTURE CALC
"""
import numpy as np

from aiida.engine import submit

from aiida.orm import Bool, Dict, KpointsData, Code, Str

from aiida_vasp_ext_mlauer.util.structures import gen_aiida_structure
from aiida.common.extendeddicts import AttributeDict

from aiida.plugins import DataFactory

from aiida_vasp_ext_mlauer.workflows.vol_opt import VolOptWorkChain


def perform_volopt(code_string, incar, kmesh, structure, structure_generation, relax, potential_family, potential_mapping, options):
    workchain = VolOptWorkChain

    # Set Up Settings
    settings = AttributeDict()

    eos_settings = AttributeDict()

    relax_settings = AttributeDict()
    relax_settings.parser_settings = AttributeDict()
    relax_settings.ADDITIONAL_RETRIEVE_LIST = ['CHGCAR']


    # set up inputs
    inputs = AttributeDict()
    inputs.workchain_metadata = AttributeDict()
    inputs.workchain_metadata.create_plot = True

    # General
    # Structures 
    inputs.structure = structure 
    
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh(kmesh)

    # Step 1

    inputs.step_1 = AttributeDict()
    # VolOptWorkchain inputs
    inputs.step_1.structure_generation = structure_generation

    # EoSWorkchain inputs
    inputs.step_1.minimum_mode = Str("Murnaghan")

    # VASP INPUTS
    inputs.step_1.code = Code.get_from_string(code_string)
    inputs.step_1.kpoints = kpoints 

    inputs.step_1.parameters = Dict(dict=incar)
    inputs.step_1.potential_family = Str(potential_family)
    inputs.step_1.potential_mapping = Dict(dict=potential_mapping)
    inputs.step_1.options = Dict(dict=options)
    inputs.step_1.settings = Dict(dict=eos_settings)
    inputs.step_1.verbose = Bool(True)



    # Step 2
    inputs.step_2 = AttributeDict()

    
    # VASP INPUTS
    inputs.step_2.code = Code.get_from_string(code_string)
    inputs.step_2.kpoints = kpoints 

    inputs.step_2.parameters = Dict(dict=incar)
    inputs.step_2.potential_family = Str(potential_family)
    inputs.step_2.potential_mapping = Dict(dict=potential_mapping)
    inputs.step_2.options = Dict(dict=options)
    inputs.step_2.settings = Dict(dict=relax_settings)
    inputs.step_2.verbose = Bool(True)

    inputs.step_2.relax = relax


    submit(workchain, **inputs)


if __name__ == '__main__':
    # Code_string is chosen among the list given by 'verdi code list'
    CODE_STRING = 'vasp6.4@jhpc'

    # INCAR equivalent
    # Set input parameters
    INCAR = {'incar': {'istart': 0, 'icharg': 2, 'encut': 240, 'ismear': 0, 'sigma': 0.1}}

    # KPOINTS equivalent
    # Set kpoint mesh
    KMESH = [11, 11, 11]

    # POTCAR equivalent
    # Potential_family is chosen among the list given by
    # 'verdi data vasp-potcar listfamilies'
    POTENTIAL_FAMILY = 'PBE.64'
    # The potential mapping selects which potential to use, here we use the standard
    # for silicon, this could for instance be {'Si': 'Si_GW'} to use the GW ready
    # potential instead
    POTENTIAL_MAPPING = {'Si': 'Si'}

    # jobfile equivalent
    # In options, we typically set scheduler options.
    # See https://aiida.readthedocs.io/projects/aiida-core/en/latest/scheduler/index.html
    # AttributeDict is just a special dictionary with the extra benefit that
    # you can set and get the key contents with mydict.mykey, instead of mydict['mykey']
    OPTIONS = AttributeDict()
    OPTIONS.resources = {'num_machines': 1, 'num_mpiprocs_per_machine': 4}
    OPTIONS.queue_name = 'debug'
    # OPTIONS.qos = "short"
    OPTIONS.max_wallclock_seconds = 1800

    # POSCAR equivalent
    # Set up Silicon Structure
    LAT_CONST = 3.9 # Angstrom
    LATTICE = np.array([[.5, .5, 0], [0, .5, .5], [.5, 0, .5]])
    SPECIES = ["Si"]
    POSITION = np.array([[0,0,0]])

    STRUCTURE = gen_aiida_structure(LAT_CONST, LATTICE, SPECIES, POSITION)

    STRUCTURE_GENERATION = AttributeDict()
    STRUCTURE_GENERATION.parameters = DataFactory("core.dict")(dict={"mode": "uniform"})
    STRUCTURE_GENERATION.max_volume_scaling = DataFactory("core.float")(0.1)
    STRUCTURE_GENERATION.lattice_points = DataFactory("core.int")(5)

    RELAX = {
    "perform": Bool(True),
    "force_cutoff": DataFactory('core.float')(1e-3),
    "steps": DataFactory('core.int')(15),
    "positions": Bool(True),
    "shape": Bool(True),
    "volume": Bool(True),
    }

    perform_volopt(CODE_STRING, INCAR, KMESH, STRUCTURE, STRUCTURE_GENERATION, RELAX, POTENTIAL_FAMILY, POTENTIAL_MAPPING, OPTIONS)