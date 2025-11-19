from aiida import load_profile

from aiida.orm import load_node, load_code, load_group, Str, Bool, Int, Float, Dict, KpointsData
from aiida.engine import submit

from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import WorkflowFactory

from aiida_vasp_ext_mlauer.util.jhpc_computer_options import get_jhpc_options
from aiida_vasp_ext_mlauer.workflows.vol_opt import VolOptWorkChain

load_profile('lauerm-prod')

group = load_group('nitride_exploration/blk_materials')


def perform_initial_vasp_calc(code_str, structure_pk, parameters, pwcut, kmesh, potential_family, potential_mapping, options):


    # load nodes
    code = load_code(code_str)
    structure = load_node(structure_pk)

    workchain = WorkflowFactory('vasp.vasp')
    

    # Set Up Settings
    settings = AttributeDict()

    # make copy of parameters and insert encut
    parameters = parameters.copy()
    parameters["incar"].update({"encut": pwcut})


    # set up inputs
    inputs = AttributeDict()

    # General
    inputs.structure = structure 
    
    inputs.kpoints = KpointsData()
    inputs.kpoints.set_kpoints_mesh(kmesh)



    node = submit(workchain, **inputs)
    return node


if __name__ == '__main__':
    CODE_STRING = "vasp6.4@jhpc"

    PARAMETERS = AttributeDict()
    PARAMETERS.vasp = {
        'incar': {
            "istart": 0,
            "icharg": 2,
            "ismear": -5,
            "gga": "PS"
            }}
    
    PARAMETERS.band = {
        'incar': {
            "istart": 0,
            "icharg": 2,
            "ismear": -5,
            "gga": "PS"
            }}

    POTENTIAL_FAMILY = "PBE.64"
    POTENTIAL_MAPPING = {
        "Al": "Al",
        "Sc": "Sc",
        "Ga": "Ga",
        "In": "In",
        "N": "N"
                         }
    
    OPTIONS = AttributeDict()
    OPTIONS.step_1 = get_jhpc_options(1, 'debug')
    OPTIONS.step_2 = get_jhpc_options(1, 'debug')


    STRUCTURE_PKs = [
    688, # wz-AlN
    689, # wz-ScN
    690, # wz-GaN
    691, # wz-InN
    ]

    STRUCTURE_GENERATION = AttributeDict()
    STRUCTURE_GENERATION.parameters = Dict(dict={"mode": "uniform"})
    STRUCTURE_GENERATION.max_volume_scaling = Float(0.1)
    STRUCTURE_GENERATION.lattice_points = Int(5)

    KMESH = [10, 10, 6]
    PWCUTOFF = 550

    for STRUCTURE_PK in STRUCTURE_PKs:
        pass 
        # node = perform_volopt(CODE_STRING, STRUCTURE_PK, PARAMETERS, PWCUTOFF, KMESH, RELAX, POTENTIAL_FAMILY, POTENTIAL_MAPPING, STRUCTURE_GENERATION, OPTIONS)

        # group.add_nodes(node)

    for pk in [4525, 4554, 4626]:
        node = restart_failed_vol_opt(pk, relax=RELAX, options=OPTIONS)
        group.add_nodes(node)

