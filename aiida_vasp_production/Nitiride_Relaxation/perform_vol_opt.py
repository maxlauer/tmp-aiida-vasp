from aiida import load_profile

from aiida.orm import load_node, load_code, load_group, Str, Bool, Int, Float, Dict, KpointsData
from aiida.engine import submit

from aiida.common.extendeddicts import AttributeDict

from aiida_vasp_ext_mlauer.util.jhpc_computer_options import get_jhpc_options
from aiida_vasp_ext_mlauer.workflows.vol_opt import VolOptWorkChain

load_profile('lauerm-prod')

group = load_group('nitride_exploration/blk_materials')


def perform_volopt(code_str, structure_pk, parameters, pwcut, kmesh, relax, potential_family, potential_mapping, structure_generation, options):


    # load nodes
    code = load_code(code_str)
    structure = load_node(structure_pk)

    workchain = VolOptWorkChain
    

    # Set Up Settings
    settings = AttributeDict()

    eos_settings = AttributeDict()

    relax_settings = AttributeDict()
    relax_settings.parser_settings = AttributeDict()

    # make copy of parameters and insert encut
    step_1_para = parameters.step_1.copy()
    step_2_para = parameters.step_1.copy()
        
    step_1_para["incar"].update({"encut": pwcut})
    step_2_para["incar"].update({"encut": pwcut})


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
    inputs.step_1.code = code
    inputs.step_1.kpoints = kpoints 

    inputs.step_1.parameters = Dict(dict=step_1_para)
    inputs.step_1.potential_family = Str(potential_family)
    inputs.step_1.potential_mapping = Dict(dict=potential_mapping)
    inputs.step_1.options = Dict(dict=options.step_1)
    inputs.step_1.settings = Dict(dict=eos_settings)
    inputs.step_1.verbose = Bool(True)



    # Step 2
    inputs.step_2 = AttributeDict()

    
    # VASP INPUTS
    inputs.step_2.code = code
    inputs.step_2.kpoints = kpoints 

    inputs.step_2.parameters = Dict(dict=step_2_para)
    inputs.step_2.potential_family = Str(potential_family)
    inputs.step_2.potential_mapping = Dict(dict=potential_mapping)
    inputs.step_2.options = Dict(dict=options.step_2)
    inputs.step_2.settings = Dict(dict=relax_settings)
    inputs.step_2.verbose = Bool(True)

    inputs.step_2.relax = relax


    node = submit(workchain, **inputs)
    return node

def restart_failed_vol_opt(process_pk, relax, options):

    failed_process = load_node(process_pk)
    builder = failed_process.get_builder_restart()

    builder.step_1.options = Dict(dict=options.step_1)

    builder.step_2.options = Dict(dict=options.step_2)
    builder.step_2.relax = relax 

    new_calc = submit(builder)
    return new_calc



if __name__ == '__main__':
    CODE_STRING = "vasp6.4@jhpc"

    PARAMETERS = AttributeDict()
    PARAMETERS.step_1 = {
        'incar': {
            "istart": 0,
            "icharg": 2,
            "ismear": -5,
            "gga": "PS"
            }}
    PARAMETERS.step_2 = {
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

    RELAX = {
        'perform': Bool(True),
        'force_cutoff': Float(1e-2),
        'steps': Int(60),
        'positions': Bool(True),
        'shape': Bool(True),
        'volume': Bool(True)
    }

    for STRUCTURE_PK in STRUCTURE_PKs:
        pass 
        # node = perform_volopt(CODE_STRING, STRUCTURE_PK, PARAMETERS, PWCUTOFF, KMESH, RELAX, POTENTIAL_FAMILY, POTENTIAL_MAPPING, STRUCTURE_GENERATION, OPTIONS)

        # group.add_nodes(node)

    for pk in [4525, 4554, 4626]:
        node = restart_failed_vol_opt(pk, relax=RELAX, options=OPTIONS)
        group.add_nodes(node)

