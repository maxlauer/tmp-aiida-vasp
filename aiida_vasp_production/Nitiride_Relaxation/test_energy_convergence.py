
### INCAR
"""
ISTART = 0
ICHARGE = 2
ENCUT = ?
GGA = PS
ISMEAR = -5
ALGO = FAST
"""
from aiida import load_profile

from aiida.orm import load_node, load_code, load_group, Str, Bool, Int, Float, Dict
from aiida.plugins import WorkflowFactory
from aiida.engine import submit

from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.workchains.converge import ConvergeWorkChain
from aiida_vasp_ext_mlauer.util.jhpc_computer_options import get_jhpc_options


load_profile('lauerm-prod')

group = load_group('nitride_exploration/blk_materials/convergence')


def convergence_test(code_str, structure_pk, parameters, pot_fam, pot_map, options):

    # load nodes
    code = load_code(code_str)
    structure = load_node(structure_pk)

    # load workchain
    workchain = WorkflowFactory('vasp.converge')

    # Set up settings and code
    settings = AttributeDict()
    inputs = AttributeDict()

    # Define inputs
    inputs.code = code 
    inputs.structure = structure

    # INCAR and convergence setup
    inputs.parameters = parameters.incar
    inputs.converge = parameters.converge

    # POTCAR setup
    inputs.potential_family = Str(pot_fam)
    inputs.potential_mapping = Dict(dict=pot_map)

    inputs.options = Dict(dict=options)

    inputs.settings = Dict(dict=settings)

    inputs.verbose = Bool(True)

    inputs.relax = {
        "perform": Bool(False)
    }

    node = submit(workchain, **inputs)
    return node



if __name__ == '__main__':
    CODE_STRING = "vasp6.4@jhpc"

    PARAMETERS = AttributeDict()
    PARAMETERS.incar = {
        'incar': {
            "istart": 0,
            "icharg": 2,
            "ismear": -5,
            "gga": "PS"
            }}
    
    PARAMETERS.converge = AttributeDict()
    

    # These determin how the plane-wave cutoffs are sampled through
    PARAMETERS.converge.pwcutoff_start = Float(300.0) # 300 - 525 - in 25 eV steps
    PARAMETERS.converge.pwcutoff_step = Float(25.0) # 300 - 525 - in 25 eV steps
    PARAMETERS.converge.pwcutoff_samples = Int(15) # 300 - 525 - in 25 eV steps

    # These determin how the k-grids are sampled through
    PARAMETERS.converge.k_dense   = Float(0.07)
    PARAMETERS.converge.k_coarse   = Float(0.35)
    PARAMETERS.converge.k_samples = Int(14)

    PARAMETERS.converge.k_spacing = Float(0.1) #default k-spacing - only used for the energy convergence test,

            # These determine the cutoff criterion
    PARAMETERS.converge.cutoff_type  = Str('energy')
    PARAMETERS.converge.cutoff_value = Float(1e-4)
        

    POTENTIAL_FAMILY = "PBE.64"
    POTENTIAL_MAPPING = {
        "Al": "Al",
        "Sc": "Sc",
        "Ga": "Ga",
        "In": "In",
        "N": "N"
                         }
    
    OPTIONS = get_jhpc_options(1)

    STRUCTURE_PKs = [
    1186, # wz-AlN
    1187, # wz-ScN
    1188, # wz-GaN
    1189, # wz-InN
    ]

    for STRUCTURE_PK in STRUCTURE_PKs:
        node = convergence_test(CODE_STRING, STRUCTURE_PK, PARAMETERS, POTENTIAL_FAMILY, POTENTIAL_MAPPING, OPTIONS)

        group.add_nodes(node)

    
    # node = convergence_test(CODE_STRING, 5830, PARAMETERS, POTENTIAL_FAMILY, POTENTIAL_MAPPING, OPTIONS)