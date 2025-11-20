from aiida import load_profile

from aiida.orm import load_node, load_code, load_group, Str, Bool, Int, Float, Dict, KpointsData
from aiida.engine import submit

from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import WorkflowFactory

from aiida_vasp.workchains.bands import BandsWorkChain


from aiida_vasp_ext_mlauer.util.jhpc_computer_options import get_jhpc_options

load_profile('lauerm-prod')

group = load_group('nitride_exploration/blk_materials/lda_half')

def perform_inital_os(code_str, structure, parameters, pwcut, kmesh, potential_family, potential_mapping, options):

    # adjust parameters
    parameters = parameters.copy()
    parameters["incar"].update({"encut": pwcut})
    # prepare kpoints
    kpoints = KpointsData()
    kpoints.set_kpoints_mesh(kmesh)


    code = load_code(code_str)
    workchain = WorkflowFactory('vasp.vasp')

    settings = AttributeDict()
    inputs = AttributeDict()

    inputs.code = code

    inputs.structure = structure

    inputs.parameters = parameters
    inputs.kpoints = kpoints

    inputs.potential_family = Str(potential_family)
    inputs.potential_mapping = Dict(dict=potential_mapping)

    inputs.options = Dict(dict=options)

    inputs.settings = Dict(dict=settings)

    node = submit(workchain, **inputs)
    return node


def perform_vasp_bands_workchain(code_str, restart_folder, structure, parameters, pwcut, kmesh, potential_family, potential_mapping, options):


    # load nodes
    code = load_code(code_str)

    workchain =BandsWorkChain


    # Set Up Settings
    settings = AttributeDict()
    settings.parser_settings = {
        "add_fermi_level": True     # I think this doesn't do anything ...
    }

    # set up inputs
    inputs = AttributeDict()
    inputs.smearing = AttributeDict()
    inputs.bands = AttributeDict()

    # set up code
    inputs.code = code

    # Structure Setup
    inputs.structure = structure

    # set up parameters - not necessary - therefore I will not use it and instead leave it (probably takes parameters from restart_folder)
    # no it doesn't take it from restart .. .this is fishy, I will set it explicitly
    # make copy of parameters and insert encut
    parameters = parameters.copy()
    parameters["incar"].update({"encut": pwcut})

    inputs.parameters = parameters

    # k-points setup
    inputs.kpoints = KpointsData()
    inputs.kpoints.set_kpoints_mesh(kmesh)


    # POTCAR setup
    inputs.potential_family = Str(potential_family)
    inputs.potential_mapping = Dict(dict=potential_mapping)

    inputs.options = Dict(dict=options)

    inputs.settings = Dict(dict=settings)


    # Band specific inputs
    inputs.restart_folder = restart_folder

    # smearing inputs
    inputs.smearing.gaussian = Bool(True)
    inputs.smearing.sigma = Float(0.05)

    node = submit(workchain, **inputs)
    return node


def perform_bands_calculation():
    pass



if __name__ == '__main__':
    CODE_STRING = "vasp6.4@jhpc"

    PARAMETERS = {
            'incar': {
                "istart": 0,
                "icharg": 2,
                "gga": "HL"
                }}

    PARAMETERS_INITIAL = {
        'incar': {
            'istart': 0,
            'icharg': 2,
            'gga': "HL"
        }
    }

    POTENTIAL_FAMILY = "LDA-05.64"
    POTENTIAL_MAPPING = {
        "Al": "Al",
        "Sc": "Sc",
        "Ga": "Ga",
        "In": "In",
        "N": "N_05"
                         }

    OPTIONS = get_jhpc_options(1, 'debug')


    VOL_OPT_PKs = [
        1227, # wz-AlN
        1256, # wz-ScN
        1292, # wz-GaN
        1328, # wz-InN
    ]

    INITAL_VASP_PKs = [
        1806, # wz-AlN
        1863, # wz-ScN
        1828, # wz-GaN
        1838, # wz-InN
    ]

    KMESH = [10, 10, 6]
    PWCUTOFF = 550

    for idx in range(len(VOL_OPT_PKs)):
        if INITAL_VASP_PKs[idx] is not None:
            continue

        wc_pk = VOL_OPT_PKs[idx]
        wc = load_node(wc_pk)
        relax_wc = [sub_wc for sub_wc in wc.called if sub_wc.process_label == 'RelaxWorkChain'][-1]

        relaxed_str = relax_wc.outputs.relax.structure
        remote_folder = relax_wc.outputs.remote_folder

        node = perform_inital_os(
            code_str=CODE_STRING,
            structure=relaxed_str,
            parameters=PARAMETERS,
            pwcut=PWCUTOFF,
            kmesh=KMESH,
            potential_family=POTENTIAL_FAMILY,
            potential_mapping=POTENTIAL_MAPPING,
            options=OPTIONS
        )

        group.add_nodes(node)


    for idx in range(len(INITAL_VASP_PKs)):

        if INITAL_VASP_PKs[idx] is None:
            continue

        vol_opt_wc_pk = VOL_OPT_PKs[idx]
        init_wc_pk = INITAL_VASP_PKs[idx]

        vol_opt_wc = load_node(vol_opt_wc_pk)
        init_wc = load_node(init_wc_pk)

        relax_wc = [sub_wc for sub_wc in vol_opt_wc.called if sub_wc.process_label == 'RelaxWorkChain'][-1]

        relaxed_str = relax_wc.outputs.relax.structure
        remote_folder = init_wc.outputs.remote_folder

        node = perform_vasp_bands_workchain(
            code_str=CODE_STRING,
            restart_folder=remote_folder,
            structure=relaxed_str,
            parameters=PARAMETERS,
            pwcut=PWCUTOFF,
            kmesh=KMESH,
            potential_family=POTENTIAL_FAMILY,
            potential_mapping=POTENTIAL_MAPPING,
            options=OPTIONS
        )

        group.add_nodes(node)





