import numpy as np

from aiida import load_profile

from aiida.engine import WorkChain, calcfunction
from aiida.plugins import DataFactory, WorkflowFactory

from aiida.common.extendeddicts import AttributeDict

from aiida_vasp_ext_mlauer.util.structures import gen_aiida_structure
from aiida_vasp_ext_mlauer.util.functions import powerrange

load_profile('lauerm-test')

"""
todo:
- too much

questions:
- should things like max_lattice_scaling and lattice_points be part of the step_1 namespace?
    - are they part of metadata? or just inputs?
- should workchain metadata be exposed, or focused together with the workchain_metadata of the higher level workchain?

- what exit codes do I need?
- what outputs do I want to expose?
- are there better names for step_1 and step_2 - eos and relax maybe? but realx would then have relax.relax ..
- should I enforce odd numbered lattice points ? WHY? I'm not sure if that is good .. for now I won't
"""


class VolOptWorkChain(WorkChain):

    _verbose = False
    _step_1_workchain_string = "vasp_ext_ml.eos"
    _step_1_workchain = WorkflowFactory(_step_1_workchain_string)

    _step_2_workchain_string = "vasp.relax"
    _step_2_workchain = WorkflowFactory(_step_2_workchain_string)

    @classmethod
    def define(cls, spec):
        super(VolOptWorkChain, cls).define(spec)
        spec.expose_inputs(cls._step_1_workchain, namespace='step_1', exclude=['structures', 'workchain_metadata'])
        spec.expose_inputs(cls._step_2_workchain, namespace='step_2', exclude=['structure'])

        # Define general Inputs
        spec.input("structure",
                   valid_type=DataFactory("core.structure"),
                   help="Input structure to be volume optimized."
                   )

        # Define general Metadata Inputs - I am doing this to group metadata together .. don't know if that is a good idea, or exposing it would be beter #QUESTION
        spec.input_namespace(name='workchain_metadata',
                             dynamic = False,
                             help='Custom namespace for WorkChain-specific metadata (e.g. verbosity, plotting, tags).'
                            )

        # Define Inputs realted to EoS WorkChain (Step 1)
        # ARE THESE METADATA OR JUST INPUTS? #QUESTION - prolly not, as they would be needed in order to redo the calculation
        spec.input("step_1.structure_generation.parameters",
                   valid_type=DataFactory("core.dict"),
                   default=lambda: DataFactory("core.dict")(dict={'mode': "uniform"}),      # I don't know if this is good or not ... #QUESTION
                   validator = lambda x: isinstance(x, DataFactory("core.dict")) # Write Validator to check, if all parameters of the selected mode are present #TODO
                   help="Parameters for structure generation for volume optimization. - <mode> - required key: Modes: 'uniform', 'power'; 'power needs key <power_constant> as well'"
                   )

        spec.input("step_1.structure_generation.max_volume_scaling",  # QUESTION
                   valid_type=DataFactory("core.float"),
                   default=lambda: DataFactory("core.float")(0.1),
                   help="Maximum volume scaling factor for volume optimization. in percent (e.g. 0.1 = 10%) of input.structure volume."
                   )
        spec.input("step_1.structure_generation.lattice_points",  # QUESTION
                   valid_type=DataFactory("core.int"),
                #    validator= lambda x: x.value >= 3 and x.value % 2 == 1,  #  enforce minimum value and odd number #QUESTION
                   default=lambda: DataFactory("core.int")(5),
                   help="Number of lattice points to be calculated for the equation of state. Depending on mode used for structure generation might be + 1 (to include original structure)"
                   )

        # Define Metadata Inputs for Step 1

        spec.input_namespace(name='workchain_metadata.eos_metadata',
                             dynamic = False,
                             help='Custom namespace for EOS WorkChain-specific metadata (e.g. verbosity, plotting, tags).'
                            )
        spec.input('workchain_metadata.create_plot',
                   valid_type = bool,
                   default = False,
                   help = f"Toggle whether to create a plot from {cls._step_1_workchain}, or not. Default is False",
                   non_db = True
                   )


        # Define Exit Codes
        spec.exit_code(0, "NO_ERROR", message="Calculation ended without Error ... well done I guess")
        spec.exit_code(420, "ERROR_NO_CALLED_WORKCHAIN", message="No called WorkChain detected")
        spec.exit_code(500, "ERROR_UNKOWN", message="Unkown Error detected in the EOS WorkChain")


        # Define Outline
        spec.outline(
            cls.initialize
        )


    def initialize(self):
        """
        Set up context and inputs for the workchain steps
        """
        self._init_context()
        self._init_inputs()


    def _init_context(self):
        # set exit code to error - I don't understand quite why yet .. but ok
        self.ctx.exit_code = self.exit_codes.ERROR_UNKOWN

        # Set up structures (for now I will keep provenance here - don't know if necessary) #QUESTION
        self.ctx.structures = generate_eos_structures(self.inputs.step_1.structure_generation)

        ## What do I need to set up in the context
        # - I need to generate the structures - for that I need a function, and modes  ... probably a calcfunction? cause it creates the structures from data-nodes which are inserted ..
        # self.ctx should have
        #  - structures - list of structures to be used in step 1
        #  - eos_finished - bool to indicate if step 1 is finished (Needed?)
        #  - eos_results - to store results from step 1
        #  - eos_total_energies
        #  - eos_mode - to store the mode used for Equation of State
        #  - workchain_metadata
        #  - inputs
        #  -
        pass

    def _init_inputs(self):
        """
        Set up verbosity value from inputs
        """

        try:
            self._verbose = self.inputs.verbose.value

        except AttributeError:
            pass

        # lastly the outline - what do I have to do?
        # 1 - prepare step_1 - init step 1 - run step 1 - verify step 1
        # 2 - prepare step_2 - init step 2 - run step 2 - verify step 2
        # 3 - finalize the results at the highest level


        # What should the outputs be? - relaxed structure obviously - anything else? final energy/volume/bulk modulus from step 2?
        # anything else can be found in step 2 .. should I just expose the step 2 outputs, just add like the plot from step 1? prolly a good idea
        # maybe I should write a new plot, where the step_2 energy is also shown .. good Idea I guess






# I don't know if I might require other parameters .. probably - how should they be provided ... I guess as a namespace of structure_generation.mode_parameters
# @calcfunction
def generate_eos_structures(inp_structure, structure_generation):
    structure = inp_structure.get_pymatgen()
    parameters = structure_generation.parameters#.get_dict()
    mode = parameters.pop('mode')
    structure_name = parameters.get('str_name', 'Structure')

    max_vol_dev = structure_generation.max_volume_scaling
    num_lat_pts = structure_generation.lattice_points

    match mode:
        case 'uniform':
            volume_deviations = 1 + np.linspace(-max_vol_dev, max_vol_dev, num_lat_pts)
        case 'power':
            volume_deviations = 1 + powerrange(-max_vol_dev, max_vol_dev, num_lat_pts, power=parameters.get('power_constant', 1.5))
        

    structures = {}
    for vol_dev in volume_deviations:
        new_structure = structure.copy()
        new_structure.lattice = new_structure.lattice.matrix * vol_dev

        structures[f"{structure_name}_at_{vol_dev}_V0"] = DataFactory('core.structure')(pymatgen_structure = new_structure)

    return structures



if __name__ == '__main__':
    workchain = WorkflowFactory('vasp_ext_ml.vol_opt')
    builder = workchain.get_builder()

    from pymatgen.core.structure import Structure, Lattice

    # Set up structure for testing - hexagonal GaN
    a = 3.33
    c_over_a = 1.602

    lat = np.array([
        [1/2, -np.sqrt(3)/2,         0],
        [1/2,  np.sqrt(3)/2,         0],
        [  0,             0,  c_over_a]
        ])

    basis_pos = np.array([
        [1/3, 2/3, 0],
        [2/3, 1/3, 0],
        [1/3, 2/3, 0],
        [2/3, 1/3, 0]
        ])
    basis_atm = ["Ga", "Ga", "N", "N"]
    basis = (basis_atm, basis_pos)

    structure = gen_aiida_structure(a, lat, *basis)


    # set up builder.
    struc = structure.get_pymatgen()

    structure_generation = AttributeDict()
    structure_generation.parameters = {'mode': 'geometric'}
    structure_generation.max_volume_scaling = 0.1
    structure_generation.lattice_points = 5


    print(structure_generation.max_volume_scaling, structure_generation.parameters)
    print(generate_eos_structures(structure, structure_generation))
