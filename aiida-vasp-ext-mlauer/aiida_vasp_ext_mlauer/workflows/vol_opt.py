import numpy as np

from aiida import load_profile

from aiida.engine import WorkChain, calcfunction
from aiida.plugins import DataFactory, WorkflowFactory

from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.workchains import compose_exit_code

from aiida_vasp_ext_mlauer.util.functions import powerrange

load_profile('lauerm-test')

"""
todo:
- too much
- think about maybe making a multi step workchain an option - idea - first scan, second more targeted eos scan and then relax - later

questions:
- should things like max_lattice_scaling and lattice_points be part of the step_1 namespace? - yes - they are definitly inputs
    - are they part of metadata? or just inputs?
- should workchain metadata be exposed, or focused together with the workchain_metadata of the higher level workchain? - no - keep together at the highest level, and supply them downwards

- what exit codes do I need? - I think the NO_CALLED_xxx_WORKCHAIN + the error of the underlying workchains should be enough
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
        spec.input('workchain_metadata.volume_precision', 
                   valid_type = int,
                   default = lambda: 3,
                   help = "Precision decimal places for volume convergence in \u213B\u00B3 in VolOptWorkChain. Default is 3",
                   non_db = True
                   )


        # Define Inputs realted to EoS WorkChain (Step 1)
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
        
        spec.input('workchain_metadata.eos_metadata.create_plot',
                   valid_type = bool,
                   default = False,
                   help = f"Toggle whether to create a plot from {cls._step_1_workchain}, or not. Default is False",
                   non_db = True
                   )
        
        


        # Define Exit Codes
        spec.exit_code(0, "NO_ERROR", message="Calculation ended without Error ... well done I guess")
        spec.exit_code(420, "ERROR_NO_CALLED_EOS_WORKCHAIN", message="No called EoSWorkChain detected")
        spec.exit_code(421, "ERROR_NO_CALLED_RELAX_WORKCHAIN", message="No called vasp.relax WorkChain detected")
        spec.exit_code(500, "ERROR_UNKOWN", message="Unkown Error detected in the VolOptWorkChain")


        # Define Outline
        spec.outline(
            cls.initialize,
            cls.init_and_run_step_1,
            cls.verify_step_1,
            cls.init_and_run_step_2,
            cls.verify_step_2,
            cls.finalize
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

        # Load Inital Structure in Context
        self.ctx.structure = self.inputs.structure

        ## Step 1
        self.ctx.eos = AttributeDict()
        self.ctx.eos.structure_generation = self.inputs.step_1.structure_generation
        self.ctx.eos_inputs = AttributeDict()

        self.ctx.eos.inputs.workchain_metadata = self.inputs.workchain_metadata.eos_metadata 


        # set up the result and finished flag
        self.ctx.eos.minimum = None
        # self.ctx.eos.finished = False #QUESTION - Is a finished flag necessary? doubt it


        ## Step 2
        self.ctx.relax = AttributeDict()
        self.ctx.relax.inputs = AttributeDict()


    def _init_inputs(self):
        """
        Set up verbosity value from inputs
        """

        try:
            self._verbose = self.inputs.verbose.value

        except AttributeError:
            pass

    
    def init_and_run_step_1(self):

        # Check existence of context inputs
        try:
            self.ctx.eos.inputs
        except AttributeError:
            raise ValueError("No input dictionary was defined in self.ctx.eos.inputs")
        
        # set up inputs for step 1
        self.ctx.eos.inputs.update(self.exposed_inputs(self._step_1_workchain))
        
        # Set up structures - #QUESTION - should I check if the parameters supplied contain a mode here again? if I have a validator that should be enough I guess
        self.report(f"Generating structures for {self._step_1_workchain.__name__} based on Structure {self.ctx.structure.pk} - Mode: {self.ctx.eos.structure_generation.parameters.get_dict().get('mode', 'uniform')}")
        self.ctx.eos.inputs.structures = generate_eos_structures(self.ctx.structure, self.ctx.eos.structure_generation)

        # run step 1
        running_eos = self.submit(self._step_1_workchain, **self.ctx.eos.inputs)
        self.report(f"Launching {self._step_1_workchain.__name__}<{running_eos.pk}>")
        self.ctx.eos.workchain = running_eos

    
    def verify_step_1(self):

        # Check if Step 1 WorkChain exists
        try:
            workchain = self.ctx.eos.workchain 
        except NameError:
            return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN
        
        step_1_exit_status = workchain.exit_status
        step_1_exit_message = workchain.exit_message

        if not step_1_exit_status:
            self.ctx.exit_code = self.exit_codes.NO_ERROR
            self.ctx.eos.minimum = workchain.outputs.eos_minimum.get_dict()

        else:
            self.ctx.exit_code = compose_exit_code(step_1_exit_status, step_1_exit_message)
            self.report(f"The called {workchain.__class__.__name__}<{workchain.pk}> returned a non-zero exit status.\n The exit status {self.ctx.exit_code} is inherited.")
            
        return self.ctx.exit_code


    def init_and_run_step_2(self):
        
        # Check existence of context inputs
        try:
            self.ctx.relax.inputs
        except AttributeError:
            raise ValueError("No input dictionary was defined in self.ctx.relax.inputs")
        
        # set up inputs for step 2
        self.ctx.relax.inputs.update(self.exposed_inputs(self._step_2_workchain))

        # Generate predicted Minimum Structure from Step 1
        eos_min_volume = np.round(self.ctx.eos.minimum.get('volume'), self.ctx.workchain_metadata.volume_precision)
        self.report(f"Generating predicted minimum structure for {self._step_2_workchain.__name__} based on inital Structure <{self.ctx.structure.pk}>, with Minimum Volume {eos_min_volume} \u213B\u00B3  from {self.ctx.eos.workchain.__class__.__name__} <{self.ctx.eos.workchain.pk}>")
        eos_relax_structure = self.ctx.structure.get_pymatgen().copy()
        
        eos_relax_structure.lattice = eos_relax_structure.lattice.matrix * (eos_min_volume / eos_relax_structure.volume)
        self.ctx.relax.inputs.structure = DataFactory("core.structure")(pymatgen_structure = eos_relax_structure)

        # run step 2
        running_relax = self.submit(self._step_2_workchain, **self.ctx.relax.inputs)
        self.report(f"Launching {self._step_2_workchain.__name__}<{running_relax.pk}>")
        self.ctx.relax.workchain = running_relax

    def verify_step_2(self):
        
        try:
            workchain = self.ctx.relax.workchain
        except NameError:
            return self.exit_codes.ERROR_NO_CALLED_RELAX_WORKCHAIN
        
        step_2_exit_status = workchain.exit_status
        step_2_exit_message = workchain.exit_message

        if not step_2_exit_status:
            self.ctx.exit_code = self.exit_codes.NO_ERROR
        else:
            self.ctx.exit_code = compose_exit_code(step_2_exit_status, step_2_exit_message)
            self.report(f"The called {workchain.__class__.__name__}<{workchain.pk}> returned a non-zero exit status.\n The exit status {self.ctx.exit_code} is inherited.")

        return self.ctx.exit_code

    def finalize(self):
        # What do I need to finalize? - is there anything, if I just pass the outputs from the sub workchains?
        # as outputs I probaly just pass along the outputs of the step 2 workchain?

        # Ooooh a plot - optionally 
        pass


        # What should the outputs be? - relaxed structure obviously - anything else? final energy/volume/bulk modulus from step 2?
        # anything else can be found in step 2 .. should I just expose the step 2 outputs, just add like the plot from step 1? prolly a good idea
        # maybe I should write a new plot, where the step_2 energy is also shown .. good Idea I guess



@calcfunction
def generate_eos_structures(inp_structure, structure_generation):
    structure = inp_structure.get_pymatgen()
    parameters = structure_generation.parameters.get_dict()
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
    pass