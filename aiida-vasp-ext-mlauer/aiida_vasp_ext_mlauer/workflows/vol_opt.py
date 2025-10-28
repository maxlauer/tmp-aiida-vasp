from aiida.engine import WorkChain
from aiida.plugins import DataFactory, WorkflowFactory

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

    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(cls._step_1_workchain, namespace='step_1', exclude=['structures', 'workchain_metadata'])
        spec.expose_inputs(cls._step_2_workchain, namespace='step_2', exclude=['structure'])
        
        # Define general Inputs
        spec.input("structure",
                   valid_type=DataFactory("structure"),
                   help="Input structure to be volume optimized."
                   )

        # Define general Metadata Inputs - I am doing this to group metadata together .. don't know if that is a good idea, or exposing it would be beter #QUESTION
        spec.input_namespace(name='workchain_metadata',
                             dynamic = False,
                             help='Custom namespace for WorkChain-specific metadata (e.g. verbosity, plotting, tags).'
                            )   

        # Define Inputs realted to EoS WorkChain (Step 1)
        # ARE THESE METADATA OR JUST INPUTS? #QUESTION - prolly not, as they would be needed in order to redo the calculation
        spec.input("step_1.structure_generation.mode",
                   valid_type=DataFactory("core.str"),
                   default=DataFactory("core.str")("uniform"),
                   help="Mode for structure generation for volume optimization."
                   )
        spec.input("step_1.structure_generation.max_lattice_scaling",  # QUESTION
                   valid_type=DataFactory("core.float"),
                   default=DataFactory("core.float")(0.1),
                   help="Maximum lattice scaling factor for volume optimization. in percent (e.g. 0.1 = 10%) of input.structure lattice constant."
                   )
        spec.input("step_1.structure_generation.lattice_points",  # QUESTION
                   valiid_type=DataFactory("core.int"),
                #    validator= lambda x: x.value >= 3 and x.value % 2 == 1,  #  enforce minimum value and odd number #QUESTION
                   default=DataFactory("core.int")(5),
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