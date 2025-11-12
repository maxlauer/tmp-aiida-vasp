import numpy as np

from aiida import load_profile

from aiida.engine import WorkChain, calcfunction, append_
from aiida.plugins import DataFactory, WorkflowFactory

from aiida.common.extendeddicts import AttributeDict

from aiida_vasp.utils.workchains import prepare_process_inputs
from aiida_vasp.utils.workchains import compose_exit_code

from aiida_vasp_ext_mlauer.util.functions import powerrange
from aiida_vasp_ext_mlauer.workflows.eos import create_eos_plot
from aiida_vasp_ext_mlauer.util.orm import namespace_to_atrdict

"""
todo:
- make CHGCAR retriveal optional
- add ISIF 3 relax enforcement


- think about maybe making a multi step workchain an option - idea - first scan, second more targeted eos scan and then relax - later
- maybe add the option of a third step after the relaxation, in order to calculate certain properties, such as the band-structure, DOS
    - and also make optional outputs of the second step as outputs
- maybe include a optional step 0, which performs a convergence workchain to determine k-points and e-cutoff
- maybe wrap around BaseRestartWorkChain instead
questions:
- should things like max_lattice_scaling and lattice_points be part of the step_1 namespace? - yes - they are definitly inputs
    - are they part of metadata? or just inputs?
- should workchain metadata be exposed, or focused together with the workchain_metadata of the higher level workchain? - no - keep together at the highest level, and supply them downwards

- what exit codes do I need? - I think the NO_CALLED_xxx_WORKCHAIN + the error of the underlying workchains should be enough
- what outputs do I want to expose?
- are there better names for step_1 and step_2 - eos and relax maybe? but realx would then have relax.relax ..
- should I enforce odd numbered lattice points ? WHY? I'm not sure if that is good .. for now I won't
"""

_NECESSARY_PARAMETERS = {
    'uniform': [],
    'power': ['power_constant']
}



class VolOptWorkChain(WorkChain):

    _verbose = False
    _step_1_workchain_string = "vasp_ext_ml.eos"
    _step_1_workchain = WorkflowFactory(_step_1_workchain_string)
    _step_1_namespaces = ['dynamics', 'structures', 'workchain_metadata']

    _step_2_workchain_string = "vasp.relax"
    _step_2_workchain = WorkflowFactory(_step_2_workchain_string)
    _step_2_workchain_default_parser_settings = {"add_energies": True, "add_forces": True, "add_stress": True}
    _step_2_namespaces = ['dynamics', 'relax']

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
                   default = 3,
                   help = "Precision decimal places for volume convergence in \u213B\u00B3 in VolOptWorkChain. Default is 3",
                   non_db = True
                   )
        spec.input('workchain_metadata.create_plot',
                   valid_type = bool,
                   default = False,
                   help = "Toggle whether to create a plot, or not. Default is False",
                   non_db = True)



        # Define Inputs realted to EoS WorkChain (Step 1)
        spec.input("step_1.structure_generation.parameters",
                   valid_type=DataFactory("core.dict"),
                   default=lambda: DataFactory("core.dict")(dict={'mode': "uniform"}),      # I don't know if this is good or not ... #QUESTION
                   validator = validate_structure_generation_parameters, # Write Validator to check, if all parameters of the selected mode are present #TODO
                   help="Parameters for structure generation for volume optimization. - <mode> - required key: Modes: 'uniform', 'power'; 'power needs key <power_constant> as well'"
                   )

        spec.input("step_1.structure_generation.max_volume_scaling",  # QUESTION
                   valid_type=DataFactory("core.float"),
                   default=lambda: DataFactory("core.float")(0.05),
                   help="Maximum volume scaling factor for volume optimization. in percent (e.g. 0.5 = 5%) of input.structure volume."
                   )
        spec.input("step_1.structure_generation.lattice_points",  # QUESTION
                   valid_type=DataFactory("core.int"),
                   validator= lambda x, y: "Value 'lattice parameters' have to be greater than 3" if x.value < 3 and x % 2 == 0 else None,  #  enforce minimum value and odd number #QUESTION - Should I change this to check for oddness, or just report that lattice_points value is changed and use +1
                   default=lambda: DataFactory("core.int")(5),
                   help="Number of lattice points to be calculated for the equation of state. If even value supplied - value + 1 is used" #implement this
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

        # spec.expose_outputs(cls._step_2_workchain, namespace='relax', include=["relax", "misc", "retrieved"])
        spec.output("vol_opt_plot", valid_type=DataFactory('core.singlefile'), required = False, help='A plot visualizing the Volume-Total Energy relationship, with the fitted EoS')
        spec.output("vol_opt", valid_type=DataFactory('core.dict'), required = True, help='Minimum Volume-Energy pair.')

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
        self.ctx.workchain_metadata = namespace_to_atrdict(self.inputs.workchain_metadata, exclude = ["eos_metadata"])

        ## Step 1
        self.ctx.eos = AttributeDict()
        # change namespace to attribute dict, as the namespace is not a node in the aiida provenance graph, so this will not change provenance, and eliminate problems
        self.ctx.eos.structure_generation = namespace_to_atrdict(self.inputs.step_1.structure_generation)
        self.ctx.eos.inputs = AttributeDict()

        self.ctx.eos.inputs.workchain_metadata = namespace_to_atrdict(self.inputs.workchain_metadata.eos_metadata) #could I just give it a

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
        self.ctx.eos.inputs.update(self.exposed_inputs(self._step_1_workchain, namespace='step_1'))

        # Set up structures - #QUESTION - should I check if the parameters supplied contain a mode here again? if I have a validator that should be enough I guess
        self.report(f"Generating structures for {self._step_1_workchain.__name__} based on Structure {self.ctx.structure.pk} - Mode: {self.ctx.eos.structure_generation.parameters.get_dict().get('mode', 'uniform')}")
        self.ctx.eos.inputs.structures = generate_eos_structures(self.ctx.structure, **self.ctx.eos.structure_generation) #calcfunction can't deal with attributedicts or namespaces - have to give input nodes as dictionary, or directly

        self.ctx.eos.inputs = prepare_process_inputs(self.ctx.eos.inputs, namespaces=self.__class__._step_1_namespaces)
        # run step 1
        running_eos = self.submit(self._step_1_workchain, **self.ctx.eos.inputs)
        self.report(f"Launching {self._step_1_workchain.__name__}<{running_eos.pk}>")
        self.to_context(workchains=append_(running_eos))    # works, aiida vasp tutorial does it as well - I will have to check out how .run works later


    def verify_step_1(self):

        # Check if Step 1 WorkChain exists
        try:
            workchain = self.ctx.workchains[-1]
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
        self.ctx.relax.inputs.update(self.exposed_inputs(self._step_2_workchain, namespace='step_2'))       # Should I adjust here, the relax - yeah probably - TODO - add check for relax - or should that be a validator, no - add a check ,that changes so it is ISIF 3 and reports this change
        self.report(f"XXXXX\n{self.ctx.relax.inputs.settings.get_dict()}, {type(self.ctx.relax.inputs.settings)}")
        # Generate predicted Minimum Structure from Step 1
        eos_min_volume = np.round(self.ctx.eos.minimum.get('volume'), self.ctx.workchain_metadata.volume_precision)
        self.report(f"Generating predicted minimum structure for {self._step_2_workchain.__name__} based on inital Structure <{self.ctx.structure.pk}>, with Minimum Volume {eos_min_volume} \u213B\u00B3  from {self._step_1_workchain.__name__} <{self.ctx.workchains[0].pk}>")
        eos_relax_structure = self.ctx.structure.get_pymatgen().copy()

        eos_relax_structure.lattice = eos_relax_structure.lattice.matrix * (eos_min_volume / eos_relax_structure.volume)
        self.ctx.relax.inputs.structure = DataFactory("core.structure")(pymatgen_structure = eos_relax_structure)

        # Make sure that the CHGCAR is getting retrieved - TODO - Make optional
        self.report(f"overwriting parser settings for {", ".join([key for key in self._step_2_workchain_default_parser_settings.keys()])} - WIP!")
        relax_settings = self.ctx.relax.inputs.settings.get_dict()      # HERE TODO GO NEXT - I HAVE TO CORRECTLY ACCESS settings.parser_settings, change them and then put them back as full Dict() settings
        relax_settings['parser_settings'] |=  DataFactory("core.dict")(self._step_2_workchain_default_parser_settings)
        self.ctx.relax.inputs.settings = relax_settings

        self.report(f"XXXXX\n{self.ctx.relax.inputs.settings}")
        self.ctx.relax.inputs = prepare_process_inputs(self.ctx.relax.inputs, namespaces=self.__class__._step_2_namespaces)
        self.report(f"XXXXX\n{self.ctx.relax.inputs.settings}")
        # run step 2
        running_relax = self.submit(self._step_2_workchain, **self.ctx.relax.inputs)
        self.report(f"Launching {self._step_2_workchain.__name__}<{running_relax.pk}>")
        self.to_context(workchains=append_(running_relax))  # could replace assign with append and have self.ctx.relax.workchain = assign_(running_relax)

    def verify_step_2(self):

        try:
            workchain = self.ctx.workchains[-1]
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


        relaxed_pts = DataFactory("core.dict")(dict={
            "volume": self.ctx.workchains[1].outputs.relax.structure.get_pymatgen().volume,
            "energy": self.ctx.workchains[1].outputs.energies.get_array('energy_extrapolated')[0]
            })

        self.out('vol_opt', relaxed_pts)

        # for now I just want to expose some of the outputs of the relaxation workchain, and add the plot (opt)
        if self.ctx.workchain_metadata.create_plot:
            total_energies = self.ctx.workchains[0].outputs.eos
            parameters = self.ctx.workchains[0].outputs.eos_parameter

            self.report(f"Generating Plot of Equation of State Calculation <{self.ctx.workchains[0].pk}>, including relaxed structure <{self.ctx.workchains[1].outputs.relax.structure.pk}> from {self.ctx.workchains[1].__class__.__name__}<{self.ctx.workchains[1].pk}>")
            plot_file = create_eos_plot(total_energies  = total_energies,
                                        parameters      = parameters,
                                        relaxed_pts     = relaxed_pts)

            self.out('vol_opt_plot', plot_file)




@calcfunction
def generate_eos_structures(inp_structure, max_volume_scaling, lattice_points, parameters):
    structure = inp_structure.get_pymatgen()
    parameters = parameters.get_dict()
    mode = parameters.pop('mode')
    structure_name = parameters.get('str_name', 'Structure')

    max_vol_dev = max_volume_scaling.value
    num_lat_pts = lattice_points.value

    match mode:
        case 'uniform':
            volume_deviations = 1 + np.linspace(-max_vol_dev, max_vol_dev, num_lat_pts)
        case 'power':
            volume_deviations = 1 + powerrange(-max_vol_dev, max_vol_dev, num_lat_pts, power=parameters.get('power_constant', 1.5))


    structures = {}
    for vol_dev in volume_deviations:
        new_structure = structure.copy()
        new_structure.lattice = new_structure.lattice.matrix * vol_dev

        # IMPORTANT - Don't use '.' in the name - that will create sub-dictionaries .. as they are interpreted as namespaces
        structures[f"{structure_name}_at_{str(vol_dev).replace(".", "_")}_V0"] = DataFactory('core.structure')(pymatgen_structure = new_structure)

    return structures


def validate_structure_generation_parameters(para_node, port):
    parameter = para_node.get_dict()
    mode = parameter.get('mode', None)
    keys = parameter.keys()


    if mode in _NECESSARY_PARAMETERS.keys():
        if any([True for para in _NECESSARY_PARAMETERS[mode] if para not in keys]):
            nec_para = [para for para in _NECESSARY_PARAMETERS[mode] if para not in keys]
            return f"necessary parameters {nec_para} - were not supplied"

    elif mode is None:
        return "key 'type' is not provided"

    else:
        return f"value for key 'type' {mode} is invalid"
