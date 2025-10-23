import numpy as np
import random

from aiida.engine import WorkChain, calcfunction, while_, append_
from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import DataFactory, WorkflowFactory

from aiida_vasp.utils.workchains import prepare_process_inputs # Should check to see that the inputs given to vasp.vasp works well .. not too sure
from aiida_vasp.utils.workchains import compose_exit_code # combines exit codes together

from scipy.interpolate import interp1d
from scipy.optimize import minimize

## THIS IS A STUB .. DOESN"T WORK ...

class EoSWorkChain_old(WorkChain):
    """
    The Equation of State Work Chain
    ----------------------------------

    The WorkChain accepts a dictionary of structures and extracteds
    """

    _verbose = False 
    _next_workchain_string = "vasp.vasp"
    _next_workchain = WorkflowFactory(_next_workchain_string)

    
    @classmethod
    def define(cls, spec):
        super(EoSWorkChain_old, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=["structure"])
        spec.input_namespace(name = "structures", 
                             valid_type = DataFactory("core.structure"),
                             dynamic = True,
                             help = "Dictionary of Structures used for the Equation of State")
        
        # Specify Exit Codes
        spec.exit_code(0, "NO_ERROR", message="Calculation ended without Error ... well done I guess")
        spec.exit_code(420, "ERROR_NO_CALLED_WORKCHAIN", message="No called WorkChain detected")
        spec.exit_code(500, "ERROR_UNKOWN", message="Unkown Error detected in the EOS WorkChain")

        spec.outline(
            cls.initialize,
            while_(cls.run_next_workchains)(
                cls.init_next_workchain,
                cls.run_next_workchain,
                cls.verify_next_workchain,
                cls.extract_volume_and_energy
            ),
            cls.finalize
        )

        spec.output('eos', valid_type=DataFactory('core.array'), help='A List of Cell Volumes and Total Energies')
        # spec.output('eos_plot', valid_type=DataFactory('core.singlefile'), help='A plot visualizing the Volume-Total Energy relationship, with the fitted EoS')
        # spec.output('eos_parameter', valid_type=DataFactory('core.dict'), help="A Dict containing the all parameters and the 'type' of the fitted EoS")
        spec.output('eos_minimum', valid_type=DataFactory('core.dict'), help="A cell volume containing the cell volume and total energy at energy minimum")


    def initialize(self):
        """ """
        self._init_context()
        self._init_inputs()


    def _init_context(self):
        # set exit code to error - I don't understand quite why yet .. but ok
        self.ctx.exit_code = self.exit_codes.ERROR_UNKOWN

        # Load Structures in Context
        self.ctx.structures = dict(self.inputs.structures)

        self.ctx.is_finished = False

        self.ctx.iteration = 0

        # define context inputs
        self.ctx.inputs = AttributeDict()

        self.ctx.total_energies = []

    
    def _init_inputs(self):
        
        try:
            self._verbose = self.inputs.verbose.value
        
        except AttributeError:
            pass 



    def run_next_workchains(self):
        """
        Return whether a new WorkChain should be run
        
        Returns True, if the last WorkChain has not finished successfully 
        """
        return not self.ctx.is_finished


    def init_next_workchain(self):

        self.ctx.iteration += 1

        # Check existence of context inputs
        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError("No input dictionary was defined in self.ctx.inputs")
        

        # Take exposed inputs (inputs passed to the WorkChain) and add them to context
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Only missing is  `structure`, which we didn't expose
        # choose randomly from `structures` entry point 
        item = random.choice(list(self.inputs.structures.keys()))
        self.ctx.inputs.structure = self.ctx.structures.pop(item)

        self.ctx.inputs = prepare_process_inputs(self.ctx.inputs, namespaces=['dynamics'])

    
    def run_next_workchain(self):
        inputs = self.ctx.inputs
        running = self.submit(self._next_workchain, **inputs)
        self.report(f'launching {self._next_workchain.__name__}<{running.pk}> iteration #{self.ctx.iteration}')
        self.to_context(workchains=append_(running))


    def verify_next_workchain(self):

        try:
            workchain = self.ctx.workchains[-1]
        except IndexError:
            self.report(f"There is no {self._next_workchain.__name__} in the called workchain list.")
            return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN
        
        # Inherit exit status from last workchain
        next_workchain_exit_status = workchain.exit_status
        next_workchain_exit_message = workchain.exit_message

        # if last workchain successful
        if not next_workchain_exit_status:
            self.ctx.exit_code = self.exit_codes.NO_ERROR

        else:
            self.ctx.exit_code = compose_exit_code(next_workchain_exit_status, next_workchain_exit_message)
            self.report(f"The called{workchain.__class__.__name__}<{workchain.pk}> returned a non-zero exit status.\nThe exit status {self.ctx.exit_code} is inherited.")

        
        if not self.ctx.structures:
            self.ctx.is_finished = True 

        return self.ctx.exit_code
    


    def extract_volume_and_energy(self):
        """
        Extract cell volume and total energy from the respective structure
        """

        workchain = self.ctx.workchains[-1]
        misc = workchain.outputs.misc.get_dict()
        total_energy = misc['total_energies']['energy_extrapolated']

        # Get Volume
        volume = self.ctx.inputs.structure.get_cell_volume()

        # Store
        self.ctx.total_energies.append([volume, total_energy])

    
    def finalize(self):
        # To get DataProvenance, the AiiDA List container has to be prepared by a calcfunction
        total_energies = store_total_energies(DataFactory('list')(list=self.ctx.total_energies))

        energy = locate_minimum(total_energies)

        self.out('eos', total_energies)
        self.out('eos_minimum', energy)



@calcfunction
def store_total_energies(total_energies):
    total_energies_list = total_energies.get_list()

    # sort by volume
    total_energies_array = np.array(total_energies_list)
    total_energies_array_sorted = total_energies_array[total_energies_array[:, 0].argsort()]

    array_data = DataFactory('array')()
    array_data.set_array('eos', total_energies_array_sorted)

    return array_data


@calcfunction
def locate_minimum(total_energies):
    total_energies_array = total_energies.get_array('eos')

    volumes = total_energies_array[:,0]
    energies = total_energies_array[:,1]

    ## get minimum energy guess
    min_energy_guess_idx = energies.argmin()


    # for now cubic interpolation - later on implement murnaghan equation fitting, or other EOS or other interpolation schemes
    new_energies = interp1d(volumes, energies, kind='cubic')
    min_energy_point = minimize(new_energies, volumes[min_energy_guess_idx], tol=1e-3)

    # Store Result
    dict_data = DataFactory('dict')(dict = {'volume': min_energy_point.x[0], 'energy': min_energy_point.fun})

    return dict_data 