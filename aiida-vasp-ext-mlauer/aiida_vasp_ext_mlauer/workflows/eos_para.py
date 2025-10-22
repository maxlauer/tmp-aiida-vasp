import numpy as np
import io 

from scipy.interpolate import interp1d
from scipy.optimize import minimize, curve_fit
from matplotlib import pyplot as plt

from aiida.engine import WorkChain, calcfunction, while_, append_
from aiida.common.extendeddicts import AttributeDict
from aiida.plugins import DataFactory, WorkflowFactory

from aiida_vasp.utils.workchains import prepare_process_inputs # Should check to see that the inputs given to vasp.vasp works well .. not too sure
from aiida_vasp.utils.workchains import compose_exit_code # combines exit codes together

from aiida.orm import Dict

from aiida_vasp_ext_mlauer.util.functions import Murnaghan

_DEFAULT_K0 = 10.0
_DEFAULT_K1 = 4.5



class EoSParaWorkChain(WorkChain):
    """
    The Parallel Equation of State Work Chain 
    ----------------------------------

    The WorkChain accepts a dictionary of structures and extracteds.
    Compared to the normal EoSWorkChain, this WorkChain runs all calculations in parallel.
    """

    _verbose = False 
    _next_workchain_string = "vasp.vasp"
    _next_workchain = WorkflowFactory(_next_workchain_string)

    
    @classmethod
    def define(cls, spec):
        super(EoSParaWorkChain, cls).define(spec)
        spec.expose_inputs(cls._next_workchain, exclude=["structure"])
        spec.input_namespace(name = "structures", 
                             valid_type = DataFactory("core.structure"),
                             dynamic = True,
                             help = "Dictionary of Structures used for the Equation of State")
        # spec.input('structure', 
        #            valid_type = DataFactory("core.structure",
        #            dynamic=True,
        #            help ="Base Structure for which to perform the calculation")
        #         )
        # spec.input('init_alat', 
        #            valid_type = DataFactory("core.float"), 
        #            dynamic = True, 
        #            default = None, 
        #            optional = True, 
        #            help = "Lattice Constant, around which to perform the scan. If not given, use structure lattice constant"
        #         )
        # spec.input('deviation', 
        #            valid_type = DataFactory("core.float"), 
        #            dynamic=True, 
        #            help = "Amount by which to vary the lattice constant. If > 1 considered as %"
        #         )

        spec.input('minimum_mode',
                   valid_type = DataFactory("core.str"),
                   default = lambda: DataFactory("core.str")('Interpolate'), 
                   help = "Determines how the energy minimum is evaluated - Interpolate: 1D cubic Interpolation; Murnaghan: fit Murnaghan EoS"
                )
        

        spec.input_namespace(name='wc_metadata',
                             dynamic = False,
                             help='Custom namespace for WorkChain-specific metadata (e.g. verbosity, plotting, tags).'
                )
        
        spec.input('wc_metadata.create_plot',
                   valid_type = bool,
                   default = True,
                   help = "Toggle whether to create a plot, or not. Default is False",
                   non_db = True
                )
        

        # Specify Exit Codes
        spec.exit_code(0, "NO_ERROR", message="Calculation ended without Error ... well done I guess")
        spec.exit_code(420, "ERROR_NO_CALLED_WORKCHAIN", message="No called WorkChain detected")
        spec.exit_code(500, "ERROR_UNKOWN", message="Unkown Error detected in the EOS WorkChain")

        spec.outline(
            cls.initialize,
            cls.init_and_run_next_workchains,
            cls.verify_next_workchain,
            cls.extract_volume_and_energy,
            cls.finalize
        )

        spec.output('eos', valid_type=DataFactory('core.array'), help='A List of Cell Volumes and Total Energies')
        spec.output('eos_plot', valid_type=DataFactory('core.singlefile'), required = False, help='A plot visualizing the Volume-Total Energy relationship, with the fitted EoS')
        spec.output('eos_parameter', valid_type=DataFactory('core.dict'), help="A Dict containing the all parameters and the 'type' of the fitted EoS")
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

        self.ctx.iteration = 0

        self.ctx.mode = self.inputs.minimum_mode

        self.ctx.wc_metadata = self.inputs.wc_metadata

        # define context inputs
        self.ctx.inputs = AttributeDict()

        self.ctx.total_energies = []

    
    def _init_inputs(self):
        
        try:
            self._verbose = self.inputs.verbose.value
        
        except AttributeError:
            pass 



    def init_and_run_next_workchains(self):


        # Check existence of context inputs
        try:
            self.ctx.inputs
        except AttributeError:
            raise ValueError("No input dictionary was defined in self.ctx.inputs")
        

        # Take exposed inputs (inputs passed to the WorkChain) and add them to context
        self.ctx.inputs.update(self.exposed_inputs(self._next_workchain))

        # Only missing is  `structure`, which we didn't expose
        # choose randomly from `structures` entry point 
        for structure_key, structure in self.ctx.structures.items():
            print(structure_key, structure)

            # Increase iteration index, update structure of ctx.inputs and make sure the inputs are correctly formatted
            self.ctx.iteration += 1
            self.ctx.inputs.structure = structure
            self.ctx.inputs = prepare_process_inputs(self.ctx.inputs, namespaces=['dynamics'])

            # Submit VaspWorkChain for Structure
            running = self.submit(self._next_workchain, **self.ctx.inputs)

            self.report(f'launching {self._next_workchain.__name__}<{running.pk}> iteration #{self.ctx.iteration}')
            self.to_context(workchains=append_(running))


    def verify_next_workchain(self):

        try:
            workchain = self.ctx.workchains[-1]
        except IndexError:
            self.report(f"There is no {self._next_workchain.__name__} in the called workchain list.")
            return self.exit_codes.ERROR_NO_CALLED_WORKCHAIN
        
        for workchain in self.ctx.workchains:
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
        """Extract the cell volume and total energy for this structure."""

        for workchain in self.ctx.workchains:
            # Fetch the total energy
            misc = workchain.outputs.misc.get_dict()
            total_energy = misc['total_energies']['energy_extrapolated']

            # Fetch the volume
            self.report(f"{workchain.inputs}")
            volume = workchain.inputs.structure.get_cell_volume()

            # Store both in a list
            self.ctx.total_energies.append([volume, total_energy])

    
    def finalize(self):
        # To get DataProvenance, the AiiDA List container has to be prepared by a calcfunction
        total_energies = store_total_energies(DataFactory('list')(list=self.ctx.total_energies))


        if self.ctx.mode == "Interpolate":
            analysis_dict = locate_minimum_interpolate(total_energies)
            self.report(f"Minimum Energy Determination using 1D cubic Interpolation")
            self.out('eos_parameter', DataFactory('core.dict')({'type': "cubic 1D Interpolation"}))

        elif self.ctx.mode == "Murnaghan":
            analysis_dict = locate_minimum_murnaghan(total_energies)
            self.report(f"Minimum Energy Determination fitting the Murnaghan Equation of State")
            

        self.out('eos', total_energies)
        self.out('eos_minimum', analysis_dict['min_energy'])
        self.out('eos_parameter', analysis_dict['parameters'])



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
def locate_minimum_interpolate(total_energies):
    total_energies_array = total_energies.get_array('eos')

    volumes = total_energies_array[:,0]
    energies = total_energies_array[:,1]

    ## get minimum energy guess
    min_energy_guess_idx = energies.argmin()

    new_energies = interp1d(volumes, energies, kind='cubic')
    min_energy_point = minimize(new_energies, volumes[min_energy_guess_idx], tol=1e-3)

    # Store Result
    dict_data = {'volume': min_energy_point.x[0], 'energy': min_energy_point.fun}

    parameters = {'type': '1D cubic Interpolation'}

    return {
        "min_energy": DataFactory("core.str")("A"), # Dict(dict = dict_data),
        "parameters": DataFactory("core.str")("A") # Dict(dict = parameters)
    }



@calcfunction
def locate_minimum_murnaghan(total_energies):
    total_energies_array = total_energies.get_array('eos')

    volumes = total_energies_array[:,0]
    energies = total_energies_array[:,1]

    # Document Type of Analysis
    parameters = {'type': 'Murnaghan EoS Fit'}

    ## get minimum energy guess
    min_energy_guess_idx = energies.argmin()

    ## Fit Murnaghan Equation of State
    params, covar = curve_fit(Murnaghan, volumes, energies, p0=[volumes[min_energy_guess_idx], energies[min_energy_guess_idx], _DEFAULT_K0, _DEFAULT_K1])
    parameters |= {'E0': params[0], 
                   'V0': params[1], 
                   'K0': params[2], 
                   'K1': params[3],
                   'covariance': covar
                   }

    min_energy_point = Murnaghan(params[1], *params)

    # Store Result
    dict_data = {'volume': params[0], 'energy': min_energy_point}

    return {
        "min_energy": DataFactory("core.dict")(dict = dict_data), # Dict(dict = dict_data),
        "parameters": DataFactory("core.dict")(dict = {"A": "B"}) # Dict(dict = parameters)
    }


    # if create_plot:
    #     fig, ax = plt.subplots()
    #     ax.set_title('Murnaghan Plot')

    #     ax.scatter(volumes, energies)

    #     x_fit = np.linspace(volumes.min() - 1, volumes.max() + 1, 500)    # I should make these numbers optional .. later
    #     y_fit = Murnaghan(x_fit, *params)
    #     ax.plot(x_fit, y_fit, label='Murnaghan Fit', color='orange')

    #     ax.axvline(x=params[1], color='grey', linestyle='--', label=fr'$V_0 = {params[1]:.3f}$')

    #     ax.set_xlabel(r'Volume / $10^\mathrm{-30}$ m^\mathrm{3}')
    #     ax.set_ylabel(r'Energy / eV')
    #     ax.legend()

    #     buffer = io.BytesIO()
    #     fig.savefig(buffer, format='png')
    #     plt.close(fig)
    #     buffer.seek(0)

    #     plot_file = DataFactory('core.singlefile')(file=buffer, filename='murnaghan_plot.png')
    # else:
    #     plot_file = None

    # return dict_data, params, covar, plot_file



# @calcfunction
# def plot_EoS(total_energies, params, type):
#     total_energies_array = total_energies.get_array('eos')

#     volumes = total_energies_array[:,0]
#     energies = total_energies_array[:,1]


