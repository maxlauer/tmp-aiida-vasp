from aiida_vasp_ext_mlauer.workflows.eos import EoSRelaxWorkChain

from aiida.common.extendeddicts import AttributeDict
from aiida import load_profile

load_profile('lauerm-test')

builder = EoSRelaxWorkChain.get_builder()
print(builder)