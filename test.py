from aiida.plugins import WorkflowFactory
from aiida.engine import ProcessBuilder

builder = ProcessBuilder(WorkflowFactory("vasp.vasp"))
builder_relax = ProcessBuilder(WorkflowFactory("vasp.relax"))
builder_eos = ProcessBuilder(WorkflowFactory("vasp_ext_ml.eos"))

print(builder)
print(builder_relax)
print(builder_eos)