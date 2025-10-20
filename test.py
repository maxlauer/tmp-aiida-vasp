import sys, importlib.util
print("python:", sys.executable)
print("find_spec:", importlib.util.find_spec('aiida_vasp_ext_mlauer'))