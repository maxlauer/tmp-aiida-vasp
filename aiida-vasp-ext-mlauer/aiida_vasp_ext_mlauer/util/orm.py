
from aiida.common.extendeddicts import AttributeDict
from plumpy.utils import AttributesFrozendict

"""
## TODO 
# namespace_to_atrdict
 - maybe keep option to append namespaces, as namespaces, and not as Attribute Dict
 - maybe add option to rename subnamespaces
"""

def namespace_to_atrdict(namespace, exclude=[]):
    atr_dict = AttributeDict()
    
    for key, val in namespace.items():
        
        if key not in exclude:
            
            if isinstance(val, AttributesFrozendict):
                sub_exclude = [exc.strip(f"{key}.") for exc in exclude if exc.split('.')[0] == key]
                atr_dict[key] = namespace_to_atrdict(val, exclude=sub_exclude)

            else:
                atr_dict[key] = val 

    return atr_dict
