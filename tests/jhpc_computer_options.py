from aiida.orm import Dict



def get_jhpc_options(num_machines, queue=None):
    if queue == 'debug':
        options = Dict(dict={'resources': {'num_machines': num_machines, },
                             'max_wallclock_seconds': 1800,
                             'queue_name': 'debug',
                             })
    elif num_machines == 1:  # single partition
        options = Dict(dict={'resources': {'num_machines': num_machines, },
                             'max_wallclock_seconds': 432000,
                             'queue_name': 'single',
                             # 'prepend_text': PREPEND_TEXT_SINGLE,
                             # 'append_text': APPEND_TEXT_SINGLE,
                             })
    else:
        options = Dict(dict={'resources': {'num_machines': num_machines, },
                             'max_wallclock_seconds': 432000,
                             'queue_name': 'regular',
                             # 'custom_scheduler_commands': '#SBATCH --exclude=node1-[001-040]',
                             })
        
    if queue in ['short', 'long', 'bigmem_long']:
        options.update_dict({'qos': queue,
                             'max_wallclock_seconds': 43200 if queue == 'short' else 4*864000})
        # wallclock for long is in reality unlimited on jhpc, can be raised
    else:
        options.update_dict({'qos': 'normal'})
    
    return options

