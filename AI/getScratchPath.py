
# System-specific, adjust according to your system.
# For optimal performance, return the path to a very-fast storage space
# with ~200GB temporarily available space.

import os

def getScratchPath ():
    hostname = os.popen('hostname').readlines()[0].strip()
    print('Hostname: %s' % hostname)
    
    if 'slurm' in hostname or 'raven' in hostname:
        # username = os.getlogin()
        username = os.environ['USER']
        scratch_dir = '/hpc/scratch/%s/ai' % username
    else:
        scratch_dir = '/projects/covidtrakr/projects_scratch/ai'
    
    if not os.path.exists(scratch_dir):
        os.mkdir(scratch_dir)
    
    return scratch_dir
    
    