r""" 

Sage script for computing quadric singularity dimensions.
This should be used in tandem with sage_launcher.py.

"""

# Modify only these next few lines
F.<t> = GF(2)
R.<v,w,x,y,z> = F[]
steps_until_print = 2**17

################################################################################

import argparse
import search_tools

parser = argparse.ArgumentParser()
parser.add_argument('num_jobs',type=int,help='break up work into this many jobs')
parser.add_argument('job',type=int,help='which job is the current one?')
parser.add_argument('filename',type=str,help='where to write output')
args = parser.parse_args()

num_jobs = args.num_jobs
job = args.job
filename = args.filename

kwargs = {}
kwargs['steps_until_print'] = steps_until_print
kwargs['filename']  = args.filename
kwargs['num_jobs'] = args.num_jobs
kwargs['job'] = args.job
kwargs['projective'] = True
kwargs['status'] = True

search_tools.quadric_singularity_dimension(R, **kwargs)



