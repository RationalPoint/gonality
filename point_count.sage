r""" 

Sage script for computing quadric point counts.
This should be used in tandem with sage_launcher.py.

"""

# Modify only these next few lines
F.<t> = GF(2)
R.<v,w,x,y,z> = F[]
sing_dimension_file = '../data/genus5_GF2/gonality5_redo/proj_singularity_dimension.data'
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
kwargs['singularity_dimension_data'] = sing_dimension_file
kwargs['num_jobs'] = args.num_jobs
kwargs['job'] = args.job
kwargs['projective'] = True
kwargs['status'] = True

search_tools.quadric_point_count(R, **kwargs)



