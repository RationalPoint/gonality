r"""
Sage script for looking for genus-4 canonical curves on one of our standard
quadric surface (xy + z^2, xy + zw, xy + N(z,w)) with a specified number of
rational points. This script is designed to locate *all of them*. It will 
be fastest if the number of rational points requested is close to the 
number of rational points on the surface {Q = 0}.

"""

# Modify only up to the next row of hashes
FF.<t> = GF(4)
R.<x,y,z,w> = FF[]
Q = x*y + z^2 + t*z*w + w^2
num_points_to_target = 14

################################################################################

import argparse
import itertools
import math
import progress
import search_tools
import sys


parser = argparse.ArgumentParser()
parser.add_argument('num_jobs',type=int,help='break up work into this many jobs')
parser.add_argument('job',type=int,help='which job is the current one?')
parser.add_argument('filename',type=str,help='where to write output')
args = parser.parse_args()

num_jobs = args.num_jobs
job = args.job
filename = args.filename


################################################################################

q = FF.cardinality()

# Space of cubic monomials. We can modify a cubic by a
# product of Q with a linear polynomials, which allows
# us to eliminate the x^2y, xy^2, z^3, and w^3 variables. 
monomials = search_tools.homogeneous_monomials(R,3)
if Q == x*y + z^2:
  kill = [x^2*y, x*y^2, z^3, z^2*w]
elif Q == x*y + z*w:
  kill = [x^2*y, x*y^2, z^2*w, z*w^2]
else:
  # We assume that Q = x*y + N(z,w), which always has a z^2 and w^2 term
  kill = [x^2*y, x*y^2, z^3, w^3]  
for term in kill:
  monomials.remove(term)
assert len(monomials) == 16

# Open the file before any serious computation is done
fp = open(filename,'w')

# Outer loop will be over choices of points on Q = 0
Q_pts = search_tools.proj_point_search(Q)
first_row = [Q.derivative(x) for x in R.gens()] # Jacobian matrix row

# Estimate total work this job will do.
num_pt_choices = binomial(len(Q_pts),num_points_to_target)
expected_dimension = len(monomials) - num_points_to_target
expected_number_of_pols = (q**expected_dimension - 1)//(q-1)
start,stop = search_tools.start_and_stop_work(num_pt_choices,num_jobs,job)
expected_work = (stop-start)*expected_number_of_pols
print('Targeting {} points: {} choices'.format(num_points_to_target,stop-start))
print('Searching vector spaces of expected dimension {}'.format(expected_dimension))
print('Total polynomials to search for this job: {} = 2^{:.2f}'.format(expected_work,math.log(expected_work,2)))

stats = {'dimension-1':0,'nonsingular':0,'irreducible':0}
prog = progress.Progress(num_pt_choices)
for i,pts in enumerate(itertools.combinations(Q_pts,num_points_to_target)):
  rows = []
  for pt in pts:
    rows.append([pol(pt) for pol in monomials])
  M = matrix(FF,rows)
  bb = M.right_kernel().basis()
  N = len(bb)
  num_pols = (q**N - 1)//(q-1)
  for cc in search_tools.proj_space_point_iterator(FF,N-1):
    coefs = sum(c*v for c,v in zip(cc,bb))
    F = sum(c*pol for c,pol in zip(coefs,monomials))
    num_pts = len(search_tools.proj_point_search(F,points=Q_pts))
    if num_pts != num_points_to_target: continue
    
    I = R.ideal(Q,F)
    # Check that Q = F = 0 defines a curve
    if search_tools.ideal_dimension(I) != 2: continue
    stats['dimension-1'] += 1
    # check for a singularity
    second_row = [F.derivative(x) for x in R.gens()]
    J = matrix(R,2,4,[first_row,second_row])
    sing_ideal = R.ideal([F,Q] + J.minors(2))
    if search_tools.ideal_dimension(sing_ideal) != 0: continue
    stats['nonsingular'] += 1    
    # check for irreducibility
    if len(I.primary_decomposition()) > 1: continue
    stats['irreducible'] += 1
    fp.write('{}\n'.format(F))
  prog()
prog.finalize()
print(stats)
