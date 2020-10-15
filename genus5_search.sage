r""" 

Sage script for looking for genus-5 curves with a specified number of points or
with no g15. Each instance either finds a curve and exits with code 1 or fails
to find a curve and exits with code 0. As a result, we can either find a curve
as desired or certify that none exists.

To run this script, we must have the following files already:

q1_file             = file containing string representation of Q1
q2_file             = file containing all choices of Q2
invariants_file     = file containing allowed invariants for the search
invariant_mask_file = file containing the bitmask for an array of allowed quadrics 

Modify the lines below to give the desired target information. The variable
num_points_to_target can be an integer or a list of integers. If
check_for_gonality_six is True, then num_points_to_target should be set to 0
(as only pointless curves can have gonality 6).

This should be used in tandem with sage_launcher.py.

"""

# Modify only these next few lines
F.<t> = GF(4)
R.<v,w,x,y,z> = F[]
num_points_to_target = range(5,18)
check_for_gonality_six = False

data_dir = '../data/genus5_GF4/'
quad_dir = 'vw_xx_txy_yy/'

q1_file = data_dir + quad_dir + 'Q1.data'
q2_file = data_dir + quad_dir + 'Q2s.data'
invariants_file = data_dir + quad_dir + 'invariants.data'
invariant_mask_file = data_dir + quad_dir + 'invariant_mask.data'

steps_until_print = 2**20


################################################################################

import argparse
import itertools
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

if len(R.gens()) != 5:
  raise ValueError('Polynomialring does not have 5 generators')
FF = R.base_ring()
if not FF.is_finite():
  raise ValueError('Base ring is not finite')
q = FF.cardinality()

# open output file
gp = open(filename,'w')  

# get input from files
fp = open(q1_file,'r')
Q1 = R(fp.readline())
fp.close()

Q2s = []
fp = open(q2_file,'r')
for Q in fp:
  Q2s.append(R(Q))
fp.close()
print('Computing nontrigonal genus-5 curves lying on {} = 0'.format(Q1))
if isinstance(num_points_to_target,(list,tuple)):
  num_pts_str = ''
  m = len(num_points_to_target)
  for i in xrange(m-1):
    num_pts_str += '{}, '.format(num_points_to_target[i])
  num_pts_str += 'or {}'.format(num_points_to_target[m-1])
else:
  num_pts_str = str(num_points_to_target)
if check_for_gonality_six:
  if num_points_to_target != 0:
    raise ValueError('There is no curve with gonality six and {} rational points'.format(num_pts_str))
  print('Targeting curves with gonality six.')
else:
  print('Targeting curves with {} points.'.format(num_pts_str))  
print('Using {} forms for Q2.'.format(len(Q2s)))

allowed_invariants = set()
fp = open(invariants_file,'r')
for line in fp:
  a,b = [int(t) for t in line.strip().split(',')]
  allowed_invariants.add((a,b))
fp.close()
print('Allowing invariants {} for quadratic_forms.'.format(list(allowed_invariants)))

q1_pts = search_tools.proj_point_search(Q1)
if check_for_gonality_six:
  q1_cubic_pts = set()
  seen = set(q1_pts)
  S = PolynomialRing(FF,names=('T',))
  while True:
    minpol = S.random_element(3)
    if minpol.is_irreducible(): break
  E = FF.extension(minpol,name='T')
  QE = Q1.change_ring(E)
  kwds = {'status':True,'steps_until_print':steps_until_print}
  for pt in search_tools.hypersurface_point_search(QE,**kwds):
    if pt in seen: continue
    q1_cubic_pts.add(pt)
    for i in range(3):
      seen.add(tuple(u**(q**i) for u in pt))
  print('There are {} orbits of cubic points to test.'.format(len(q1_cubic_pts)))

# Construct an array of correct length and
# set the bit for each quadric with allowed invariants
quads = search_tools.homogeneous_monomials(R,2)
num_pols = (q**len(quads) - 1)//(q-1)
allowed_polys = [0]*num_pols # 1 means allowed
num_allowed = 0
fp = open(invariant_mask_file,'r')
s = fp.readline()
fp.close()
print('Constructing array of allowed quadrics.')
prog = progress.Progress(num_pols,steps_until_print=steps_until_print)
for i in xrange(num_pols):
  prog()
  val = int(s[i])
  if val:
    allowed_polys[i] = val
    num_allowed += 1
prog.finalize()

num_megs = float(sys.getsizeof(allowed_polys)*1.0 / 2**20)
msg = 'Found {} allowed quadrics ({:.2f}MB)'
print(msg.format(num_allowed,num_megs))
del s
sys.stdout.flush()

# Get other invariant data
# sdp = open(singularity_file,'r')
# pcp = open(count_file,'r')

# # Construct an array of correct length and
# # set the bit for each quadric with allowed invariants
# old_allowed_polys = [0]*num_pols # 1 means allowed
# old_num_allowed = 0
# print('Constructing array of allowed quadrics the old way.')
# prog = progress.Progress(num_pols,steps_until_print=steps_until_print)
# for i in xrange(num_pols):
#   prog()
#   a,b = int(sdp.readline()), int(pcp.readline())
#   if (a,b) in allowed_invariants:
#     old_allowed_polys[i] = 1
#     old_num_allowed += 1
# prog.finalize()
# from sys import getsizeof
# num_megs = float(getsizeof(old_allowed_polys)*1.0 / 2**20)
# msg = 'Found {} allowed quadrics ({:.2f}MB)'
# print(msg.format(old_num_allowed,num_megs))

# # Clean up files
# sdp.close()
# pcp.close()

# assert allowed_polys == old_allowed_polys


start,stop = search_tools.start_and_stop_work(num_pols, num_jobs, job)
if num_jobs > 1:
  print('Start index = {}, Stop index = {}'.format(start,stop))  

num_pairs = len(Q2s) * (stop-start)
print('Searching {} forms for Q3.'.format(stop-start))
print('Total pairs to search: {}\n'.format(num_pairs))

# Loop over Q2,Q3 pairs and test
AA2_pts = list(search_tools.affine_space_point_iterator(FF,2))
prog = progress.Progress(num_pairs,steps_until_print=steps_until_print)
if isinstance(num_points_to_target,(list,set,tuple)):
  min_points_to_target = min(num_points_to_target)
else:
  min_points_to_target = num_points_to_target
found_one = False
for Q2 in Q2s:
  if found_one: break
  if len(search_tools.proj_point_search(Q2,points=q1_pts)) < min_points_to_target:
    prog(stop-start)
    continue
  # Compute cubic points on Q1 = Q2 = 0 if looking for gonality 6
  if check_for_gonality_six:
    q12_cubic_pts = search_tools.proj_point_search(Q2,points=q1_cubic_pts)
  for Q3 in search_tools.poly_iterator(quads,projective=True,start=start,stop=stop):
    prog()
    Qs = [Q1,Q2,Q3]
    good_triple = True
    for aa in AA2_pts:
      Q = Q3 + aa[0]*Q1 + aa[1]*Q2
      if Q == 0: continue
      j = search_tools.poly_index(Q,quads,projective=True)
      if not allowed_polys[j]:
        good_triple = False
        break
    if not good_triple: continue
    # Now the invariants for all Q's in the span are ok.

    # Test number of points
    right_num_pts = search_tools.proj_variety_has_n_points([Q2,Q3],num_points_to_target,points=q1_pts)
    if not right_num_pts: continue

    # Check if this curve has gonality 6    
    if check_for_gonality_six:
      # If the curve has a rational point, then it has a g_1^5.
      if num_points_to_target != 0: continue
      # If the curve has a cubic point, then it has a g_1^5. 
      has_cubic_pt = not search_tools.proj_variety_has_n_points([Q3],0,points=q12_cubic_pts)
      if has_cubic_pt: continue

    # Test irreducibility and smoothness
    I = R.ideal(Qs)
    if search_tools.ideal_dimension(I) != 2: continue

    # check for a singularity other than the affine cone point
    inds = itertools.product(xrange(3),xrange(5))
    D = {(i,j): Qs[i].derivative(R.gens()[j]) for i,j in inds}
    J = matrix(R,3,5,D)
    sing_ideal = R.ideal(Qs+J.minors(3))
    if search_tools.ideal_dimension(sing_ideal) != 0: continue
    # Now V(I) is a smooth curve, given as the intersection of 3 quadrics.
    # The only thing left that could fail is that C is not irreducible.
    if not search_tools.ideal_is_prime(I): continue
    gp.write('{}, {}, {}\n'.format(*Qs))
    found_one = True
    break
prog.finalize()
gp.close()

if found_one:
  # This tells sage_launcher.py to kill all the other jobs and quit.
  if check_for_gonality_six:
    print('Found a smooth (pointless) curve of genus 5 with gonality 6.')
  else:
    print('Found a smooth curve of genus 5 with {} points'.format(num_pts_str))
  sys.exit(int(1))
else:
  print('Found no curve satisfying desired conditions.')
  sys.exit(int(0))
  
