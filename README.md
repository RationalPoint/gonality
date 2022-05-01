### Search for Curves of Small Genus over Finite Fields



We give a description of the sequence of steps required to perform the computations in the papers

1. Faber, Xander and Grantham, Jon. "Binary Curves of small fixed genus and gonality with many rational points." _Journal of Algebra_ 597 (2022) pp. 24--26. [arXiv:2005.07054 [math.NT]](https://arxiv.org/abs/2005.07054)

2. Faber, Xander and Grantham, Jon. "Ternary and Quaternary Curves of small fixed genus and gonality with many rational points." _Experimental Mathematics_ (2021). [arXiv:2010.07992 [math.NT]](https://arxiv.org/abs/2010.07992)

CODE AUTHOR: Xander Faber

LICENSE: GPL-3.0 License

A command that begins with "sage:" is intended to be run in an interactive Sage session.  
A command that begins with "$" is intended to be run on the command line.

---

#### Procedure described in Lemma 5.6 of Reference 2.

These instructions assume that you have all of the Sage and Python scripts in a single working directory. Write Q = xy + z^2 + tzw + w^2 for these instructions.

1. Run the search for genus-4 canonical curves on V(Q) with 15 rational points.

   At the top of the `genus4_optimal_search.sage` script, adjust the value of the variable `num_points_to_target` to be 15. Then run the script with  
 
   ```
   $ sage genus4_optimal_search.sage 1 0 15pts.data
   ```  

   The file `15pts.data` will contain all cubic forms F such that V(F) \cap V(Q) is a smooth curve with 15 rational points.

2. Run the search for genus-4 canonical curves on V(Q) with 14 rational points.
 
   At the top of the `genus4_optimal_search.sage` script, adjust the value of the variable `num_points_to_target` to be 14. Then run the script with

   ```
    $ sage genus4_optimal_search.sage 1 0 14pts.data
   ```
   
   The file `14pts.data` will contain all cubic forms F such that V(F) \cap V(Q) is a smooth curve with 14 rational points.

---

#### Algorithms 2 and 3 in Reference 2.

These instructions assume that you have a working directory containing all of the Sage and Python scripts, as well as a subdirectory called `genus4` containing all of the C source files.

1. Run Algorithm 2 to find cubic forms that meet V(xy + N(z,w)) in a curve with no rational or quadratic point.

   Navigate to the genus4 subdirectory. Adjust the Makefile to point to your copies of the MPIR and FLINT libraries, which are part of your local Sage build. Compile and run the search over GF(4) with
  
   ```
   $ make search4
   ```
  
   This will create a file called `genus4_GF4.data`, each line of which gives a cubic form over GF(4). (Replacing `search4` with `search2` or `search3` will perform the corresponding searches over GF(2) or GF(3), respectively.)

   
2. Run Algorithm 3 to find smooth curves among the intersections of the cubic surfaces output by Algorithm 2 and the quadric surface V(xy + N(z,w))

   ```
   $ sage genus4_gonality5_finish.py 4 [path to genus4_GF4.data]
   ```
  
   (Replace all of the 4s with 2s or 3s if you would like to look at the case q = 2 or 3, respectively.)
  
   The output will be written to the file `genus4_GF4.data.smooth`.  It gives a complete list of isomorphism classes of cubic forms F that where identified such that V(F) \cap V(xy + N(z,w)) is a smooth curve of genus 4 and gonality 5. 

---

#### Procedures from Section 6.5 of Reference 2, including Algorithm 4

We describe the necessary commands for the computations over GF(3). If they differ substantially from the commands needed for the case of GF(2) or GF(4), we indicate these differences. These instructions assume that you have all of the Sage and Python scripts in a single working directory.

1. Compute orthogonal groups for the relevant quadratic forms.

   To write the projective orthogonal group of vw + x^2 + y^2 in GF(3)[v,w,x,y,z] to the file `POQ_vw_xx_yy.data`:
  
   ```
   sage: import search_tools
   sage: POQ = search_tools.odd_char_dim5_nonsplit_rank4_orthogonal_group(GF(3), 1)
   sage: search_tools.write_list_of_mats_to_file(POQ,'POQ_vw_xx_yy.data')
   ```
   
   To write the projective orthogonal group of vw + xy + z^2 in GF(3)[v,w,x,y,z] to the file `POQ_vw_xy_zz.data`:
  
   ```
   sage: import search_tools
   sage: POQ = search_tools.odd_char_dim5_smooth_orthogonal_group(GF(3))
   sage: search_tools.write_list_of_mats_to_file(POQ,'POQ_vw_xy_zz.data')
   ```
    
    To write the (projective) orthogonal group of vw + x^2 + txy + y^2 in GF(4)[v,w,x,y,z] to the file `POQ_vw_xx_txy+yy.data`:
  
   ```
   sage: import search_tools
   sage: FF.<t> = GF(4)
   sage: POQ = search_tools.char2_dim5_nonsplit_rank4_orthogonal_group(FF, t)
   sage: search_tools.write_list_of_mats_to_file(POQ,'POQ_vw_xx_txy_yy.data')
   ```
     
   To write the (projective) orthogonal group of vw + xy + z^2 in GF(4)[v,w,x,y,z] to the file `POQ_vw_xy_zz.data`:
  
   ```
   sage: import search_tools
   sage: FF.<t> = GF(4)     
   sage: POQ = search_tools.char2_dim5_smooth_orthogonal_group(FF)
   sage: search_tools.write_list_of_mats_to_file(OQ,'POQ_vw_xy_zz.data')
   ```
 
2. Compute ancillary information that determines the set B(Q1).

   First compute the singularity dimensions of each of the quadric surfaces over GF(3). We will use the `sage_launcher.py` script to launch and manage multiple instances of Sage. To see its syntax, run the script with no argument.
  
   Modify the top of the script `singularity_dimension.sage` in order to address quadratic forms over our field of interest. For GF(3), this will be:
  
   ```
   # Modify only these next few lines
   FF = GF(3)
   R.<v,w,x,y,z> = FF[]
   steps_until_print = 2**17
   ```
   
   Now launch all of the jobs with:
  
   ```
   $ python ./sage_launcher.py -o sing_dim.data singularity_dimension.sage sing_dim 24
   ```
  
   We need to concatenate the output files we just created in the correct order to ensure they will be used correctly in what follows:
  
   ```
   sage: import search_tools
   sage: filename = './sing_dim/sing_dim.data'
   sage: search_tools.cat_files(24,filename,filename)
   ```
  
   Next compute the number of rational points on each of the quadric surfaces over GF(3). You will need to modify the top of the script `point_count.sage` to say where you put the singularity dimension file `sing_dim.data`. It should read as follows for the computation we have described thus far:
  
   ```
   FF = GF(3)
   R.<v,w,x,y,z> = FF[]
   sing_dimension_file = 'sing_dim/sing_dim.data'
   steps_until_print = 2**17
   ``` 
  
   ```
   $ python ./sage_launcher.py -o point_count.data point_count.sage point_count 24
   ```
  
   Again, we need to concatenate the output files we just created in the correct order:
  
   ```
   sage: import search_tools
   sage: filename = './point_count/point_count.data'
   sage: search_tools.cat_files(24,filename,filename)
   ```
  
   Now we are getting into a part of the computation where it pays to keep things separate for our distinct choices of quadratic form Q1. Make a directory called `vw_xx_yy` and another directory called `vw_xy_zz`. (Change this appropriately for the quadratic forms involved in computing over GF(4).)
  
   Next create a file invariants.data in each of these subdirectories with the allowed values of (dimension of singular locus, number of rational points), one per line.  
   
   * For example, if Q1 = vw + x^2 + y^2 over GF(3), then `invariants.data` would look like  
     
     ```
     -1, 40  
     0, 31  
     ```
   
   * If Q1 = vw + xy + z^2 over GF(3), then `invariants.data` would look like
  
     ```
     -1, 40
     ```
   
   * If Q1 = vw + x^2 + txy + y^2 over GF(4), then `invariants.data` would look like
  
     ```
     -1, 85  
     0, 69
     ```
  
   * If Q1 = vw + xy + z^2 over GF(4), then `invariants.data` would look like
    
      ```
      -1, 85
      ```
  
   Now we combine this information into a bit mask for later use, which we write to the file `invariant_mask.data` in the appropriate subdirectory.
    
   ```
   sage: import search_tools
   sage: R.<v,w,x,y,z> = GF(3)[]  
   sage: sing_file = 'sing_dim/sing_dim.data'  
   sage: point_file = 'point_count/point_count.data'  
   sage: inv_data = 'vw_xx_yy/invariants.data'
   sage: mask_file = 'vw_xx_yy/invariant_mask.data'
   sage: args = (R,inv_data,mask_file,sing_file,point_file)
   sage: search_tools.allowed_invariant_quadric_mask(*args)
   ```
  
   Do the same for vw + xy + z^2 by replacing `vw_xx_yy` by `vw_xy_zz` in the above commands.

3. Compute the set A(Q1) -- orbit representatives for quadrics modulo the orthogonal group.

   ```
   sage: import search_tools
   sage: R.<v,w,x,y,z> = GF(3)[]
   sage: Q1 = v*w + x^2 + y^2
   sage: allowed_invariants = [(0,31),(-1,40)]
   sage: args = (Q1,allowed_invariants)
   sage: kwargs = {}
   sage: kwargs['quadric_sing_data'] = 'sing_dim/sing_dim.data'
   sage: kwargs['quadric_point_count_data'] = 'point_count/point_count.data'
   sage: kwargs['orthog_grp_data'] = 'POQ_vw_xx_yy.data'
   sage: kwargs['sage_outfile'] = 'vw_xx_yy/Q2s.data'
   sage: kwargs['steps_until_print'] = 2**15 # Progress printing
   sage: Q2s = search_tools.second_quadric_class_reps(*args,**kwargs)
  
   sage: Q1 = v*w + x*y + z^2
   sage: allowed_invariants = [(-1,40)]
   sage: args = (Q1,allowed_invariants)
   sage: kwargs = {}
   sage: kwargs['quadric_sing_data'] = 'sing_dim/sing_dim.data'
   sage: kwargs['quadric_point_count_data'] = 'point_count/point_count.data'
   sage: kwargs['orthog_grp_data'] = 'POQ_vw_xy_zz.data'
   sage: kwargs['sage_outfile'] = 'vw_xy_zz/Q2s.data'
   sage: kwargs['steps_until_print'] = 2**15 # Progress printing
   sage: Q2s = search_tools.second_quadric_class_reps(*args,**kwargs)
   ```

4. Find a smooth intersection of 3 quadrics with given number of points (and possibly gonality 6), or detect if no such intersection exists.

   Suppose that we want to target curves lying on Q1 = vw + x^2 + y^2 with 2 rational points. First, create a file `Q1.data` in your subdirectory `vw_xx_yy` with the single line
  
   ```
   v*w + x^2 + y^2
   ```
   
   Then, modify the top of the script `genus5_search.sage`.  Here is what it should look like in our case:
  
   ```
   # Modify only these next few lines
   FF = GF(3)
   R.<v,w,x,y,z> = FF[]
   num_points_to_target = 2
   check_for_gonality_six = False
  
   data_dir = './'
   quad_dir = 'vw_xx_yy/'
  
   q1_file = data_dir + quad_dir + 'Q1.data'
   q2_file = data_dir + quad_dir + 'Q2s.data'
   invariants_file = data_dir + quad_dir + 'invariants.data'
   invariant_mask_file = data_dir + quad_dir + 'invariant_mask.data'
  
   steps_until_print = 2**20
   ```
  
   Now we run the script with the command
  
   ```
   $ python ./sage_launcher.py -c -o vw_xx_yy.data genus5_search.sage vw_xx_yy/curves.2 48
   ```
  
   Each line of the file `vw_xx_yy/curves.2/vw_xx_yy.data` is of the form
  
   ```
   Q1, Q2, Q3
   ```
  
   where the variety V(Q1,Q2,Q3) is a smooth genus-5 curve with gonality 5 and 2 rational points. (There may be more than one line in the file because multiple Sage instances found a curve at almost the same time.) If the file is empty, then there is no smooth genus-5 curve with gonality 5 and 2 rational points lying on V(vw + x^2 + y^2). 
  
   To look for curves lying on Q1 = vw + xy + z^2 with gonality 6, we start by creating the file `vw_xy_zz/Q1.data` with the single line
  
   ```
   v*w + x*y + z^2
   ```
  
   Then we modify the top of `genus5_search.sage` as follows:
  
   ```
   # Modify only these next few lines
   FF = GF(3)
   R.<v,w,x,y,z> = FF[]
   num_points_to_target = 0
   check_for_gonality_six = True
  
   data_dir = './'
   quad_dir = 'vw_xy_zz/'
  
   q1_file = data_dir + quad_dir + 'Q1.data'
   q2_file = data_dir + quad_dir + 'Q2s.data'
   invariants_file = data_dir + quad_dir + 'invariants.data'
   invariant_mask_file = data_dir + quad_dir + 'invariant_mask.data'
  
   steps_until_print = 2**20
   ```
  
   Now run the script with the command
  
   ```
   $ ./sage_launcher.py -c -o vw_xy_zz.data genus5_search.sage vw_xy_zz/curves.gonality6 48
   ```
  
   Each line of the file `vw_xx_yy/curves.gonality6/vw_xy_zz.data` is of the form
  
   ```
   Q1, Q2, Q3
   ```
   
   where the variety V(Q1,Q2,Q3) is a smooth genus-5 curve with gonality 6 (and no rational point). This file was empty when we ran the computation, which means there is no smooth genus-5 curve with gonality 6 lying on V(vw + xy + z^2).
