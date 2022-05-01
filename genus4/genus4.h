/*----------------------------------------------------------------------------*/
/*		      Genus-4 Excessive Curve Search Code                     */
/*----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

An excessive curve is a curve over a finite field of genus g and gonality g+1.
We only search for excessive curves when q = 2,3,4 because Lauter's method rules
out the case q = 5, 7, and Weil's inequality rules out q > 7.

Every such curve is canonically embedded in P^3 as the intersection
of a unique quadric surface and a cubic surface, and the quadric
surface must be smooth and nonsplit over the ground field. 
q = 2: Q = x*y + z^2 + z*w + w^2
q = 3: Q = x*y + z^2 + w^2
q = 4: Q = x*y + z^2 + t*z*w + w^2, where t^2 + t + 1 = 0

TABLE OF CONTENTS

- MACROS

  - MEMORY_ERROR

- FORMAT_TIME (format_time.c)

  - format_time

- FIELDS (fields.c)

  - field_elements_array
  - subfield_elements_array
  - construct_fields
  - construct_add_table
  - fq_add_table
  - construct_mul_table
  - fq_mul_table
  - fq_vec_dot_table

- POINTS (points.c)

  - point_t
  - point_init
  - point_clear
  - point_set
  - point_pretty_print
  - point_normalize
  - point_eq_raw
  - point_eq
  - point_frobenius
  - point_quadratic_frobenius
  - next_vector
  - next_projective_vector
  - test_next_projective_vector

- CURVES (curves.c)

  - cubic_coef_vec_t (struct)
  - cubic_coef_vec_init
  - cubic_coef_vec_clear
  - cubic_coef_vec_pretty_print
  - cubic_coef_vec_eval
  - cubic_monomial_evaluations
  - print_cubic_to_file_old
  - print_cubic_to_file
  - print_cubic_to_screen
  - nonsplit_quadric_surface_points
  - cubic_coef_vec_to_index_vec

------------------------------------------------------------------------------*/

#ifndef GENUS4_H
#define GENUS4_H

#include <fmpz.h>
#include <fq.h>
#include <fq_vec.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*----------------------------------------------------------------------------*/
/*---------------------------------  MACROS  ---------------------------------*/
/*----------------------------------------------------------------------------*/

#define MEMORY_ERROR {{fprintf(stdout,"ERROR: Cannot allocate memory at line %d\n",__LINE__); exit(1);}}

/*----------------------------------------------------------------------------*/
/*-------------------------------  FORMAT TIME  ------------------------------*/
/*----------------------------------------------------------------------------*/

// Format time t as hrs,mins,secs and write to character array time_str
void format_time(char *time_str, time_t t);

/*----------------------------------------------------------------------------*/
/*---------------------------------  POINTS  ---------------------------------*/
/*----------------------------------------------------------------------------*/

// Point of P^3(GF(q))
typedef struct {
  fq_t x;
  fq_t y;
  fq_t z;
  fq_t w;  
} point_t;

// Initialize and clear P in PP^3(F)
void point_init(point_t *P, fq_ctx_t F);
void point_clear(point_t *P, fq_ctx_t F);

// Set P = Q
void point_set(point_t *P, point_t *Q, fq_ctx_t F);

// Print P to stdout in a nice format (no line break)
void point_pretty_print(point_t *P, fq_ctx_t F);

// Rescale P so that first nonzero coordinate is 1
void point_normalize(point_t *P, fq_ctx_t F);

// Return 1 if P = Q as points in AA^4(F), 0 otherwise.
// (This can be used to test equality of projective points
// if we already know they are normalized.)
int point_eq_raw(const point_t *P, const point_t *Q, fq_ctx_t F);

// Return 1 if P = Q, 0 otherwise.
// This normalizes P and Q in the process.
int point_eq(point_t *P, point_t *Q, fq_ctx_t F);

// Set Q = f(P), where f is the p^d-Frobenius map on F = GF(p^*)
void point_frobenius(point_t *Q, const point_t *P, long d, fq_ctx_t F);

// Set Q = f(P), where f is the q-Frobenius map on F = GF(q^2).
// This function only works when q is prime or the square of a prime.
void point_quadratic_frobenius(point_t *Q, const point_t *P, fq_ctx_t GFq_squared);

// View the array *vv as a vector of integers in [0,bound] of length len_vv.
// Modify the vector to be the next in the lex ordering if such exists and
// return 1. The first vector in this ordering is [0,...,0]. If
// *vv = [bound,...,bound], then there is no next vector; return 0.
// No validity testing is done on the array.
int next_vector(int8_t *vv, int8_t len_vv, int8_t bound);

// View the array *vv as a vector of integers in [0,bound] of length len_vv.
// Modify the vector to be the next in the lex ordering subject to the
// additional restriction that the leftmost nonzero entry must be 1, if such
// exists, and return 1. The first vector in this ordering is [0,...,0,1].
// If *vv = [1,bound,...,bound], then there is no next vector; return 0.
// No validity testing is done on the array.
// NOTE: We treat 0 and 1 like they are field elements, and in order to apply
//   this to actual fields, we will need to ensure that our lists of field
//   elements start with 0,1. 
int next_projective_vector(int8_t *vv, int8_t len_vv, int8_t bound);

// Simple code for testing 
void test_next_projective_vector(int8_t m, int8_t n);

/*----------------------------------------------------------------------------*/
/*---------------------------------  FIELDS  ---------------------------------*/
/*----------------------------------------------------------------------------*/

// Allocate memory for and a list of the field elements in GFq
void field_elements_array(fq_t **elts, fq_ctx_t GFq);

// Allocate memory for and make a list of the field elements in GFq
// that are fixed by x -> x^(p^d), where q = p^degree.
void subfield_elements_array(fq_t **elts, fq_ctx_t GFq, long d);

// Construct the field GF(q^2) and make an array of its elements, with
// the elements of GF(q) coming first. This code only works for q = 2, 3, 4.
void construct_fields(long q, fq_ctx_t GFq_squared, fq_t **quad_elts);

// construct addition table for field, with elements ordered in fld_elts.
// Space will be allocated for the table. The sum of fld_elts[i] and fld_elts[j]
// is stored at position i + q*j in the table. 
void construct_add_table(int8_t **add_table, fq_t *fld_elts, fq_ctx_t GFq);

// add field indices for GF(q) using the addition table
int8_t fq_add_table(int8_t *op1, int8_t *op2, int8_t *add_table, long q);

// construct multiplciation table for field, with elements ordered in fld_elts.
// Space will be allocated for the table. The product of fld_elts[i] and fld_elts[j]
// is stored at position i + q*j in the table. 
void construct_mul_table(int8_t **mult_table, fq_t *fld_elts, fq_ctx_t GFq);

// multiply field indices for GF(q) using the multiplication table
int8_t fq_mul_table(int8_t *op1, int8_t *op2, int8_t *mul_table, long q);

// dot product of vectors of field indices using addition/multiplication tables
int8_t fq_vec_dot_table(int8_t *op1, int8_t *op2, long len,
		      int8_t *add_table, int8_t *mul_table, long q);

/*----------------------------------------------------------------------------*/
/*---------------------------------  CURVES  ---------------------------------*/
/*----------------------------------------------------------------------------*/

// Monomial ordering:
// [x^3, x^2*y, x^2*z, x^2*w, x*y^2, x*y*z, x*y*w, x*z^2, x*z*w, x*w^2,
//  y^3, y^2*z, y^2*w, y*z^2, y*z*w, y*w^2, z^3, z^2*w, z*w^2, w^3].
typedef fq_struct* cubic_coef_vec_t;

void cubic_coef_vec_init(cubic_coef_vec_t *vec, fq_ctx_t GFq);
void cubic_coef_vec_clear(cubic_coef_vec_t *vec, fq_ctx_t GFq);

// Print vv to standard out
void cubic_coef_vec_pretty_print(cubic_coef_vec_t vv, fq_ctx_t GFq);

// Evaluate the initialized cubic monomial vector vv at point P.
void cubic_coef_vec_eval(cubic_coef_vec_t vv, point_t *P, fq_ctx_t GFq);

// Allocate memory and store the vector of cubic monomial evaluations
// v_P = (m(P) : m \in cubic monomials), for P in pts
void cubic_monomial_evaluations(cubic_coef_vec_t **vecs, point_t *pts, long num_pts,
				fq_t *field_elts, fq_ctx_t GFq);


// Print this cubic to fp in a format that can be read by Sage, with a newline
// character at the end.
void print_cubic_to_file_old(fq_struct *cubic, FILE *fp, fq_ctx_t F);
  
// Print this cubic to fp in a format that can be read by Sage, with a newline
// character at the end. Drop terms with zero coefficient.
void print_cubic_to_file(fq_struct *cubic, FILE *fp, fq_ctx_t F);


// Print this cubic to screen in a format that can be read by Sage,
// with a newline character at the end. 
void print_cubic_to_screen(fq_struct *cubic, fq_ctx_t F);

// Find all quadratic points up to Galois conjugacy on the standard nonsplit
// quadric surface over GF(q) given by xy + z^2 + azw + w^2 = 0, where
// q = 2: a = 1
// q = 3: a = 0
// q = 4: a = T^2 + T, where T^4 + T + 1 = 0.
// This function will allocate memory, initialize the points, and
// return the number of points it finds (which we expect to be
// q^2 + q + binomial(q^2+1,2)).
long nonsplit_quadric_surface_points(long q, point_t **quadric_pts,
				     fq_t *field_elts, fq_ctx_t GFq_squared);

// translate a cubic_coef_vec_t object into a vector of indices,
// given by the ordering in fld_elts. This is not a particularly fast
// implementation, so it shouldn't be used in the inner loop. 
void cubic_coef_vec_to_index_vec(int8_t *ind_vec, cubic_coef_vec_t cc_vec,
				 fq_t *fld_elts, fq_ctx_t GFq);
					      

/*----------------------------------------------------------------------------*/

#endif /* GENUS4_H */
