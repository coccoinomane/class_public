#ifndef __MESH__
#define __MESH__

#include "common.h"

/* Header and definitions needed to compile as standalone */
// #include "math.h"
// #include "stdio.h"
// #define MIN(a,b) (((a)<(b)) ? (a) : (b) ) /**< the usual "min" function */
// #define MAX(a,b) (((a)<(b)) ? (b) : (a) ) /**< the usual "max" function */
// #define _SUCCESS_ 0
// #define _FAILURE_ 0


/** Weights for the error function approximation */
//@{
#define erf_a1 0.278393
#define erf_a2 0.230389
#define erf_a3 0.000972
#define erf_a4 0.078108
//@}


struct interpolation_mesh {

  long int n_points;
  long int n_boxes;
  double l_max;
  double link_length;
  double soft_coeff;
  double group_length;
  
  int *** grid_3D;
  double ***** mesh_3D;
  
  int ** grid_2D;
  double **** mesh_2D;
  
  /** Should we compute the grid? If not, the user has to overwrite by hand the grid field
  with a precomputed grid. */
  short compute_grid;
  
  /** Counters to keep track of the memory usage */
  //@{
  long int n_allocated_in_mesh;
  long int n_allocated_in_grid;
  //@}

  /** Counters to keep track of which meshes are used */
  //@{
  long int count_interpolations;
  long int count_range_extensions;
  //@}

  ErrorMsg error_message; /**< Area to write error messages */

};




/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  // ====================================================================================
  // =                                3D interpolation                                  =
  // ====================================================================================

  int mesh_3D_sort (
      struct interpolation_mesh * mesh,
      int *** grid,
      double (*values)[4]
      );

  int mesh_3D_int (
      struct interpolation_mesh * mesh,
      double x,
      double y,
      double z,
      double * interpolated_value
      );

  int mesh_3D_free (
      struct interpolation_mesh * mesh
      );

  double distance_3D (
      double * vec1,
      double * vec2
  );


  // ====================================================================================
  // =                                2D interpolation                                  =
  // ====================================================================================

  int mesh_2D_sort (
      struct interpolation_mesh * mesh,
      int ** grid,
      double (*values)[3]
      );

  int mesh_2D_int (
      struct interpolation_mesh * mesh,
      double x,
      double y,
      double * interpolated_value
      );

  int mesh_2D_free (
      struct interpolation_mesh * mesh
      );

  double distance_2D (
      double * vec1,
      double * vec2
  );

#ifdef __cplusplus
}
#endif


#endif