#ifndef __MESH__
#define __MESH__

#include "common.h"

/* Uncomment to use this module outside CLASS */
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


/**
 * Type of functions meant to compute the distance between two points
 */
typedef double distance_function_type (double * vec1, double * vec2);


/**
 * Structure containing the interpolation table & weights to interpolate
 * a function defined on an arbitrary domain in 2D or 3D space (the mesh).
 *
 * The interpolation is done via the functions mesh_2D_int() and mesh_3D_int();
 * the structure needs to be initialised first by calling either mesh_2D_sort()
 * or mesh_3D_sort().
 */

struct interpolation_mesh {

  int n_dim; /**< Dimensionality of the domain of the function to interpolate. So far,
             only 2D (n_dim=2) and 3D (n_dim=3) interpolation is supported. */

  int n_nodes; /**< Number of support points of the function to be interpolated */

  double link_length; /**< Extent of the local region influencing the interpolation; for
                      a grid, this should correspond roughly to the largest distance of
                      two neighbouring points. Nodes farther away than the interpolated
                      point will not influence the interpolated value. */

  double group_length; /**< Grouping length, used to downweight clusters of many close
                       points, so that they count as one and do not bias the interpolation.
                       For a grid, the grouping length should correspond to the shortest
                       distance between two neighbouring points. */
  
  double soft_coeff;  /**< Parameter to softens the linking length; the larger it is, the
                      larger the linking length will be. The default value on most meshes
                      should be 0.5 */

  double max;         /**< Maximum value for the arguments of the function to be
                      interpolated. Nodes with any coordinate larger than this
                      value will be ignored. */

  int n_bins;   /**< Number of linear regions in which the domain will be split.
                It is given by l_max/bin_size. */

  double bin_size;   /**< Length of a bin side. It is given by link_length*(1+soft_coeff) */

  int n_bins_z;   /**< Number of linear regions in which the third dimension will be split.
                  It is different from n_bins only for 2D interpolation, in which case it is
                  equal to one. That is, we treat the 2D interpolation like a 3D interpolation
                  with a single bin in the 3rd dimension. */

  distance_function_type * distance;  /**< Function used to compute the distance between two
                                      points */

  int * ix;  /**< List of x-bin where the node points belong */

  int * iy;  /**< List of y-bin where the node points belong */

  int * iz;  /**< List of z-bin where the node points belong */

  int *** grid;   /**< Grid in 3D space of size n_bins^3. It associates to the (ix,iy,iz)
                  bin the number of nodes contained in the bin. For 2D interpolation, the
                  iz level is ignored. The grid is completely determined by the domain
                  (x,y,z) and is thus independent from the function values f(x,y,z).*/

  int **** id;    /**< Assigns an identifier to each binned node. The identifier is just the index
                  of the node in the values array. The id is needed to establish a correspondence
                  between the binned nodes and their coordinates & values. Like the grid, the ID
                  is completely determined by the domain (x,y,z) and is thus independent from the
                  function values f(x,y,z). */

  double (*values)[5]; /**< Array containing for each node its coordinates, the value at the node
                       of the function to be interpolated, and the local density of nodes. The
                       node density is computed with respect to the nodes in the same bin; the
                       bin size is determined using the linking length. This array is indexed
                       as follows:
                          
                        values[n][0] -> value of f in the node n
                        values[n][1] -> x coordinate of the node n
                        values[n][2] -> y coordinate of the node n
                        values[n][3] -> z coordinate of the node n
                        values[n][4] -> density of nodes around the node n.
    
                      For 2D interpolation, the z coordinate is ignored. */

  int ready; /**< If this flag is true, the mesh (either 2D or 3D) has been filled
             and is ready to be interpolated */
  
  short compute_grid; /**< If this flag is true, the grid will be recomputed for each
                      call of mesh_2D_sort() and mesh_3D_sort(). If false, the grid
                      has been provided by the user as an argument. */

  /** Counters to keep track of the memory usage */
  //@{
  int n_allocated_in_mesh;
  int n_allocated_in_grid;
  //@}

  /** Counters to keep track of which meshes are used */
  //@{
  int count_interpolations;
  int count_range_extensions;
  //@}

  ErrorMsg error_message; /**< Area to write error messages */

};




/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int mesh_init (
        int n_dim,
        int n_nodes,
        double (*values)[4],
        double max,
        double link_length,
        double group_length,
        double soft_coeff,
        int *** grid,
        int **** id,
        struct interpolation_mesh * mesh
        );

  int mesh_interpolate (
        struct interpolation_mesh * mesh,
        double x,
        double y,
        double z,
        double * interpolated_value
        );

  int mesh_free (
        struct interpolation_mesh * mesh
        );

  double distance_3D (
           double * vec1,
           double * vec2
           );

  double distance_2D (
           double * vec1,
           double * vec2
           );
 
#ifdef __cplusplus
}
#endif


#endif