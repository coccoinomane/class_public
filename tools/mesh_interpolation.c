#include "mesh_interpolation.h"

/**
 * Prepare a function f(x,y) or f(x,y,z) for interpolation with the mesh algorithm.
 *
 * The interpolation mesh thus initialised can be used to interpolate the function
 * with the mesh_interpolate() function. The function can be either defined on
 * a bidimensional (n_dim=2) or tridimensional (n_dim=3) domain.
 *
 * In detail, this function does:
 *
 * -# Determine the size of each bin in 3D space based on the input linking
 *    length. Nodes farther than the linking length won't influence the
 *    interpolation, so link_length is effectively a correlation scales.
 *
 * -# Create a grid, that is count the number of nodes in each bin of
 *    side equal to the linking length. This is stored in mesh->grid[ix][iy][iz].
 *
 * -# Assign to each node in the grid a unique identifier, in order to establish
 *    a correspondence between the binned nodes and their values (x,y,z,f).
 * 
 * -# Compute for each binned node its local density, ie. the density of nodes
 *    around it in the same bin.
 *
 * For each other node in the bin, the density gets a contribution of
 *  e^(-distance squared/grouping_length) so that a group of nodes clustered
 * on a scale smaller than the grouping length counts as one particle. The density
 * is needed only to downweight those points that are clustered so that they count
 * as one (see mesh_interpolate()).
 *
 * Note that so far the mesh only supports x>0, y>0, z>0.
 */

int mesh_init (
      int n_dim, /**< Input: Dimensionality of the function domain; set n_dim=2 for bidimensional
                (x,y) or n_dim=3 for tridimensional (x,y,z) */
      int n_nodes, /**< Input: Number of points where the function is known */
      double (*values)[4], /**< Input: Array containing the coordinates and the function value
                           for each of the nodes. Build it in the following way:
                             values[n][0] -> value of f in the node n
                             values[n][1] -> x coordinate of the node n
                             values[n][2] -> y coordinate of the node n
                             values[n][3] -> z coordinate of the node n
                           The ordering of the nodes (ie. which n goes first) is arbitrary. For
                           2D interpolation, the z coordinate is ignored. */
      double max, /**< Input: Nodes with a coordinate larger than this will be ignored */
      double link_length, /**< Input: Linking length (see header file) */
      double group_length, /**< Input: Grouping length (see haeder file) */
      double soft_coeff, /**< Input: Softening coefficients (see header file) */
      int *** grid, /**< Input: Use this grid instead of computing it from scratch; set to NULL to ignore */
      int **** id, /**< Input: Use this id array instead of computing it from scratch; set to NULL to ignore */
      struct interpolation_mesh * mesh /**< Output: The mesh to initialise */
      )
{

  class_test (mesh->ready == _TRUE_,
    mesh->error_message,
    "you are trying to initialise a mesh that is already ready to use");

  /* Update the mesh structure with the input parameters */
  mesh->n_dim = n_dim;
  mesh->n_nodes = n_nodes;
  mesh->max = max;
  mesh->link_length = link_length;
  mesh->group_length = group_length;
  mesh->soft_coeff = soft_coeff;

  /* Uncomment to adjust the linking length automatically. The criterion
  is to have at least n points in each interpolation box. An interpolation
  box is a volume centred on a bin and including the adjacent bins. This
  adjustment is likely to increase the linking length for n ~ 1. Large
  linking lengths will create large bins and slow down the interpolation
  considerably, because more nodes are used, but it will also reduced memory
  consumption because less bins are needed. */
  // double n_points_per_box = 1;
  // double link_length_min = n_points_per_box * max / pow (n_nodes, 1.0/n_dim) / (1+soft_coeff) / n_dim;
  // if (link_length < link_length_min)
  //   mesh->link_length = link_length_min;

  /* Initialize counters and flags */
  mesh->n_allocated_in_grid = 0;
  mesh->n_allocated_in_mesh = 0;
  mesh->ready = _FALSE_;

  /* Do not compute the grid if the user provided one */
  mesh->compute_grid = (grid==NULL);

  /* Bin the domain of the function in bins according to the linking length
  and on the global maximum */
  mesh->bin_size = mesh->link_length * (1 + mesh->soft_coeff);
  mesh->n_bins = ceil(mesh->max/mesh->bin_size);

  /* If the function has no nodes or somehow we could not bin the domain, we
  flag the mesh as incomplete and move on */
  if ((mesh->n_nodes <= 0) || (mesh->n_bins <= 0)) {
    mesh->ready = _FALSE_;
    return _SUCCESS_;
  }



  // ====================================================================================
  // =                                 Dimensionality                                   =
  // ====================================================================================

  class_test ((mesh->n_dim != 2) && (mesh->n_dim != 3),
    mesh->error_message,
    "only 2D (n_dim=2) and 3D (n_dim=3) interpolation supported");

  /* In three dimensions, the z direction is not special and has the same number of
  bins as the x and y directions */
  if (mesh->n_dim == 3) {
    mesh->n_bins_z = mesh->n_bins;
    mesh->distance = distance_3D;
  }

  /* In two dimensions, we treat the z directions as containing a single bin with
  all the nodes */
  else if (mesh->n_dim == 2) {
    mesh->n_bins_z = 1;
    mesh->distance = distance_2D;
  }


  // ====================================================================================
  // =                                   Copy values                                    =
  // ====================================================================================

  mesh->values = calloc (5 * mesh->n_nodes, sizeof(*mesh->values));

  #pragma omp parallel for
  for (int n = 0; n < mesh->n_nodes; ++n) {

    mesh->values[n][0] = values[n][0]; /* function value at the node */
    mesh->values[n][1] = values[n][1]; /* x coordinate of the node */
    mesh->values[n][2] = values[n][2]; /* y coordinate of the node */
    if (n_dim > 2)
      mesh->values[n][3] = values[n][3]; /* z coordinate of the node */

  }
  


  if (mesh->compute_grid==_TRUE_) {

    // ====================================================================================
    // =                                   Create grid                                    =
    // ====================================================================================

    /* Bin the support points in a grid based on the linking length. The grid contains
    for each (ix,iy,iz) bin the number of nodes that belong to that bin. */

    /* Allocate grid */

    mesh->grid = (int***) malloc(mesh->n_bins*sizeof(int**));
    for (int i=0; i<mesh->n_bins; i++) {
      mesh->grid[i] = (int**) malloc(mesh->n_bins*sizeof(int*));
      for (int j=0; j<mesh->n_bins; j++) {
        mesh->grid[i][j] = (int*) calloc(mesh->n_bins_z, sizeof(int));
      }
    }

    mesh->n_allocated_in_grid =
      mesh->n_bins +
      pow (mesh->n_bins, 2) +
      pow (mesh->n_bins, 2) * mesh->n_bins_z * sizeof(int)/(double)sizeof(int*);


    /* Fill grid */

    int abort = _FALSE_;
    #pragma omp parallel for
    for (int n=0; n<mesh->n_nodes; n++) {

      double x = values[n][1];
      double y = values[n][2];
      double z = (mesh->n_dim<3)?0:values[n][3];
      
      class_test_parallel (x<0 || y<0 || z<0,
        mesh->error_message,
        "mesh interpolation so far only supports positive (x,y,z)");

      int ix = floor(x/mesh->bin_size);
      int iy = floor(y/mesh->bin_size);
      int iz = floor(z/mesh->bin_size);

      /* Skip the data that is larger than max */
      if ((ix >= mesh->n_bins) || (iy >= mesh->n_bins) || (iz >= mesh->n_bins_z))
        continue;

      #pragma omp atomic
      mesh->grid[ix][iy][iz]++;

    }

    if (abort == _TRUE_)
      return _FAILURE_;



    // ====================================================================================
    // =                                  ID the nodes                                    =
    // ====================================================================================

    /* Assign an identifier to each binned node. The identifier is just the
    index of the node in the mesh->values array. The id is needed to establish a
    correspondence between the binned nodes and their coordinates & values.
    Like the grid, the ID is completely determined by the domain (x,y,z) and is
    thus independent from the function values f(x,y,z). */


    /* Allocate id */
  
    int count = 0;
  
    mesh->id = (int****) malloc(mesh->n_bins*sizeof(int***));
    for (int i=0; i<mesh->n_bins; i++) {
      mesh->id[i] = (int***) malloc(mesh->n_bins*sizeof(int**));
      for (int j=0; j<mesh->n_bins; j++) {
        mesh->id[i][j] = (int **) malloc(mesh->n_bins_z*sizeof(int*));
        for (int k=0; k<mesh->n_bins_z; k++) {
          mesh->id[i][j][k] = (int *) malloc((mesh->grid[i][j][k])*sizeof(int));
        }
      }
    }

    mesh->n_allocated_in_mesh = 
      mesh->n_bins +
      pow (mesh->n_bins, 2) +
      pow (mesh->n_bins, 2) * mesh->n_bins_z +
      mesh->n_nodes * sizeof(int)/(double)sizeof(int*);

    class_test (mesh->n_allocated_in_mesh == 0,
      mesh->error_message,
      "empty mesh, something went wrong. No nodes?");


    /* Assign id. Note that this loop cannot be parallelised naively with:
      #pragma omp parallel for private (i,ix,iy,iz)
    because it relies on incrementing a quantity that, in the same loop, is 
    used to index an array (counter). */

    int (*counter)[mesh->n_bins][mesh->n_bins_z] =
      calloc(mesh->n_bins*mesh->n_bins*mesh->n_bins_z, sizeof(int));

    for (int n=0; n<mesh->n_nodes; n++) {

      int ix = floor(values[n][1]/mesh->bin_size);
      int iy = floor(values[n][2]/mesh->bin_size);
      int iz = (mesh->n_dim<3)?0:floor(values[n][3]/mesh->bin_size);

      /* Skip the data that is larger than mesh->max */
      if ((ix >= mesh->n_bins) || (iy >= mesh->n_bins) || (iz >= mesh->n_bins_z))
        continue;

      /* Place the current node in the (ix,iy,iz) bin */
      mesh->id[ix][iy][iz][counter[ix][iy][iz]] = n;
      counter[ix][iy][iz]++;

    }

    /* Verify that grid = counter */
    for (int i = 0; i<mesh->n_bins; i++)
      for (int j = 0; j<mesh->n_bins; j++)
        for (int k = 0; k<mesh->n_bins_z; k++)
          class_test (counter[i][j][k] != mesh->grid[i][j][k],
            mesh->error_message,
            "%d!=%d\n", __func__, counter[i][j][k], mesh->grid[i][j][k]);

    free (counter);

  }

  /* Use the user-provided grid */

  else {
    
    mesh->id = id;
    mesh->grid = grid;
    
  }

  
  // ====================================================================================
  // =                                 Local density                                    =
  // ====================================================================================
  
  /* For each node, compute the local density of nodes, and store it in the values
  array. This is a loop over the nodes, which means that it can be quick if there are
  not too many. */

  #pragma omp parallel for
  for (int i = 0; i<mesh->n_bins; i++) {
    for (int j = 0; j<mesh->n_bins; j++) {
      for (int k = 0; k<mesh->n_bins_z; k++) {
        for (int n = 0; n<mesh->grid[i][j][k]; n++) {
          
          int id_n = mesh->id[i][j][k][n];
          
          for (int m = 0; m<mesh->grid[i][j][k]; m++) { 

            int id_m = mesh->id[i][j][k][m];

            double dist = (*mesh->distance) (mesh->values[id_n], mesh->values[id_m]);
            double density = exp(-dist*dist/pow(mesh->group_length,2));

            #pragma omp atomic
            mesh->values[id_n][4] += density;

          }
        }
      }
    }
  }
  
  /* The mesh is finally ready and can be used to interpolate the function
  using mesh_interpolate() */
  mesh->ready = _TRUE_;
  
  return _SUCCESS_;
  
}



/**
 * Interpolate the function f in (x,y,z) using its pre-computed interpolation mesh.
 *
 * The interpolation mesh has to be initialised first with the mesh_init() function.
 * The interpolation is first attempted using the nodes in a volume centred on the bin
 * where the requested point belongs and including all the adjacent bins. If no nodes
 * are found, the volume is enlarged to include another layer of adjacent bins.
 */

int mesh_interpolate (
    struct interpolation_mesh * mesh,
    double x,
    double y,
    double z,
    double * interpolated_value
    )
{

  class_test (!mesh->ready,
    mesh->error_message,
    "cannot interpolate, mesh is not ready! No nodes? max too small?");


  // ====================================================================================
  // =                                   Bracket point                                  =
  // ====================================================================================

  #ifdef DEBUG
  class_test (x>mesh->max || y>mesh->max || z>mesh->max || x<0 || y<0 || z<0,
    mesh->error_message,
    "(x,y,z)=(%g,%g,%g) out of bounds", x, y, z);
  #endif

  /* Locate the bin where (x,y,z) belongs */
  int ix = floor(x/mesh->bin_size);
  int iy = floor(y/mesh->bin_size);
  int iz = (mesh->n_dim<3)?0:floor(z/mesh->bin_size);

  /* Start by considering only the bins adjacent to the one with (x,y,z). Note that
  for 2D interpolation, the third direction z only has one bin (izmin=izmax=0). */
  int ixmin = MAX (0,ix-1);
  int iymin = MAX (0,iy-1);
  int izmin = MAX (0,iz-1);
  int ixmax = MIN (mesh->n_bins-1,ix+1);
  int iymax = MIN (mesh->n_bins-1,iy+1);
  int izmax = MIN (mesh->n_bins_z-1,iz+1);
  
  /* Debug: print the brackets of (x,y,z) */
  // printf("%d:%d:%d %d:%d:%d \n",ixmin,ix,ixmax,iymin,iy,iymax);
  
  double r[4];
  r[1] = x;
  r[2] = y;
  r[3] = z;


  // ====================================================================================
  // =                                    Interpolate                                   =
  // ====================================================================================

  double norm = 0;
  double result = 0;
  double density = 0;
  int n_usable_nodes = 0;

loop:
  result = 0;
  density = 0;
  n_usable_nodes = 0;

  /* Consider only the boxes close to  */
  for (int ix = ixmin; ix<=ixmax; ix++) {
    for (int iy = iymin; iy<=iymax; iy++) {
      for (int iz = izmin; iz<=izmax; iz++) {
        for (int n = 0; n < mesh->grid[ix][iy][iz]; n++) {

          /* Uncomment to compute the true density around the node. This will give better
          accuracy, but also slow the code down a lot */
          // int status;
          // for (double ix2 = ixmin; ix2<=ixmax; ix2++) {
          //   for (double iy2 = iymin; iy2<=iymax; iy2++) {
          //     for (double iz2 = izmin; iz2<=izmax; iz2++) {
          //       for (double n2 = 0; n2<grid[ix2][iy2][iz2]; n2++) {
          //         status = gsl_sf_exp_e(-pow(mesh->distance(values[ix][iy][iz][n],values[ix2][iy2][iz2][n2]),2)/pow(mesh->group_length,2),&jresult);
          //         if (status == GSL_ERANGE || jresult.val != jresult.val)
          //           density = 0.0;
          //         else
          //           density = jresult.val;
          //         values[ix][iy][iz][n][4] += density;
          // }}}}

          #pragma omp atomic
          n_usable_nodes++;

          /* Identifier of the current node */
          int id_n = mesh->id[ix][iy][iz][n];

          /* We weight the contribution from each node with a 1/distance law, so that the
          closer the node the stronger the influence on the point (x,y,z). We also include
          a 1/density factor so that the contribution from nodes in high-density regions
          is penalised. The objective is to make these clusters of points to count as one. */
          double value = mesh->values[id_n][0];
          double dist = mesh->distance (mesh->values[id_n], r);
          double density = mesh->values[id_n][4];

          /* Compute the error function with the fast approximation in Abramowitz & Stegun,
          Sec. 7.1. The maximum error is 5e-4. */
          double x = (dist-mesh->link_length)/(mesh->link_length*mesh->soft_coeff);
          double x2 = x*x;
          double erfx = 1 - 1/(1+erf_a1*x+erf_a2*x2+erf_a3*x*x2+erf_a4*x2*x2);
     
          /* Uncomment to use the more accurate erf function */
          // double erfx = erf((dist-mesh->link_length)/(mesh->link_length*mesh->soft_coeff))
        
          double weight = 1/(density*(dist+10e-7*mesh->link_length)) * (1.00000001 - erfx);
  
          /* Each node contributes to the interpolation proportionally to its weight */
          #pragma omp atomic
          result += weight * value;

          /* The weights introduce a scale that we shall eliminate with a normalisation factor */
          #pragma omp atomic
          norm += weight;

          /* Debug: print the contribution to the interpolation from this node */
          // printf("norm:%f, density:%f, val:%f, dist:%f, val@dist:%g\n",
          //   weight, mesh->values[id_n][4], result/norm, dist, value);

        }
      }
    }
  }

  /* Irregular meshes might have very isolated nodes. If this is the case, we
  extend the region of interpolation to the adjacent bins. */
  if (n_usable_nodes < 1) {
    
    class_test (ixmin==0 && ixmax==(mesh->n_bins-1) &&
                iymin==0 && iymax==(mesh->n_bins-1) &&
                izmin==0 && izmax==(mesh->n_bins_z-1),
      mesh->error_message,
      "reached edges of the box without finding a node; n_nodes=0?");
    
    // printf("Increasing Range!!! \n");
    #pragma omp atomic
    mesh->count_range_extensions++;
    ixmin = MAX (0,ixmin-1);
    iymin = MAX (0,iymin-1);
    izmin = MAX (0,izmin-1);
    ixmax = MIN (mesh->n_bins-1,ixmax+1);
    iymax = MIN (mesh->n_bins-1,iymax+1);
    izmax = MIN (mesh->n_bins_z-1,izmax+1);
    goto loop; 
  }

  /* Update counter */
  #pragma omp atomic
  mesh->count_interpolations++;

  /* Return normalised value */
  *interpolated_value = result/norm;
  
  return _SUCCESS_;
    
}


/**
 * Free the memory associated to the mesh and set mesh->ready=false.
 */

int mesh_free (
    struct interpolation_mesh * mesh
    )
{
  
  if (mesh->ready) {

    /* Free grid */
    if (mesh->compute_grid==_TRUE_) {

      for (int i = 0; i<mesh->n_bins; i++) {
        for (int j = 0; j<mesh->n_bins; j++) {

          free (mesh->grid[i][j]);

        } free (mesh->grid[i]);
      } free (mesh->grid);


      /* Free mesh */
      for (int i = 0; i<mesh->n_bins; i++) {
        for (int j = 0; j<mesh->n_bins; j++) {
          for (int k = 0; k<mesh->n_bins_z; k++) {

            free (mesh->id[i][j][k]);

          } free (mesh->id[i][j]);
        } free (mesh->id[i]);
      } free (mesh->id);

    }

    /* Free values */
    free (mesh->values);

    mesh->ready = _FALSE_;

  }

  return _SUCCESS_;
  
}



/**
 * Compute the Euclidean distance between two points in 2D space
 */

double distance_2D (double * vec1, double * vec2) {
  
  return  sqrt ( (vec1[1]-vec2[1])*(vec1[1]-vec2[1])
                +(vec1[2]-vec2[2])*(vec1[2]-vec2[2]) );
}



/**
 * Compute the Euclidean distance between two points in 3D space
 */

double distance_3D (double * vec1, double * vec2) {
  
  return  sqrt ( (vec1[1]-vec2[1])*(vec1[1]-vec2[1])
                +(vec1[2]-vec2[2])*(vec1[2]-vec2[2])
                +(vec1[3]-vec2[3])*(vec1[3]-vec2[3]) );
}


