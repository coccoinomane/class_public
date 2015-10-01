/** @file song_tools.h Documented header file for song_tools.c */

#ifndef __SONG_TOOLS__
#define __SONG_TOOLS__

#include "common.h"

#ifndef _SPLINE_NATURAL_
#define _SPLINE_NATURAL_ 0 /**< natural spline: ddy0=ddyn=0 */
#endif

#ifndef _SPLINE_EST_DERIV_
#define _SPLINE_EST_DERIV_ 1 /**< spline with estimation of first derivative on both edges */
#endif



// ======================================================================================
// =                                     Binary files                                   =
// ======================================================================================

/**
 * Structure representing a binary file made of contiguous data blocks.
 *
 * Each block in the data file is homogeneous in that it contains the
 * same data type.
 */

struct binary_file {
  
  long int n_blocks; /**< Number of blocks in the binary file */

  struct binary_block ** block_array; /**< Array of binary blocks in the file */
  
  long int size_bytes; /**< Size of the binary file in bytes */

  short has_header; /**< Should the binary file have a human readable ASCII header? */

  int header_size; /**< Size of the header string */

  char * header; /**< Header of the binary file. This is the human readable (ASCII) part of 
                 the binary file. It can be read by opening the binary file with a text
                 editor. */

  int header_size_max; /**< Maximum size of the header string */
  
  short convert_to_float; /**< Should we convert all blocks with double precision data (double) to single
                          precision (float), in order to reduce the size of the binary file? */

  char path[_FILENAMESIZE_]; /**< String with the path of the binary file */

  FILE * stream; /**< File stream pointing to the binary file (input to fwrite) */
  
  ErrorMsg error_message; /**< Area to write error messages */
  
};


/**
 * Structure representing a data block in a binary file.
 *
 * A binary block is a homogeneous sequence of data of the same type inside a
 * binary file.
 *
 * For example, an array of N integers represents a data block with:
 * - size = N
 * - type_size = sizeof(int)
 * - size_bytes = N * sizeof(int)
 * - type = "int"
 * - description = "an array of N integers"
 *
 * A block can contain any number of data entries, each representing a
 * contiguous chunk of memory.
 * 
 */

struct binary_block {

  long int id; /**< Unique sequential identifier for this block in the parent file */

  int n_data; /**< Number of data entries in this block */

  struct binary_data ** data_array; /**< List of structure each representing a data entry in the block */
  
  int size; /**< Length of the block, in units of the data type of the block. It is given by the sum
            of the sizes of the single data entries. */  

  long int size_bytes; /**< Size of the block in bytes */

  int type_size; /**< Size of the data type in this block */

  long int start_byte; /**< Location of the block in the parent binary file, in bytes */

  char name[256]; /**< String with the name of the variable used to store the data in the program calling this tool */

  char description[512]; /**< String with a short description of the data contained in the block */

  char type[62]; /**< String with the data type of the block */
  
  short convert_to_float; /**< Should we convert this block data from double precision to single precision? */
  
};


/**
 * Structure representing a data entry inside a block of a binary file.
 *
 * A data entry is a contiguous chunk of memory inside a block. It is
 * completely characterised by a pointer to a memory location and a size.
 */

struct binary_data {
  
  void * pointer; /**< Pointer to the data location in memory */

  int size; /**< Number of elements to write to file for this data */  
  
};

/**
 * Maximum length of a line in the header file.
 */
#define _MAX_LINE_LENGTH_ 1024




/**
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

// ======================================================================================
// =                                   3J and 6J symbols                                =
// ======================================================================================

  int threej_single(
         int l1, int l2, int l3, int m2, int m3, // In
         double *threej,                         // Out
         ErrorMsg errmsg       
         );

  int threej_ratio_L_recursive (
        int l1, int l2, int l3, int N,     // In
        double *result,                    // Out, should be allocated with M+1 elements
        ErrorMsg errmsg
        );
        
  int threej_ratio_L (
        int l1, int l2, int l3,            // In
        int N1, int N2, int N3,            // In
        double *result,                    // Out
        ErrorMsg errmsg
        );

  int threej_ratio_L1 (
        int l1, int l2, int l3, int N,     // In
        double *result,                    // Out
        ErrorMsg errmsg
        );

  int threej_ratio_M_recursive (
        int l1, int l2, int l3, int M,     // In
        double *result,                    // Out, should be allocated with M+1 elements
        ErrorMsg errmsg
        );
        
  int threej_ratio_M (
        int l1, int l2, int l3, int M,     // In
        double *result,                    // Out
        ErrorMsg errmsg
        );
        
  int threej_A (
        int l1, int l2, int l3,
        int m1,
        double *result,
        ErrorMsg errmsg
        );

  int threej_B (
        int l1, int l2, int l3,
        int m1, int m2, int m3,
        double *result,
        ErrorMsg errmsg
        );

  int threej_C (
        int l1, int l2, int l3,
        int m2, int m3,
        double *result,
        ErrorMsg errmsg
        );

  int threej_D (
        int l1, int l2, int l3,
        int m2, int m3,
        double *result,
        ErrorMsg errmsg
        );

  int threej_A_factor (
        int l1, int l2, int l3,
        int n,
        double *result,
        ErrorMsg errmsg
        );

    

  // ======================================================================
  // =                           Bessel functions                         =
  // ======================================================================

  double spherical_bessel_j(
         int l,
         double x
         );
         

  // ====================================================================================
  // =                                Coupling factors                                  =
  // ====================================================================================

  double coupling_c_plus (int l, int m1, int m);
  double coupling_c_minus (int l, int m1, int m);
  double coupling_d_zero (int l, int m1, int m);
  double coupling_d_minus (int l, int m1, int m);
  double coupling_d_plus (int l, int m1, int m);    

  int coupling_general (
    int l2, int l3, int m1, int F,
    double * three_j_000, /* should be preallocated with at least l2_max doubles */
    int three_j_000_size,
    double * three_j_mmm, /* should be preallocated with at least m1_max doubles */
    int three_j_mmm_size,
    int * l1_min, int * l1_max,
    int * m2_min, int * m2_max,
    double ** result,     /* should be preallocated with at least l2_max*m1_max doubles */
    ErrorMsg errmsg 
    );
    

  // ============================================================================================
  // =                                 Legendre polynomials                                     =
  // ============================================================================================

  double plegendre_lm(int l, int m, double x);
  double plegendre_lm_rescaled(int l, int m, double x);
  double plegendre_lm_rescaled_analytically (int l, int m, double x);
  double plegendre (int n, double x);


  // ============================================================================================
  // =                             Multipole related functions                                  =
  // ============================================================================================
  
  int multipole2offset_l_m(int l, int m, int m_max);
  int size_l_m(int l_max, int m_max);
	int multipole2offset_l_indexm(int L, int M, int * m_vec, int m_size);
  int offset2multipole_l_indexm (int offset, int l_max, int * m_vec, int m_size,
                                 int * L, int * index_M);
	int size_l_indexm(int l_max, int * m_vec, int m_size);
  int multipole2offset_indexl_indexm(int L, int M, int * l_vec, int l_size, int * m_vec, int m_size);
	int size_indexl_indexm(int * l_vec, int l_size, int * m_vec, int m_size);
  int offset2multipole_indexl_indexm(int offset, int * l_vec, int l_size, int * m_vec, int m_size,
                       int * index_L, int * index_M);
  int multipole2offset_unconstrained_n_l_m(int n, int l, int m, int l_max, int m_max);
  int multipole2offset_n_l_m(int n, int l, int m, int l_max, int m_max);
  int size_n_l_m(int n_max, int l_max, int m_max);
  int multipole2offset_n_l_indexm(int N, int L, int M, int l_max, int * m_vec, int m_size);
  int size_n_l_indexm (int n_max, int l_max, int * m_vec, int m_size);
  

  // =====================================================================================
  // =                                     Interpolation                                 =
  // =====================================================================================

  int spline_sources_derivs(
			     double * x, /* vector of size tau_size */
			     int tau_size,
			     double *** y_array, /* array of size tau_size*tp_size with elements 
						  										y_array[index_tau*tp_size+index_tp] */
			     int tp_size,   
			     double *** ddy_array, /* array of size tau_size*tp_size */
			     short spline_mode,
           int index_mode,
           int index_ic,           
           int index_k,
           int k_size,
			     ErrorMsg errmsg
           );
           
  int spline_derivs_two_levels(
  			     double * x, /* vector of size tau_size */
  			     int tau_size,
  			     double ** y_array,
  			     int tp_size,   
  			     double ** ddy_array,
  			     short spline_mode,
  			     ErrorMsg errmsg
             );

  int spline_interpolate_two_levels(
  			     double * x_array,
  			     int tau_size,
  			     double ** y_array,
  			     double ** ddy_array,
  			     int tp_size,
  			     double x,
  			     int * last_index,
  			     double * result,
  			     int result_size, /** from 1 to tp_size */
  			     ErrorMsg errmsg
             );

  int spline_sources_interpolate_growing_closeby(
  			     double * x_array,
  			     int tau_size,
  			     double *** y_array,
  			     double *** ddy_array,
  			     int tp_size,
  			     double x,
  			     int * last_index,
  			     double * result,
  			     int result_size, /** from 1 to tp_size */
             int index_mode,
             int index_ic,           
             int index_k,
             int k_size,
  			     ErrorMsg errmsg
  			     );
           
  int spline_interpolate_two_levels_growing_closeby(
  			     double * x_array,
  			     int tau_size,
  			     double ** y_array,
  			     double ** ddy_array,
  			     int tp_size,
  			     double x,
  			     int * last_index,
  			     double * result,
  			     int result_size, /** from 1 to tp_size */
  			     ErrorMsg errmsg
             );

  int spline_sources_interpolate(
  			     double * x_array,
  			     int x_size,
  			     double *** y_array,
  			     double *** ddy_array,
  			     int tp_size,
  			     double x,
  			     int * last_index,
  			     double * result,
  			     int result_size, /** from 1 to tp_size */
             int index_mode,
             int index_ic,           
             int index_k,
             int k_size,
  			     ErrorMsg errmsg
             );

  int array_interpolate_linear_nospline(
			       double * x_array,
			       int n_lines,
			       double * array,
			       double * array_splined,
			       int n_columns,
			       double x,
			       int * last_index,
			       double * result,
			       int result_size, /** from 1 to n_columns */
			       ErrorMsg errmsg);
	
  int array_spline_derive_table_lines(
             double * x_array,
             int x_size,
             double * y_array,
             double * ddy_array,
             int y_size,
             double * dy_array,
             ErrorMsg errmsg);


  // ====================================================================================
  // =                                Sampling related                                  =
  // ====================================================================================

  int log_space (double * xx, double x_min, double x_max, int n_points);
  int lin_space (double * xx, double x_min, double x_max, int n_points);



  // ====================================================================================
  // =                                Matrix operations                                 =
  // ====================================================================================

  double Determinant(double **a,int n);
  void CoFactor(double **a,int n,double **b);
  void Transpose(double **a,int n);
  void InverseMatrix(double **in,int n,double **out);
  void PrintMatrix(double **in,int n);
  


  // ====================================================================================
  // =                                Assert functions                                  =
  // ====================================================================================

  int is_triangular_int (int l1, int l2, int l3);
  int is_triangular_double (double l1, double l2, double l3);



  // ====================================================================================
  // =                                Fortran functions                                 =
  // ====================================================================================
  
  /* The following functions are Fortran procedures from the Slatec library
  in the file tools/slatec_3j_f90.f90 */
  
  /**
   * Function from the Slatec library to compute the 3j symbol for
   * all the allowed values of 'm2'. 
   */
  void drc3jm_ (double *l1, double *l2, double *l3,
                double *m1, double *m2_min, double *m2_max,
                double *thrcof, int *ndim, int *ier);
	       

  /** 
   * Function from the Slatec library to compute the 3j symbol for
   * all the allowed values of 'l1'.  IMPORTANT: this is a Fortran
   * function contained in the file tools/slatec_3j_f90.f90
   */
  void drc3jj_ (double *l2, double *l3,
                double *m2, double *m3, double *l1_min, double *l1_max,
                double *thrcof, int *ndim, int *ier);

  /**
   * Function from the SLATEC library to compute the J Bessel function J_l ( x )
   * for values of l going from l to l+N-1.  IMPORTANT: this is a Fortran
   * function contained in the file tools/slatec_3j_f90.f90
   */
  void dbesj_ (double *x, double *l, int *N, double *result, int *NZ); // double precision
  void besj_ (float *x, float *l, int *N, float *result, int *NZ);     // single precision



  // ====================================================================================
  // =                                  Binary files                                    =
  // ====================================================================================

  int binary_init (
    struct binary_file * file,
    FILE ** file_stream,
    char * file_path,
    char * mode,
    int header_size,
    short convert_to_float
    );

  int binary_free (
    struct binary_file * file
    );

  int binary_append (
    struct binary_file * file,
    void * data_pointer,
    int size,
    int type_size,
    char * description,
    char * type,
    char * name
    );

  int binary_append_int (
    struct binary_file * file,
    void * data_pointer,
    int size,
    char * description,
    char * name
    );

  int binary_append_double (
    struct binary_file * file,
    void * data_pointer,
    int size,
    char * description,
    char * name
    );

  int binary_append_string (
    struct binary_file * file,
    void * data_pointer,
    int size,
    char * description,
    char * name
    );

  int binary_add_block (
    struct binary_file * file,
    void * data_pointer,
    int size,
    int type_size,
    char * description,
    char * type,
    char * name,
    long int block_id
    );

  int binary_add_data (
    struct binary_file * file,
    long int block_id,
    void * data_pointer,
    int data_size
    );

  int binary_append_data (
    struct binary_file * file,
    void * data_pointer,
    int data_size
    );

  int binary_add_data_to_block (
    struct binary_block * block,
    void * data_pointer,
    int data_size,
    ErrorMsg errmsg
    );

  int binary_delete_block (
    struct binary_file * file,
    long int block_id
    );

  int binary_change_size (
    struct binary_file * file,
    long int block_id,
    long int index_data,
    int new_size
    );

  int binary_change_type (
    struct binary_file * file,
    long int block_id,
    char * new_type,
    int new_type_size
    );

  int binary_write (
    struct binary_file * file
    );

  int binary_add_header_line (
    struct binary_file * file,
    char * line
    );

  int binary_sprintf (
    struct binary_file * file,
    const char * format,
    ...
  );

  int binary_add_header_map (
    struct binary_file * file,
    int * map_size
    );



  // ====================================================================================
  // =                                      Misc                                        =
  // ====================================================================================

  double identity_double (double x);

  int sign_int (int x);

  int ordering_int (
        int * n,         
        int * ordering,  
        ErrorMsg errmsg
        );

  int reorder_int (
        int * n,        
        int * ordering, 
        ErrorMsg errmsg
        );

  int merge_arrays_double (
        double *v1,
        int v1_size,
        double *v2,
        int v2_size,
        double **out,
        int * out_size,
        int (*compar)(const void *, const void *),
        ErrorMsg errmsg
        );
  

#ifdef __cplusplus
}
#endif


// ===============================================================================
// =                              Preprocessor macros                            =
// ===============================================================================

/* Needed by spherical_bessel_j */
#define _GAMMA1_ 2.6789385347077476336556
#define _GAMMA2_ 1.3541179394264004169452


#endif
