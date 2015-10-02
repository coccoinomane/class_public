/** @file binary.h Documented header file for binary.c */

#ifndef __BINARY__
#define __BINARY__

#include "common.h"
#include <errno.h>

/**
 * Structure representing a binary file made of contiguous data blocks.
 *
 * Each block in the data file is homogeneous in that it contains the
 * same data type.
 */

struct binary_file {
  
  long int n_blocks; /**< Number of blocks in the binary file */

  struct binary_block ** block_array; /**< Array of binary blocks in the file */
  
  off_t size_bytes; /**< Size of the binary file in bytes */

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

  off_t size_bytes; /**< Size of the block in bytes */

  int type_size; /**< Size of the data type in this block */

  long int start_byte; /**< Location of the block in the parent binary file, in bytes */

  char name[256]; /**< Name of the variable used to store the data in the calling program */

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
 * Maximum length of a line in the header part of a binary file.
 */
#define _MAX_HEADER_LINE_LENGTH_ 1024




/**
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

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


#ifdef __cplusplus
}
#endif


#endif
