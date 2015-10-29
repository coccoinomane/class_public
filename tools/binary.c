/** @file binary.c 
 *
 * Library of functions to write binary files with a meaningful information
 * header.
 *
 * Created by Guido Walter Pettinari on 02.10.2015
 * Last edited by Guido Walter Pettinari on 02.10.2015
 */

#include "binary.h"
#include "math.h"


/**
 * Prepare a binary file for being written.
 *
 * Initialise accumulators and open the actual file. If the file
 * already exists, it will be overwritten.
 */

int binary_init (
  struct binary_file * file, /**< Binary file being initialised */
  FILE ** file_stream, /**< File stream to link with the file */
  char * file_path, /**< String with the path to the file */
  char * mode, /**< Mode: w for writing, a for appending */
  int header_size_max, /**< Size of the header part of the binary file */
  short convert_to_float /**< Convert double precision entries to single precision? */
  )
{

  /* Initialise blocks and bytes number */
  file->n_blocks = 0;
  file->size_bytes = 0;
  file->header_size = 0,
  file->block_array = NULL; /* for realloc to work */
  
  /* Update single precision flag */
  file->convert_to_float = convert_to_float;
  
  /* Update the header size */
  file->header_size_max = header_size_max;
  file->has_header = (header_size_max > 0) ? _TRUE_: _FALSE_;

  /* Open file for writing */
  class_test ((mode[0] != 'w') && (mode[0] != 'a'),
    file->error_message,
    "mode can be either append (a) or write (w)");
  class_open (*file_stream, file_path, mode, file->error_message);
  
  /* Update structure with file stream and path */
  file->stream = *file_stream;
  strcpy (file->path, file_path);
  
  /* Include the ASCII header as the first block of the binary file */
  if (file->has_header == _TRUE_) {

    /* Allocate the header */
    class_calloc (file->header, header_size_max, sizeof(char), file->error_message);

    /* Add the header block */
    class_call (binary_append (
                  file,
                  file->header,
                  file->header_size_max,
                  sizeof (char),
                  "the header you are reading",
                  "char",
                  "file->header"),
      file->error_message,
      file->error_message);
  }
  
  return _SUCCESS_;
  
}



/**
 * Add a block to the binary file and initialise it with a data object.
 *
 * This function does not write the data block to file; it just takes
 * note of the memory location that needs to be written. The actual
 * writing is performed at a later time by binary_write().
 *
 * This means that the data pointer will not be accessed by this function;
 * the data it points to is accessed only when binary_write() is called.
 *
 * The block will be inserted in the binary file at the position given by
 * block_id, shifting all other blocks to the right. To append it
 * at the end of the binary file, set block_id=file->n_blocks.
 *
 * If the block exists in the binary file with the same name, description
 * and data type, then we assume the user wants to add data to that existing
 * block rather than creating a new block. This feature allows to call
 * binary_add_block() in a loop in order to add non-contiguous memory locations
 * to a single block.
 */

int binary_add_block (
  struct binary_file * file, /**< Binary file to which the block will be added */
  void * data_pointer, /**< Pointer to the data to be written in the new block */
  int data_size, /**< Number of elements to store in the block */
  int type_size, /**< Size in bytes of each element belonging to the block; usually given by sizeof(data type) */
  char * description, /**< Short description for the data being saved in the new block */
  char * type, /**< String with the data type of the new block (eg. "int" or "double") */
  char * name, /**< String with the name of the variable holding the data that will be written in the new block */
  long int block_id /**< Position of the new block in the binary file; all pre-existing blocks following the new block will be shifted forward in memory */
  ) 
{

  // ====================================================================================
  // =                                      Checks                                      =
  // ====================================================================================

  class_test (data_pointer == NULL,
    file->error_message,
    "NULL pointer for block %d (name=%s)", block_id, name);

  class_test (data_size < 0,
    file->error_message,
    "wrong size (%d) for block %d (name=%s)", data_size, block_id, name);

  class_test (type_size <= 0,
    file->error_message,
    "wrong type size (%d) for block %d (name=%s)", type_size, block_id, name);

  class_test (block_id > file->n_blocks,
    file->error_message,
    "cannot add block at position %d (min=0, max=%d)", block_id, file->n_blocks);

  class_test ((strcmp(type,"int")==0) && (type_size!=sizeof(int)),
     file->error_message,
     "'int' input has not the correct data size");

  class_test ((strcmp(type,"double")==0) && (type_size!=sizeof(double)),
     file->error_message,
     "'double' input has not the correct data size");

  class_test ((strcmp(type,"float")==0) && (type_size!=sizeof(float)),
     file->error_message,
     "'float' input has not the correct data size");

  class_test ((strcmp(type,"char")==0) && (type_size!=sizeof(char)),
     file->error_message,
     "'char' input has not the correct data size");

  class_test ((strcmp(type,"short")==0) && (type_size!=sizeof(short)),
     file->error_message,
     "'short' input has not the correct data size");

  

  // ====================================================================================
  // =                                   Special case                                   =
  // ====================================================================================

  /* If the block exists in the binary file with the same name, description
  and data type, then just add a new data entry to the existing block rather
  than creating a new block */

  short block_exists = 
    (file->n_blocks > block_id) &&
    (file->block_array[block_id]->type_size == type_size) &&
    (strcmp (file->block_array[block_id]->description, description) == 0) &&
    (strcmp (file->block_array[block_id]->type, type) == 0) &&
    (strcmp (file->block_array[block_id]->name, name) == 0);
                  
  if (block_exists) {

    class_call (binary_add_data (
                  file,
                  block_id,
                  data_pointer,
                  data_size),
      file->error_message,
      file->error_message);

    return _SUCCESS_;

  }


  // ====================================================================================
  // =                                   Create block                                   =
  // ====================================================================================

  /* Create and allocate the new block */
  struct binary_block * new_block;
  class_alloc (new_block, sizeof (struct binary_block), file->error_message);
  new_block->data_array = NULL; /* for realloc to work */

  /* The block starts with no data */
  new_block->n_data = 0;
  new_block->size = 0;
  new_block->size_bytes = 0;
  
  /* Give an identity to the block */
  new_block->id = block_id;
  new_block->type_size = type_size;
  strcpy (new_block->description, description);
  strcpy (new_block->type, type);
  strcpy (new_block->name, name);

  /* If asked, mark the block for conversion to single precision */
  new_block->convert_to_float = _FALSE_;
  if ((file->convert_to_float==_TRUE_) && (strcmp (type,"double")==0))
    new_block->convert_to_float = _TRUE_;

  /* If the block is to be appended at the end of the file, its first byte corresponds
  to the last byte of the file. If the block is to be added in place of another block,
  it inherits the first byte from the replaced block. */
  if (block_id == file->n_blocks)
    new_block->start_byte = file->size_bytes;
  else 
    new_block->start_byte = file->block_array[block_id]->start_byte;
  
  /* Add the data to the block. */
  class_call (binary_add_data_to_block (
                new_block,
                data_pointer,
                data_size,
                file->error_message),
    file->error_message,
    file->error_message);
    
  /* Debug - Print details on the block */
  // printf ("\n");
  // printf ("new_block->name = %s\n", new_block->name);
  // printf ("new_block->id = %d\n", new_block->id);
  // printf ("new_block->description = %s\n", new_block->description);
  // printf ("new_block->type = %s\n", new_block->type);
  // printf ("new_block->size = %d\n", new_block->size);
  // printf ("new_block->type_size = %d\n", new_block->type_size);
  // printf ("new_block->start_byte = %d\n", new_block->start_byte);

  

  // ====================================================================================
  // =                                   Insert block                                   =
  // ====================================================================================

  /* Increase by one the number of blocks in the binary file */
  file->n_blocks++;

  /* Create an empty slot at the end of the blocks array */
  class_realloc (
    file->block_array,
    file->block_array,
    file->n_blocks*sizeof(struct binary_block *),
    file->error_message
    );

  /* Modify the other blocks to account for the new block */
  for (int index_block=(file->n_blocks-1); index_block > block_id; --index_block) {

    /* Make room for the new block by shifting to the rigth the blocks with higher ID,
    starting from the last */
    file->block_array[index_block] = file->block_array[index_block-1];
    
    /* Update the block location in the binary file */
    file->block_array[index_block]->id += 1;
    file->block_array[index_block]->start_byte += new_block->size_bytes;
    
  }

  /* Insert the block at the right position */
  file->block_array[block_id] = new_block;

  /* Increase the size counter of the binary file */
  file->size_bytes += new_block->size_bytes;
  
  return _SUCCESS_;
  
}




/**
 * Write the binary file using the content of the file structure.
 * 
 * After this step, you cannot make further modifications to the
 * block structure of the binary file.
 */

int binary_write (
  struct binary_file * file
  )
{

  /* Downgrade the block to single precision, if requested */  
  for (int index_block=0; index_block < file->n_blocks; ++index_block)
    if (file->block_array[index_block]->convert_to_float == _TRUE_)
      class_call (binary_change_type (
                    file,
                    index_block,
                    "float",
                    sizeof(float)),
        file->error_message,
        file->error_message);
  
  
  // ====================================================================================
  // =                                   Write header                                   =
  // ====================================================================================
  
  if (file->has_header == _TRUE_) {

    /* Include the ASCII header as the first block of the binary file. We do it here
    after all other blocks have been assigned because the header contains information
    on all other blocks. */

    int map_size;

    class_call (binary_add_header_map (
                  file,
                  &map_size),
      file->error_message,
      file->error_message);
      
    /* Add a termination string to the header, without the comment sign */

    char termination[] = "END_OF_HEADER\n";

    class_test (
      (file->header_size+(strlen(termination)+1)) > file->header_size_max,
      file->error_message,
      "header reached maximum size (%d), will ignore further input",
      file->header_size_max);

    sprintf (file->header, "%s%s", file->header, "END OF HEADER\n");
    
    file->header_size += strlen(termination) + 1;
      
  }
  
  
  // ====================================================================================
  // =                                    Write data                                    =
  // ====================================================================================
  
  /* Sequentially write blocks to file */

  for (int index_block=0; index_block < file->n_blocks; ++index_block) {

    struct binary_block * block = file->block_array[index_block];

    /* Downgrade the block to single precision, if requested */
    if (block->convert_to_float == _TRUE_)
      class_call (binary_change_type (
                    file,
                    index_block,
                    "float",
                    sizeof(float)),
        file->error_message,
        file->error_message);

    /* Debug - Print info on the block being written */
    // printf ("Writing block #%7ld (start_byte=%7ld, size_bytes=%7ld, name=%s) to file %s\n",
    //   block->id, block->start_byte, block->size_bytes, block->name, file->path);

    for (int index_data=0; index_data < block->n_data; ++index_data) {

        struct binary_data * data = block->data_array[index_data];
        
        /* Downgrade the data from double to float, if requested */
        if (block->convert_to_float == _TRUE_) {
          float * float_pointer;
          class_alloc (float_pointer, data->size*sizeof(float), file->error_message);
          for (int i=0; i < data->size; ++i)
            float_pointer[i] = (float)(((double *)data->pointer)[i]);
          data->pointer = float_pointer;
        }
        
        /* Write to file the data entry */
        fwrite (data->pointer,
                block->type_size,
                data->size,
                file->stream);

        /* Free the temporary array */
        if (block->convert_to_float == _TRUE_)
          free (data->pointer);

    }
  }



  // ====================================================================================
  // =                                      Checks                                      =
  // ====================================================================================

  /* Check that the right number of bytes was written to file */
  fflush (file->stream);

  struct stat st; 

  class_test ( stat(file->path, &st) != 0,
    file->error_message,
    "Cannot determine size of %s: %s\n",
    file->path, strerror(errno));

  class_test (st.st_size != file->size_bytes,
    file->error_message,
    "Byte size of binary file (%d) doesn't match internal representation (%d)",
    st.st_size, file->size_bytes);
  

  return _SUCCESS_;
  
}


/**
 * Add a data block to the end of the binary file.
 *
 * This is a wrapper to binary_add_block(); see documentation there for
 * more details.
 */

int binary_append (
  struct binary_file * file, /**< Binary file to which the block will be appended */
  void * data_pointer, /**< Pointer to the data to be written in the new block */
  int size, /**< Number of elements to store in the block */
  int type_size, /**< Size in bytes of each element belonging to the block; usually given by sizeof(data type) */
  char * description, /**< Short description for the data being saved in the new block */
  char * type, /**< String with the data type of the new block (eg. "int" or "double") */
  char * name /**< String with the name of the variable holding the data that will be writtne in the new block */
  )
{
  
  class_call ( binary_add_block (
                 file,
                 data_pointer,
                 size,
                 type_size,
                 description,
                 type,
                 name,
                 file->n_blocks),
    file->error_message,
    file->error_message);
  
  return _SUCCESS_;
  
}

/**
 * Add a data block made of integers to the end of the binary file.
 *
 * This is a wrapper to binary_add_block(); see documentation there for
 * more details.
 */

int binary_append_int (
  struct binary_file * file, /**< Binary file to which the block will be appended */
  void * data_pointer, /**< Pointer to the data to be written in the new block */
  int size, /**< Number of elements to store in the block */
  char * description, /**< Short description for the data being saved in the new block */
  char * name /**< String with the name of the variable holding the data that will be writtne in the new block */
  )
{
  
  char type[] = "int";
  int type_size = sizeof (int);
  
  class_call ( binary_add_block (
                 file,
                 data_pointer,
                 size,
                 type_size,
                 description,
                 type,
                 name,
                 file->n_blocks),
    file->error_message,
    file->error_message);
  
  return _SUCCESS_;
  
}


/**
 * Add a data block made of long integers to the end of the binary file.
 *
 * This is a wrapper to binary_add_block(); see documentation there for
 * more details.
 */

int binary_append_long_int (
  struct binary_file * file, /**< Binary file to which the block will be appended */
  void * data_pointer, /**< Pointer to the data to be written in the new block */
  int size, /**< Number of elements to store in the block */
  char * description, /**< Short description for the data being saved in the new block */
  char * name /**< String with the name of the variable holding the data that will be writtne in the new block */
  )
{
  
  char type[] = "long int";
  int type_size = sizeof (long int);
  
  class_call ( binary_add_block (
                 file,
                 data_pointer,
                 size,
                 type_size,
                 description,
                 type,
                 name,
                 file->n_blocks),
    file->error_message,
    file->error_message);
  
  return _SUCCESS_;
  
}


/**
 * Add a data block made of double precision numbers to the end of the binary
 * file.
 *
 * This is a wrapper to binary_add_block(); see documentation there for
 * more details.
 */

int binary_append_double (
  struct binary_file * file, /**< Binary file to which the block will be appended */
  void * data_pointer, /**< Pointer to the data to be written in the new block */
  int size, /**< Number of elements to store in the block */
  char * description, /**< Short description for the data being saved in the new block */
  char * name /**< String with the name of the variable holding the data that will be writtne in the new block */
  )
{
  
  char type[] = "double";
  int type_size = sizeof (double);
  
  class_call ( binary_add_block (
                 file,
                 data_pointer,
                 size,
                 type_size,
                 description,
                 type,
                 name,
                 file->n_blocks),
    file->error_message,
    file->error_message);
  
  return _SUCCESS_;
  
}


/**
 * Add a data block made of double precision numbers to the end of the binary
 * file.
 *
 * This is a wrapper to binary_add_block(); see documentation there for
 * more details.
 */

int binary_append_string (
  struct binary_file * file, /**< Binary file to which the block will be appended */
  void * data_pointer, /**< Pointer to the data to be written in the new block */
  int size, /**< Number of elements to store in the block */
  char * description, /**< Short description for the data being saved in the new block */
  char * name /**< String with the name of the variable holding the data that will be writtne in the new block */
  )
{
  
  char type[] = "char";
  int type_size = sizeof (char);
  
  class_call ( binary_add_block (
                 file,
                 data_pointer,
                 size,
                 type_size,
                 description,
                 type,
                 name,
                 file->n_blocks),
    file->error_message,
    file->error_message);
  
  return _SUCCESS_;
  
}


/**
 * Add a data entry to an existing block in a binary file.
 *
 * The following objects will be updated:
 * - block->data_array
 * - block->size
 * - block->size_bytes
 * - file->size_bytes
 */

int binary_add_data (
  struct binary_file * file, /**< Binary file to which the data will be added */
  long int block_id, /**< Identifier of the data block to which the new data will be added */
  void * pointer, /**< Location of the new data to be added to the block */
  int size /**< Number of elements in the new data */
  ) 
{
  
  /* Shortcut to the block where the data belongs */
  struct binary_block * block = file->block_array[block_id];
  
  /* Add the data entry to the block */
  class_call (binary_add_data_to_block (
                block,
                pointer,
                size,
                file->error_message),
    file->error_message,
    file->error_message);
  
  /* Increase the size counter of the binary file */
  file->size_bytes += size*block->type_size;

  return _SUCCESS_;
  
}



/**
 * Add a data entry to the last block in a binary file.
 *
 * The following objects will be updated:
 * - block->data_array
 * - block->size
 * - block->size_bytes
 * - file->size_bytes
 */

int binary_append_data (
  struct binary_file * file, /**< Binary file to which the data will be added */
  void * pointer, /**< Location of the new data to be added to the block */
  int size /**< Number of elements in the new data */
  ) 
{
  
  class_call (binary_add_data (
                file,
                file->n_blocks-1,
                pointer,
                size),
    file->error_message,
    file->error_message);

  return _SUCCESS_;

}



/**
 * Add a data entry to a block.
 *
 * The following objects will be updated:
 * - block->data_array
 * - block->size
 * - block->size_bytes
 */

int binary_add_data_to_block (
  struct binary_block * block, /**< Data block to which the new data will be added */
  void * pointer, /**< Location of the new data to be added to the block */
  int size, /**< Number of elements in the new data */
  ErrorMsg errmsg /**< Area to write error messages */
  ) 
{

  /* Check that the pointer and size given make sense */
  class_test (pointer == NULL, errmsg, "found NULL pointer");
  class_test (size < 0, errmsg, "found negative size");
  
  /* Create the new data structure */
  struct binary_data * data;
  class_alloc (data, sizeof (struct binary_data), errmsg);
  
  /* Update the new data structure with location and size of the data */
  data->pointer = pointer;
  data->size = size;
  
  /* Increment the number of data objects in the block by one */ 
  block->n_data++;
  
  /* Extend the data array by one slot */
  class_realloc (
    block->data_array,
    block->data_array,
    block->n_data*sizeof (struct binary_data *),
    errmsg);

  /* Append the location & size of the data to the block */
  block->data_array[block->n_data-1] = data;
  
  /* Increase the size counter of the block */
  block->size += size;
  block->size_bytes += size * block->type_size;
  
  return _SUCCESS_;
  
}



/**
 * Change the size of a data entry in a block and update the parent
 * block and file accordingly.
 */

int binary_change_size (
  struct binary_file * file, /**< Binary file to which the block to modify belongs */
  long int block_id, /**< Indentifier of the block to modify */
  long int index_data, /**< Index of the data entry in the block */
  int new_size /**< New size (as in number of elements) for the data entry */
  ) 
{

  /* Check that the block exists in the file */
  class_test (block_id >= file->n_blocks,
    file->error_message,
    "the block you asked for (%d) does not exist", block_id);

  struct binary_block * block = file->block_array[block_id];

  class_test (new_size <= 0,
    file->error_message,
    "wrong size (%d) for data #%d (name=%s, id=%d)", new_size, index_data, block->name, block_id);

  /* Check that the data exists in the block */
  class_test (index_data >= block->n_data,
    file->error_message,
    "the data you asked for (%d) does not exist in the block #%d", index_data, block_id);

  struct binary_data * data = block->data_array[index_data];

  /* Update sizes */
  data->size = new_size;
  block->size += new_size;
  block->size_bytes += new_size * block->type_size;
  file->size_bytes += new_size * block->type_size;

  /* Shift the byte position of the other blocks to the right */
  for (int index_block=(block_id+1); index_block < file->n_blocks; ++index_block)
    file->block_array[index_block]->start_byte += new_size * block->type_size;
  
  return _SUCCESS_;
  
}



/**
 * Change the type of a block and update the parent file accordingly.
 *
 * The actual data, including the data pointers, are not changed. This
 * function is useful to downgrade certain blocks (double -> float or
 * int -> short) in order to save disk space.
 */

int binary_change_type (
  struct binary_file * file, /**< Binary file to which the block to modify belongs */
  long int block_id, /**< Indentifier of the block to modify */
  char * new_type, /**< String with the new type, eg. float */
  int new_type_size /**< Size in bytes of the new type, eg. sizeof(float) */
  )
{

  /* Check that the block exists in the file */
  class_test (block_id >= file->n_blocks,
    file->error_message,
    "the block you asked for (%d) does not exist", block_id);

  struct binary_block * block = file->block_array[block_id];

  class_test (new_type_size <= 0,
    file->error_message,
    "wrong type size (%d) for block #%d (name=%s)", new_type_size, block_id, block->name);

  /* Update the type */
  strcpy (block->type, new_type);
  block->type_size = new_type_size;

  /* Find the new size of the block */
  int old_size_bytes = block->size_bytes;
  int new_size_bytes = 0;
  for (int index_data=0; index_data < block->n_data; ++index_data)
    new_size_bytes += new_type_size*block->data_array[index_data]->size;

  /* Update the size in bytes of the block */
  block->size_bytes = new_size_bytes;

  /* Update the size in bytes of the parent file */
  file->size_bytes += new_size_bytes - old_size_bytes;

  /* Shift the byte position of the other blocks to the right */
  for (int index_block=(block_id+1); index_block < file->n_blocks; ++index_block)
    file->block_array[index_block]->start_byte += new_size_bytes - old_size_bytes;

  return _SUCCESS_;

}



/**
 * Delete a data block with a given ID and shift the other blocks accordingly
 */

int binary_delete_block (
  struct binary_file * file,
  long int block_position /**< Identifier of the block to delete */
  )
{
  
  /* Copy the size in bytes of the block to delete */
  int size_bytes = file->block_array[block_position]->size_bytes;
  
  /* Free the memory associated to the block to delete */
  free (file->block_array[block_position]);

  /* Modify the other blocks to account one less block */
  for (int block_id=block_position; block_id < (file->n_blocks-1); ++block_id) {

    /* Make room for the new block by shifting to the rigth the blocks with higher ID,
    starting from the last */
    file->block_array[block_id] = file->block_array[block_id+1];
    
    /* Update the block location in the binary file */
    file->block_array[block_id]->id -= 1;
    file->block_array[block_id]->start_byte -= size_bytes;
    
    class_test (file->block_array[block_id]->start_byte < 0,
      file->error_message,
      "error in deleting block");
    
  }

  /* Decrease by one the number of blocks in the binary file */
  file->n_blocks--;

  /* Reduce the size of block_array by one */
  class_realloc (
    file->block_array,
    file->block_array,
    file->n_blocks*sizeof(struct binary_block *),
    file->error_message
    );
  
  return _SUCCESS_;
  
}



/**
 * Build the binary map and append it to the header file.
 *
 * The binary map is a string table that informs the user of the location of the
 * various data blocks inside the binary file. It will be the last element
 * appearing in the header file.
 *
 * This function needs to be called after all blocks have been added to the
 * binary file, that is, after the last call of binary_add_block() and
 * before the call to binary_write_file().
 * 
 * TODO: This function adds lines to the header but needs to be called after
 * the header is added to the binary file. Then, how do we tell the binary
 * file the actual size of the header without relying on header_size_max?
 * 
 */

int binary_add_header_map (
  struct binary_file * file,
  int * map_size /**< Size of the binary map string */
  ) 
{

  /* Size of header before adding the binary map */
  int header_size_without_map = file->header_size;

  /* Add the map title to the header */
  binary_sprintf (file, "");
  binary_sprintf (file, "Binary map:");

  /* Find the longest description and name among the various blocks */
  int n_longest_description = 0;
  int n_longest_name = 0;
  for (int index_block=0; index_block < file->n_blocks; ++index_block) {
    struct binary_block * block = file->block_array[index_block];
    n_longest_description = MAX (n_longest_description, strlen (block->description));
    n_longest_name = MAX (n_longest_name, strlen (block->name));
  }
  
  /* Build the first line of the map with the keys */
  char format[_MAX_HEADER_LINE_LENGTH_];
  sprintf (format, "%%17s %%17s %%17s %%17s %%17s    %%-%ds %%-%ds", n_longest_name+2, n_longest_description+2);
  binary_sprintf (file, format, "BLOCK", "TYPE", "SIZE", "SIZE (BYTES)", "POS (BYTES)", "NAME", "DESCRIPTION");

  /* Build a line for each block with information on that block */
  for (int index_block=0; index_block < file->n_blocks; ++index_block) {

    struct binary_block * block = file->block_array[index_block];

    sprintf (format, "%%17d %%17s %%17d %%17d %%17d    %%-%ds %%-%ds", n_longest_name+2, n_longest_description+2);
    binary_sprintf (file, format, index_block, block->type, block->size, block->size*block->type_size,
      block->start_byte, block->name, block->description);

  }

  /* Length of the map */
  *map_size = file->header_size - header_size_without_map;


  return _SUCCESS_;

}



/**
 * Add a line to the header file.
 *
 * This function will prepend to the line a comment character (#) and
 * append a newline character (\n).
 */

int binary_add_header_line (
  struct binary_file * file, /**< Binary file containing the header file */
  char * line /**< String to be added to the header */
  )
{

  if (file->has_header == _TRUE_) {

    /* Truncate the line if it exceeds the maximum length */
    if (strlen(line) > _MAX_HEADER_LINE_LENGTH_)
      line[_MAX_HEADER_LINE_LENGTH_] = '\0';

    /* Number of characters that will be added to the header */
    int n_to_add = strlen (line) + 3;

    /* If adding the line will make the header grow over the allocated
    boundaries, do not add anything */
    if ((file->header_size+n_to_add) > file->header_size_max) {

      printf ("WARNING %s:%d: header reached maximum size (%d), will ignore further input\n",
        __FILE__, __LINE__, file->header_size_max);
      
      return _SUCCESS_;
    }

    /* Prepend a comment symbol to the line */
    char comment[4] = _COMMENT_;

    /* Append the line to the header */
    sprintf (file->header,
      "%s%s%s\n",
      file->header,
      _COMMENT_,
      line);

    /* Update the header size */
    file->header_size += n_to_add;

  }
  
  return _SUCCESS_;
  
}


/**
 * Add a line to the header file using the sprintf format.
 */

int binary_sprintf (
  struct binary_file * file, /**< Binary file containing the header file */
  const char * format, /**< Format for the sprintf function */
  ... /**< Arguments for the sprintf function */
  )
{
  
  /* Build the line by calling vsprintf on */
  char buffer[1024];
  va_list arg;
  va_start(arg, format);
  vsprintf(buffer, format, arg);
  va_end(arg);

  /* Add the line to the header */
  class_call (binary_add_header_line (
                file,
                buffer),
    file->error_message,
    file->error_message);

  return _SUCCESS_;
  
}


/**
 * Free the memory associated to the binary file.
 */

int binary_free (
  struct binary_file * file
  )
{
  
  if (file->has_header == _TRUE_)
    free (file->header);
  
  for (int index_block=0; index_block < file->n_blocks; ++index_block) {

    for (int index_data=0; index_data < file->block_array[index_block]->n_data; ++index_data)
      free (file->block_array[index_block]->data_array[index_data]);
    free (file->block_array[index_block]->data_array);    

    free (file->block_array[index_block]);
  }
  free (file->block_array);

  fclose (file->stream);

  free (file);

  return _SUCCESS_;

}

