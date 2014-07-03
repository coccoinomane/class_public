#ifndef __PARSER__
#define __PARSER__

#include "common.h"

#define _LINE_LENGTH_MAX_ 1024 /**< size of the string read in each line of the file (extra characters not taken into account) */
#define _ARGUMENT_LENGTH_MAX_ 1024 /**< maximum size of each argument (name or value), including the final null character */

typedef char FileArg[_ARGUMENT_LENGTH_MAX_];

/* Option for the function 'parser_add_entry' */
enum entry_operation {
  REPLACE = 0,
  APPEND = 1,
  RAISE = 2
};

/* after reading a given file, all relevant information stored in this structure, in view of being processed later*/
struct file_content {
  char * filename;
  int size;
  FileArg * name;  /**< list of (size) names */
  FileArg * value; /**< list of (size) values */
  short * read;    /**< set to _TRUE_ if this parameter is effectively read */
  short * overwritten;    /**< set to _TRUE_ if this parameter has been overwritten  */
};

/**************************************************************/

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int parser_read_file(
		     char * filename,
		     struct file_content * pfc,
		     ErrorMsg errmsg
		     );

int parser_init(
		struct file_content * pfc,
		int size,
        char * filename,
		ErrorMsg errmsg
		);

int parser_free(
		struct file_content * pfc
		);

int parser_read_line(
		char * line,
		int * is_data,
		char * name,
		char * value,
		ErrorMsg errmsg
		);

int parser_read_int(
		    struct file_content * pfc,
		    char * name,
		    int * value,
		    int * found,
		    ErrorMsg errmsg
		    );

int parser_read_double(
		    struct file_content * pfc,
		    char * name,
		    double * value,
		    int * found,
		    ErrorMsg errmsg
		    );

  int parser_read_double_and_position(
                                      struct file_content * pfc,
                                      char * name,
                                      double * value,
                                      int * position,
                                      int * found,
                                      ErrorMsg errmsg
                                      );

int parser_read_string(
		       struct file_content * pfc,
		       char * name,
		       FileArg * value,
		       int * found,
		       ErrorMsg errmsg
		       );

int parser_read_list_of_doubles(
				struct file_content * pfc,
				char * name,
				int * size,
				double ** pointer_to_list,
				int * found,
				ErrorMsg errmsg
				);

int parser_read_list_of_integers(
				struct file_content * pfc,
				char * name,
				int * size,
				int ** pointer_to_list,
				int * found,
				ErrorMsg errmsg
				);

int parser_read_list_of_strings(
				struct file_content * pfc,
				char * name,
				int * size,
				char ** pointer_to_list,
				int * found,
				ErrorMsg errmsg
				);

int parser_cat(
	       struct file_content * pfc1,
	       struct file_content * pfc2,
	       struct file_content * pfc3,
	       ErrorMsg errmsg
	       );
         
int parser_add_entry (
        struct file_content * pfc,
        char * name,
        char * value,
        enum entry_operation what_to_do,
        int * index,
        ErrorMsg errmsg
        );
         
int parser_create_or_replace_entry (
        struct file_content * pfc,
        char * name,
        char * value,
        int * found,
        ErrorMsg errmsg
        );
        
int parser_create_or_append_to_entry (
        struct file_content * pfc,
        char * name,
        char * value,
        int * found,
        ErrorMsg errmsg
        );
                 
int parser_overwrite_entry (
		    struct file_content * pfc,
		    char * name,
		    char * new_value,
        int * found,
		    ErrorMsg errmsg
		    );

int parser_remove_entry (
        struct file_content * pfc,
        char * name,
        int * found,
        ErrorMsg errmsg
        );

int parser_print (
        struct file_content * pfc
        );

#ifdef __cplusplus
}
#endif

#endif
