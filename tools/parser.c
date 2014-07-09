#include "parser.h"

int parser_read_file(
		     char * filename,
		     struct file_content * pfc,
		     ErrorMsg errmsg
		     ){
  FILE * inputfile;
  char line[_LINE_LENGTH_MAX_];
  int counter;
  int is_data;
  FileArg name;
  FileArg value;

  class_open(inputfile,filename,"r",errmsg);

  counter = 0;
  while (fgets(line,_LINE_LENGTH_MAX_,inputfile) != NULL) {
    class_call(parser_read_line(line,&is_data,name,value,errmsg),errmsg,errmsg);
    if (is_data == _TRUE_) counter++;
  }

  class_test(counter == 0,
	     errmsg,
	     "No readable input in file %s",filename);

  class_call(parser_init(pfc,counter,filename,errmsg),
	     errmsg,
	     errmsg);

  rewind(inputfile);

  counter = 0;
  while (fgets(line,_LINE_LENGTH_MAX_,inputfile) != NULL) {
    class_call(parser_read_line(line,&is_data,name,value,errmsg),errmsg,errmsg);
    if (is_data == _TRUE_) {
      strcpy(pfc->name[counter],name);
      strcpy(pfc->value[counter],value);
      pfc->read[counter]=_FALSE_;
      counter++;
    }
  }

  fclose(inputfile);

  return _SUCCESS_;

}

int parser_init(
		struct file_content * pfc,
		int size,
    char * filename,
		ErrorMsg errmsg
		) {

  if (size > 0) {
    pfc->size=size;
    class_alloc(pfc->filename,(strlen(filename)+1)*sizeof(char),errmsg);
    strcpy(pfc->filename,filename);
    class_alloc(pfc->name,size*sizeof(FileArg),errmsg);
    class_alloc(pfc->value,size*sizeof(FileArg),errmsg);
    class_alloc(pfc->read,size*sizeof(short),errmsg);
    class_alloc(pfc->overwritten,size*sizeof(short),errmsg);
  }

  return _SUCCESS_;
}

int parser_free(
		struct file_content * pfc
		) {

  if (pfc->size > 0) {
    free(pfc->name);
    free(pfc->value);
    free(pfc->read);
    free(pfc->overwritten);
    free(pfc->filename);
  }

  return _SUCCESS_;
}

int parser_read_line(
		     char * line,
		     int * is_data,
		     char * name,
		     char * value,
		     ErrorMsg errmsg
		     ) {

  char * phash;
  char * pequal;
  char * left;
  char * right;

  /* check that there is an '=' */

  pequal=strchr(line,'=');
  if (pequal == NULL) {*is_data = _FALSE_; return _SUCCESS_;}

  /* if yes, check that there is not an '#' before the '=' */

  phash=strchr(line,'#');
  if ((phash != NULL) && (phash-pequal<2)) {*is_data = _FALSE_; return _SUCCESS_;}

  /* get the name, i.e. the block before the '=' */

  left=line;
  while (left[0]==' ') {
    left++;
  }

  right=pequal-1;
  while (right[0]==' ') {
    right--;
  }

  if (right-left < 0) {*is_data = _FALSE_; return _SUCCESS_;}

  class_test(right-left+1 >= _ARGUMENT_LENGTH_MAX_,
	     errmsg,
	     "name starting by '%s' too long; shorten it or increase _ARGUMENT_LENGTH_MAX_",strncpy(name,left,(_ARGUMENT_LENGTH_MAX_-1)));

  strncpy(name,left,right-left+1);
  name[right-left+1]='\0';

  /* get the value, i.e. the block after the '=' */

  left = pequal+1;
  while (left[0]==' ') {
    left++;
  }

  if (phash == NULL)
    right = line+strlen(line)-1;
  else
    right = phash-1;

  while (right[0]<=' ') {
    right--;
  }

  if (right-left < 0) {*is_data = _FALSE_; return _SUCCESS_;}

  class_test(right-left+1 >= _ARGUMENT_LENGTH_MAX_,
	     errmsg,
	     "value starting by '%s' too long; shorten it or increase _ARGUMENT_LENGTH_MAX_",strncpy(value,left,(_ARGUMENT_LENGTH_MAX_-1)));

  strncpy(value,left,right-left+1);
  value[right-left+1]='\0';

  *is_data = _TRUE_;

  return _SUCCESS_;

}

int parser_read_int(
		    struct file_content * pfc,
		    char * name,
		    int * value,
		    int * found,
		    ErrorMsg errmsg
		    ) {
  int index;
  int i;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* read parameter value. If this fails, return an error */

  class_test(sscanf(pfc->value[index],"%d",value) != 1,
	     errmsg,
	     "could not read value of parameter %s in file %s\n",name,pfc->filename);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_double(
		       struct file_content * pfc,
		       char * name,
		       double * value,
		       int * found,
		       ErrorMsg errmsg
		       ) {
  int index;
  int i;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* read parameter value. If this fails, return an error */

  class_test(sscanf(pfc->value[index],"%lg",value) != 1,
	     errmsg,
	     "could not read value of parameter %s in file %s\n",name,pfc->filename);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_double_and_position(
		       struct file_content * pfc,
		       char * name,
		       double * value,
               int * position,
		       int * found,
		       ErrorMsg errmsg
		       ) {
  int index;
  int i;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* read parameter value. If this fails, return an error */

  class_test(sscanf(pfc->value[index],"%lg",value) != 1,
	     errmsg,
	     "could not read value of parameter %s in file %s\n",name,pfc->filename);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  * position = index;

  return _SUCCESS_;

}

int parser_read_string(
		       struct file_content * pfc,
		       char * name,
		       FileArg * value,
		       int * found,
		       ErrorMsg errmsg
		       ) {
  int index;
  int i;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* read parameter value. */

  strcpy(*value,pfc->value[index]);

  /* Set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_list_of_doubles(
				struct file_content * pfc,
				char * name,
				int * size,
				double ** pointer_to_list,
				int * found,
				ErrorMsg errmsg
				) {
  int index;
  int i;

  char * string;
  char * substring;
  FileArg string_with_one_value;

  double * list;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* count number of comas and compute size = number of comas + 1 */
  i = 0;
  string = trim (pfc->value[index], ',');   /* Remove leading and trailing commas that would otherwise mess with the parsing */ 
  do {
    i ++;
    substring = strchr(string,',');
    string = substring+1;
  } while(substring != NULL);

  *size = i;

  /* free and re-allocate array of values */
  class_alloc(list,*size*sizeof(double),errmsg);
  *pointer_to_list = list;

  /* read one double between each comas */
  i = 0;
  string = trim (pfc->value[index], ',');   /* Remove leading and trailing commas that would otherwise mess with the parsing */ 
  do {
    i ++;
    substring = strchr(string,',');
    if (substring == NULL) {
      strcpy(string_with_one_value,string);
    }
    else {
      strncpy(string_with_one_value,string,(substring-string));
      string_with_one_value[substring-string]='\0';
    }
    class_test(sscanf(string_with_one_value,"%lg",&(list[i-1])) != 1,
	       errmsg,
	       "could not read %dth value of list of parameters %s in file %s\n",
	       i,
	       name,
	       pfc->filename);
    string = substring+1;
  } while(substring != NULL);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_list_of_integers(
				 struct file_content * pfc,
				 char * name,
				 int * size,
				 int ** pointer_to_list,
				 int * found,
				 ErrorMsg errmsg
				 ) {
  int index;
  int i;

  char * string;
  char * substring;
  FileArg string_with_one_value;

  int * list;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* count number of comas and compute size = number of comas + 1 */
  i = 0;
  string = trim (pfc->value[index], ',');   /* Remove leading and trailing commas that would otherwise mess with the parsing */ 
  do {
    i ++;
    substring = strchr(string,',');
    string = substring+1;
  } while(substring != NULL);

  *size = i;

  /* free and re-allocate array of values */
  class_alloc(list,*size*sizeof(int),errmsg);
  *pointer_to_list = list;

  /* read one integer between each comas */
  i = 0;
  string = trim (pfc->value[index], ',');   /* Remove leading and trailing commas that would otherwise mess with the parsing */ 
  do {
    i ++;
    substring = strchr(string,',');
    if (substring == NULL) {
      strcpy(string_with_one_value,string);
    }
    else {
      strncpy(string_with_one_value,string,(substring-string));
      string_with_one_value[substring-string]='\0';
    }
    class_test(sscanf(string_with_one_value,"%d",&(list[i-1])) != 1,
	       errmsg,
	       "could not read %dth value of list of parameters %s in file %s\n",
	       i,
	       name,
	       pfc->filename);
    string = substring+1;
  } while(substring != NULL);

  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */

  * found = _TRUE_;
  pfc->read[index] = _TRUE_;

  /* check for multiple entries of the same parameter. If another occurence is found,
     return an error. */

  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }

  /* if everything proceeded normally, return with 'found' flag equal to true */

  return _SUCCESS_;

}

int parser_read_list_of_strings(
				struct file_content * pfc,
				char * name,
				int * size,
				char ** pointer_to_list,
				int * found,
				ErrorMsg errmsg
				) {
  int index;
  int i;

  char * string;
  char * substring;
  FileArg string_with_one_value;

  char * list;

  /* intialize the 'found' flag to false */

  * found = _FALSE_;

  /* search parameter */

  index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  /* if parameter not found, return with 'found' flag still equal to false */

  if (index == pfc->size)
    return _SUCCESS_;

  /* count number of comas and compute size = number of comas + 1 */
  i = 0;
  string = trim (pfc->value[index], ',');   /* Remove leading and trailing commas that would otherwise mess with the parsing */ 
  do {
    i ++;
    substring = strchr(string,',');
    string = substring+1;
  } while(substring != NULL);

  *size = i;

  /* free and re-allocate array of values */
  class_alloc(list,*size*sizeof(FileArg),errmsg);
  *pointer_to_list = list;

  /* read one string between each comas */
  i = 0;
  string = trim (pfc->value[index], ',');   /* Remove leading and trailing commas that would otherwise mess with the parsing */ 
  do {
    i ++;
    substring = strchr(string,',');
    if (substring == NULL) {
      strcpy(string_with_one_value,string);
    }
    else {
      strncpy(string_with_one_value,string,(substring-string));
      string_with_one_value[substring-string]='\0';
    }
    strcpy(list+(i-1)*_ARGUMENT_LENGTH_MAX_,string_with_one_value);
    //Insert EOL character:
    *(list+i*_ARGUMENT_LENGTH_MAX_-1) = '\n';
    string = substring+1;
  } while(substring != NULL);
  /* if parameter read correctly, set 'found' flag to true, as well as the flag
     associated with this parameter in the file_content structure */
  * found = _TRUE_;
  pfc->read[index] = _TRUE_;
  /* check for multiple entries of the same parameter. If another occurence is
     found,
     return an error. */
  for (i=index+1; i < pfc->size; i++) {
    class_test(strcmp(pfc->name[i],name) == 0,
	       errmsg,
	       "multiple entry of parameter %s in file %s\n",name,pfc->filename);
  }
  /* if everything proceeded normally, return with 'found' flag equal to true */
  return _SUCCESS_;
}

int parser_cat(
	       struct file_content * pfc1,
	       struct file_content * pfc2,
	       struct file_content * pfc3,
	       ErrorMsg errmsg
	       ) {

  int i;

  class_test(pfc1->size < 0.,
	     errmsg,
	     "size of file_content structure probably not initialized properly\n");

  class_test(pfc2->size < 0.,
	     errmsg,
	     "size of file_content structure probably not initialized properly\n");

  if (pfc1->size == 0) {
    class_alloc(pfc3->filename,(strlen(pfc2->filename)+1)*sizeof(char),errmsg);
    sprintf(pfc3->filename,"%s",pfc2->filename);
  }
  if (pfc2->size == 0) {
    class_alloc(pfc3->filename,(strlen(pfc1->filename)+1)*sizeof(char),errmsg);
    sprintf(pfc3->filename,"%s",pfc1->filename);
  }
  if ((pfc1->size !=0) && (pfc2->size != 0)) {
    class_alloc(pfc3->filename,(strlen(pfc1->filename)+strlen(pfc2->filename)+5)*sizeof(char),errmsg);
    sprintf(pfc3->filename,"%s or %s",pfc1->filename,pfc2->filename);
  }

  pfc3->size = pfc1->size + pfc2->size;
  class_alloc(pfc3->value,pfc3->size*sizeof(FileArg),errmsg);
  class_alloc(pfc3->name,pfc3->size*sizeof(FileArg),errmsg);
  class_alloc(pfc3->read,pfc3->size*sizeof(short),errmsg);

  for (i=0; i < pfc1->size; i++) {
    strcpy(pfc3->value[i],pfc1->value[i]);
    strcpy(pfc3->name[i],pfc1->name[i]);
    pfc3->read[i]=pfc1->read[i];
  }

  for (i=0; i < pfc2->size; i++) {
    strcpy(pfc3->value[i+pfc1->size],pfc2->value[i]);
    strcpy(pfc3->name[i+pfc1->size],pfc2->name[i]);
    pfc3->read[i+pfc1->size]=pfc2->read[i];
  }

  return _SUCCESS_;

}

/**
 * Add one entry to the input file_content structure. 
 *
 * If the entry already exists, store its index in 'index' and act
 * accordingly to the option 'what_to_do', which can be one of the
 * following:
 * - REPLACE: overwrite entry with the provided value.
 * - APPEND: append to the existing value the provided one.
 * - RAISE: returns _FAILURE_ with an error message
 *
 * If the entry does not exist, the entry will be created, pfc->size
 * will be incremented by one, and 'index' will refer to the position
 * of the new entry.
 */

int parser_add_entry (
        struct file_content * pfc,
        char * name,
        char * value,
        enum entry_operation what_to_do,
        int * index,
        ErrorMsg errmsg
        )
{

  /* Make sure the input string is not too long */
  class_test (strlen(value)>_ARGUMENT_LENGTH_MAX_, errmsg, "string is too long");

  /* Search for the parameter */
  *index=0;
  while ((*index < pfc->size) && (strcmp(pfc->name[*index],name) != 0))
    (*index)++;
  int found = (*index < pfc->size);

  if (!found) {

    /* If the entry is NOT found, create it and increment the number of entries */      
    strcpy (pfc->name[pfc->size], name);
    strcpy (pfc->value[pfc->size], value);
    pfc->size++;
    return _SUCCESS_;
  }
  else {

    /* If the entry is found, act according to the variable 'what_to_do' */      
    switch (what_to_do) {
      
      case REPLACE:
        strcpy (pfc->value[*index], value);
        pfc->overwritten[*index] = _TRUE_;
        return _SUCCESS_;
        
      case APPEND:
        class_test ((strlen(value)+strlen(pfc->value[*index]))>_ARGUMENT_LENGTH_MAX_, errmsg,
          "string for name='%s' is too long: original=%d, append=%d, max_length=%d",
          name, strlen(pfc->value[*index]), strlen(value), _ARGUMENT_LENGTH_MAX_);
        strcat (pfc->value[*index], value);
        pfc->overwritten[*index] = _TRUE_;
        return _SUCCESS_;

      case RAISE:
        class_stop (errmsg, "parameter '%s' already present in input file structure", name);
        
      default:
        class_stop (errmsg, "option '%d' not contemplated", name);
    }
  }

  class_stop (errmsg, "should not be here");

}


/**
 * OBSOLETE, USE 'parser_add_entry' INSTEAD!
 *
 * Add one entry to the input file_content structure.
 *
 * This function has two behaviours according to the value passed for 'found'.
 * If 'found' is a NULL pointer and the entry corresponding to 'name' already exists,
 * print an error message and return _FAILURE_. 
 * If 'found' is a valid pointer to 'int', overwrite it with either _TRUE_ or _FALSE_
 * whether the entry was found or not. If the entry is found, overwrite it with
 * 'value'.
 */
int parser_create_or_replace_entry (
        struct file_content * pfc,
        char * name,
        char * value,
        int * found,
        ErrorMsg errmsg
        )
{

  /* make sure the input string is not too long */
  class_test (strlen(value)>_ARGUMENT_LENGTH_MAX_, errmsg, "string is too long");

  /* search parameter */

  int index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  if (found == NULL) {

    /* if the entry exists and found==NULL, return error */
    class_test (index < pfc->size,
      errmsg,
      "parameter '%s' already present in input file structure", name);

  }
  else {

    if (index < pfc->size) {

      /* if the entry exists and found!=NULL, overwrite parameter value */
      *found = _TRUE_;
      strcpy (pfc->value[index], value);
      pfc->overwritten[index] = _TRUE_;
      return _SUCCESS_;
    }
    else {
      *found = _FALSE_;
    }
  }
    
  /* at this point, we must have index==pfc->size */
  class_test (index!=pfc->size, errmsg, "should not happen");
      
  /* create the new entry and increment the number of parameters */      
  strcpy (pfc->name[pfc->size], name);
  strcpy (pfc->value[pfc->size], value);
  pfc->size++;

  /* if everything proceeded normally, return _SUCCESS_ */
  return _SUCCESS_;

}

/**
 * OBSOLETE, USE 'parser_add_entry' INSTEAD!
 *
 * Add one entry to the input file_content structure.
 *
 * This function has two behaviours according to the value passed for 'found'.
 * If 'found' is a NULL pointer and the entry corresponding to 'name' already exists,
 * print an error message and return _FAILURE_. 
 * If 'found' is a valid pointer to 'int', overwrite it with either _TRUE_ or _FALSE_
 * whether the entry was found or not. If the entry is found, append the string 'value'
 * to it.
 */
int parser_create_or_append_to_entry (
        struct file_content * pfc,
        char * name,
        char * value,
        int * found,
        ErrorMsg errmsg
        )
{

  /* make sure the input string is not too long */
  class_test (strlen(value)>_ARGUMENT_LENGTH_MAX_, errmsg, "string is too long");

  /* search parameter */
  int index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  if (found == NULL) {

    /* if the entry exists and found==NULL, return error */
    class_test (index < pfc->size,
      errmsg,
      "parameter '%s' already present in input file structure", name);

  }
  else {

    if (index < pfc->size) {

      /* if the entry exists and found!=NULL, append to the current parameter value */
      *found = _TRUE_;
      class_test ((strlen(value)+strlen(pfc->value[index]))>_ARGUMENT_LENGTH_MAX_, errmsg,
        "string for name='%s' is too long: original=%d, append=%d, max_length=%d",
        name, strlen(pfc->value[index]), strlen(value), _ARGUMENT_LENGTH_MAX_);
      strcat (pfc->value[index], value);
      pfc->overwritten[index] = _TRUE_;
      return _SUCCESS_;
    }
    else {
      *found = _FALSE_;
    }
  }
    
  /* at this point, we must have index==pfc->size */
  class_test (index!=pfc->size, errmsg, "should not happen");
      
  /* create the new entry and increment the number of parameters */      
  strcpy (pfc->name[pfc->size], name);
  strcpy (pfc->value[pfc->size], value);
  pfc->size++;

  /* if everything proceeded normally, return _SUCCESS_ */
  return _SUCCESS_;

}

// /**
//  * Add text to a parameter entry in the input file_content structure.
//  *
//  * If 'found' is a NULL pointer and the entry corresponding to 'name' does not exist,
//  * print an error message. Otherwise, overwrite 'found' with either _TRUE_ or _FALSE_
//  * whether the entry was found or not, with no error messages. If it is found, append
//  * 'append_value' to the entry's value, otherwise create it from scratch.
//  *
//  */
// int parser_create_or_append_entry (
//         struct file_content * pfc,
//         char * name,
//         char * append_value,
//         int * found,
//         ErrorMsg errmsg
//         )
// {
//
//   /* search parameter */
//
//   int index=0;
//   while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
//     index++;
//
//   if (found == NULL) {
//
//     /* if the entry does not exist and found==NULL, return error */
//     class_test (index == pfc->size,
//       errmsg,
//       "parameter '%s' not present in input file structure", name);
//
//   }
//   else {
//
//     if (index == pfc->size) {
//
//       /* if the entry does not exist, create it */
//       *found = _FALSE_;
//       strcpy (pfc->name[pfc->size], name);
//       strcpy (pfc->value[pfc->size], append_value);
//       pfc->size++;
//       return _SUCCESS_;
//     }
//     else {
//
//       *found = _TRUE_;
//     }
//   }
//
//   /* at this point, the entry must exist */
//   class_test (index >= pfc->size, errmsg, "should not happen");
//
//   /* append 'append_value' to the entry */
//   strcat (pfc->value[pfc->size], append_value);
//
//   /* if everything proceeded normally, return _SUCCESS_ */
//
//   return _SUCCESS_;
//
// }


/**
 * Modify one entry of the input file_content structure.
 *
 * This function has two behaviours according to the value passed for 'found'.
 * If 'found' is a NULL pointer and the entry corresponding to 'name' does no exist,
 * print an error message and return _FAILURE_. 
 * If 'found' is a valid pointer to 'int', overwrite it with either _TRUE_ or _FALSE_
 * whether the entry was found or not. If the entry is found, overwrite it with
 * 'value'.
 */
int parser_overwrite_entry (
        struct file_content * pfc,
        char * name,
        char * new_value,
        int * found,
        ErrorMsg errmsg
        )
{

  /* make sure the input string is not too long */
  class_test (strlen(new_value)>_ARGUMENT_LENGTH_MAX_, errmsg, "string is too long");

  /* search parameter */
  int index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;

  if (found == NULL) {
    class_test (index == pfc->size,
      errmsg,
      "parameter '%s' not found in input file structure", name);
  }
  else {
    if (index < pfc->size) {
      *found = _TRUE_;
    }
    else {
      *found = _FALSE_;
      return _SUCCESS_;
    }
  }

  /* overwrite parameter value if found. */
  strcpy (pfc->value[index], new_value);

  /* if parameter overwritten correctly, set the flag 
  associated with this parameter in the file_content structure */ 
  pfc->overwritten[index] = _TRUE_;

  /* if everything proceeded normally, return _SUCCESS_ */
  return _SUCCESS_;

}


/**
 * Remove an entry from the file_content structure.
 *
 * This function has two behaviours according to the value passed for 'found'.
 * If 'found' is a NULL pointer and the entry corresponding to 'name' does not exist,
 * print an error message and return _FAILURE_. 
 * If 'found' is a valid pointer to 'int', overwrite it with either _TRUE_ or _FALSE_
 * whether the entry was found or not. If the entry is found, replace its name with an
 * empty string, so that it cannot be found by any of the parser functions.
 */
int parser_remove_entry (
        struct file_content * pfc,
        char * name,
        int * found,
        ErrorMsg errmsg
        )
{

  /* search parameter */
  int index=0;
  while ((index < pfc->size) && (strcmp(pfc->name[index],name) != 0))
    index++;
  
  if (found == NULL) {
    class_test (index == pfc->size,
      errmsg,
      "parameter '%s' not found in input file structure", name);
  }
  else {
    if (index < pfc->size) {
      *found = _TRUE_;
    }
    else {
      *found = _FALSE_;
      return _SUCCESS_;
    }
  }

  /* removing the entry is equivalent to setting its name to an empty string,
  so that it will never be found. */
  strcpy (pfc->name[index], "");

  /* if parameter overwritten correctly, set the flag 
  associated with this parameter in the file_content structure */ 
  pfc->overwritten[index] = _TRUE_;

  /* if everything proceeded normally, return _SUCCESS_ */
  return _SUCCESS_;

}

/** 
 * Print the content of the parameter structure
 */
int parser_print (
        struct file_content * pfc
        )
{

  for (int i=0; i < pfc->size; ++i) {
    printf ("parameter[%25s] = %s\n", pfc->name[i], pfc->value[i]);
  }
  
  return _SUCCESS_;

}

/** 
 * Remove leading and trailing commmas from the parser values.
 */
int parser_remove_extra_commas (
        struct file_content * pfc
        )
{

  for (int i=0; i < pfc->size; ++i)
    strcpy (pfc->value[i], trim(((&(pfc->value[i])[0])), ','));
  
  return _SUCCESS_;

}

/** 
 * Remove from a string the leading and trailing characters 'c', if it is there.
 * Copied from http://stackoverflow.com/a/122974/2972183, credits to user indiv.
 */
char * trim(char *str, char c)
{
    size_t len = 0;
    char *frontp = str - 1;
    char *endp = NULL;

    if( str == NULL )
            return NULL;

    if( str[0] == '\0' )
            return str;

    len = strlen(str);
    endp = str + len;

    /* Move the front and back pointers to address
     * the first non-whitespace characters from
     * each end.
     */
    while( (*(++frontp)) == c );
    while( ((*(--endp)) == c) && endp != frontp );
    // while( isspace(*(++frontp)) );
    // while( isspace(*(--endp)) && endp != frontp );

    if( str + len - 1 != endp )
            *(endp + 1) = '\0';
    else if( frontp != str &&  endp == frontp )
            *str = '\0';

    /* Shift the string so that it starts at str so
     * that if it's dynamically allocated, we can
     * still free it on the returned pointer.  Note
     * the reuse of endp to mean the front of the
     * string buffer now.
     */
    endp = str;
    if( frontp != str )
    {
            while( *frontp ) *endp++ = *frontp++;
            *endp = '\0';
    }

    return str;
}
