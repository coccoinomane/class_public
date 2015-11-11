#include "common.h"

void class_protect_sprintf(char* dest, char* tpl,...) {
  va_list args;
  va_start(args,tpl);
  vsnprintf(dest, 2048,tpl,args);
  va_end(args);
}

void class_protect_fprintf(FILE* stream, char* tpl,...) {
  va_list args;
  char dest[6000];
  va_start(args,tpl);
  vsnprintf(dest, 2048,tpl,args);
  va_end(args);
  fprintf(stream,dest);
}

void* class_protect_memcpy(void* dest, void* from, size_t sz) {
  return memcpy(dest, from,sz);
}

int get_number_of_titles(char * titlestring){
  int i;
  int number_of_titles=0;

  for (i=0; i<strlen(titlestring); i++){
    if (titlestring[i] == '\t')
      number_of_titles++;
  }
  return number_of_titles;
}

#ifdef WITH_SONG1

/* Initialise the file pointer to the log file */
FILE * log_file = NULL;

/**
 * Print the same message to two streams.
 *
 * The first stream is used only if verbose_level>min1, the second only if
 * verbose_level>min2. 
 * 
 * This function can be used to simultaneously print a message to screen
 * (stream1=stdout) and to an already opened log file (stream2=FILE_pointer).
 * In this case, using min1<min2 would record more information than what is
 * printed to screen.
 */
void fprintf_2way(
  int verbose_level,
  FILE *stream1,
  int min1,
  FILE *stream2,
  int min2,
  char *fmt, ...)
{
  va_list argp;

  if ((verbose_level > min1) && (stream1 != NULL)) {
    va_start(argp, fmt);
    vfprintf(stream1, fmt, argp);
    va_end(argp);
  }
    
  if ((verbose_level > min2) && (stream2 != NULL)) {
    va_start(argp, fmt);
    vfprintf(stream2, fmt, argp);
    va_end(argp);
  }
}


/**
 * Print the same message to three streams.
 *
 * The first stream is used only if verbose_level>min1, the second only if
 * verbose_level>min2 and the third only if verbose_level>min3. 
 */
void fprintf_3way(
  int verbose_level,
  FILE *stream1,
  int min1,
  FILE *stream2,
  int min2,
  FILE *stream3,
  int min3,
  char *fmt, ...)
{
  va_list argp;

  if ((verbose_level > min1) && (stream1 != NULL)) {
    va_start(argp, fmt);
    vfprintf(stream1, fmt, argp);
    va_end(argp);
  }
    
  if ((verbose_level > min2) && (stream2 != NULL)) {
    va_start(argp, fmt);
    vfprintf(stream2, fmt, argp);
    va_end(argp);
  }

  if ((verbose_level > min3) && (stream3 != NULL)) {
    va_start(argp, fmt);
    vfprintf(stream3, fmt, argp);
    va_end(argp);
  }
}

#endif // WITH_SONG1





