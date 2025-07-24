#ifndef FHANDLE_FILE
#define FHANDLE_FILE

#include <zlib.h>
#include "../SEutils.h" // for safe_malloc

/******************
 * Functions to abstract file access interfaces
 * Note that no function is provided to replace truncate_file
 * If this is necessary it will need to be implemented
 ******************/


/*
// not currently supported
#ifdef _WIN32
  #include <bzlib.h>
#endif
*/

enum {
  UNKNOWN = 0,
  UNCOMPRESSED,
  GZ_TYPE,
  BZ2_TYPE,
  XZ_TYPE,
};

typedef struct file_t {
  void *file_ptr;
  int file_type;
} file_t;

static const int MAX_READ_RETRIES = 10;
static const int MAX_WRITE_RETRIES = 10;

/*** externally exposed methods ***/
file_t *safe_fopen(const char* fpath, const char* mode, const int type);
void safe_fclose(file_t *f);
size_t safe_fread(void *buffer, size_t size, size_t count, file_t *stream);
size_t unsafe_fread(void *buffer, size_t size, size_t count, file_t *stream);
size_t safe_fwrite(void *buffer, size_t size, size_t count, file_t *stream);
char safe_getc(file_t *f);
int safe_fseek(file_t *f, long offset, int origin);
void safe_rewind(file_t *f);
long safe_ftell(file_t *f);
int safe_feof(file_t *f);


#endif