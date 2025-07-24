#include "FileHandlers.h"

static void detect_file_type(const char *expath, int *ztype, int *subtype){
  /*
   * adapted from XVector/io_utils.c, which itself is taken from
   * do_url() in R/src/main/connections.c
   */
  FILE *fp;
  char buf[7];
  int res;

  *ztype = UNKNOWN;
  *subtype = 0;
  if ((fp = fopen(expath, "rb")) == NULL)
    return;
  memset(buf, 0, 7);
  res = fread(buf, 5, 1, fp);
  fclose(fp);
  if (res != 1)
    return;
  if (buf[0] == '\x1f' && buf[1] == '\x8b')
    *ztype = GZ_TYPE;
  else if (strncmp(buf, "BZh", 3) == 0)
    *ztype = BZ2_TYPE;
  else if (buf[0] == '\xFD' && strncmp(buf+1, "7zXZ", 4) == 0)
    *ztype = XZ_TYPE;
  else if ((buf[0] == '\xFF') && strncmp(buf+1, "LZMA", 4) == 0) {
    *ztype = XZ_TYPE;
    *subtype = 1;
  } else if (memcmp(buf, "]\0\0\200\0", 5) == 0) {
    *ztype = XZ_TYPE;
    *subtype = 1;
  } else {
    *ztype = UNCOMPRESSED;
  }
  return;
}

static void fopen_with_known_type(file_t *f, const char* path, const char *mode){
  int ftype = f->file_type;
  switch(ftype){
  case UNCOMPRESSED:
    f->file_ptr = fopen(path, mode);
    break;
  case GZ_TYPE:
    f->file_ptr = gzopen(path, mode);
    break;
  case BZ2_TYPE:
    error("Files compressed with bz2 are not yet supported (%s)", path);
  case XZ_TYPE:
    error("Files compressed with xz are not yet supported (%s)", path);
  case UNKNOWN:
  default:
    error("Internal error calling fopen on file %s in mode %s.", path, mode);
  }
}

static file_t *fopen_dispatch(const char* fpath, const char* mode, const int type){
  file_t *file_holder = safe_malloc(sizeof(file_t));
  file_holder->file_ptr = NULL;
  file_holder->file_type = UNKNOWN;

  // determine the file type
  int ftype = type, fsubtype = 0;
  if(ftype == UNKNOWN)
    detect_file_type(fpath, &ftype, &fsubtype);

  // subtype is currently ignored but may be used in the future if we add xz support
  file_holder->file_type = ftype;

  // open the file given that we now know its type
  fopen_with_known_type(file_holder, fpath, mode);

  return file_holder;
}

static int fclose_dispatch(file_t *f){
  int ftype = f->file_type;
  switch(ftype){
  case UNCOMPRESSED:
    return fclose(f->file_ptr);
  case GZ_TYPE:
    return gzclose(f->file_ptr);
  case BZ2_TYPE:
  case XZ_TYPE:
  case UNKNOWN:
  default:
    // should never get here, since we can never open files in this way
    error("Internal error attempting to call fclose on file with format %d",
      f->file_type);
  }
};

static size_t fread_dispatch(void *buf, size_t size, size_t count, file_t *f){
  int ftype = f->file_type;
  size_t bytes_read = 0;
  switch(ftype){
  case UNCOMPRESSED:
    bytes_read = fread(buf, size, count, f->file_ptr);
    break;
  case GZ_TYPE:
    bytes_read = gzread(f->file_ptr, buf, size*count);
    break;
  case BZ2_TYPE:
  case XZ_TYPE:
  case UNKNOWN:
  default:
    // should never get here, since we can never open files in this way
    error("Internal error attempting to call fread on file with format %d",
      f->file_type);
  }
  return bytes_read;
}

static size_t fwrite_dispatch(const void* buf, size_t size, size_t count, file_t *f){
  int ftype = f->file_type;
  size_t bytes_written = 0;
  switch(ftype){
  case UNCOMPRESSED:
    bytes_written = fwrite(buf, size, count, f->file_ptr);
    break;
  case GZ_TYPE:
    bytes_written = gzwrite(f->file_ptr, buf, size*count);
    break;
  case BZ2_TYPE:
  case XZ_TYPE:
  case UNKNOWN:
  default:
    // should never get here, since we can never open files in this way
    error("Internal error attempting to call fwrite on file with format %d",
      f->file_type);
  }
  return bytes_written;
}

static char getc_dispatch(file_t *f){
  int ftype = f->file_type;
  switch(ftype){
  case UNCOMPRESSED:
    return getc(f->file_ptr);
  case GZ_TYPE:
    return gzgetc((gzFile) f->file_ptr);
  case BZ2_TYPE:
  case XZ_TYPE:
  case UNKNOWN:
  default:
    // should never get here, since we can never open files in this way
    error("Internal error attempting to call getc on file with format %d",
      f->file_type);
  }
  return 0;
}

static int fseek_dispatch(file_t *f, long offset, int origin){
  int ftype = f->file_type;
  int status = -1;
  switch(ftype){
  case UNCOMPRESSED:
    status = fseek(f->file_ptr, offset, origin);
    break;
  case GZ_TYPE:
    // note that SEEK_END may not be supported for gz files
    status = gzseek(f->file_ptr, offset, origin);
    break;
  case BZ2_TYPE:
  case XZ_TYPE:
  case UNKNOWN:
  default:
    // should never get here, since we can never open files in this way
    error("Internal error attempting to seek in file");
  }
  return status;
}

static long ftell_dispatch(file_t *f){
  int ftype = f->file_type;
  switch(ftype){
  case UNCOMPRESSED:
    return ftell(f->file_ptr);
  case GZ_TYPE:
    return gztell(f->file_ptr);
  case BZ2_TYPE:
  case XZ_TYPE:
  case UNKNOWN:
  default:
    // should never get here, since we can never open files in this way
    error("Internal error attempting to call getc() on file");
  }
  return -1;
}

static int feof_dispatch(file_t *f){
  int ftype = f->file_type;
  switch(ftype){
  case UNCOMPRESSED:
    return feof(f->file_ptr);
  case GZ_TYPE:
    return gzeof(f->file_ptr);
  case BZ2_TYPE:
  case XZ_TYPE:
  case UNKNOWN:
  default:
    // should never get here, since we can never open files in this way
    error("Internal error attempting to call feof() on file");
  }
  // default should be saying the file is closed to abort as early as possible
  return 1;
}

/***************************/
/* File Function Interface */
/***************************/

file_t *safe_fopen(const char* fpath, const char* mode, const int type){
  file_t *f = fopen_dispatch(fpath, mode, type);
  if(!f->file_ptr){
    free(f);
    error("Internal error calling fopen on file %s (Failed to open).", fpath);
  }
  return f;
}

void safe_fclose(file_t *f){
  if(!f) error("Internal error calling fclose on null pointer.");

  int status = fclose_dispatch(f);
  if(status != 0)
    error("Internal error calling fclose on file with format %d.", f->file_type);
  free(f);
  f = NULL;
  return;
}

size_t safe_fread(void *buffer, size_t size, size_t count, file_t *stream){
  size_t found_values = fread_dispatch(buffer, size, count, stream);
  if(found_values != count){
    // two scenarios:

    // 1. read past the end of the file (throw error and return)
    if(safe_feof(stream))
      error("%s", "Internal error: fread tried to read past end of file.\n");

    // 2. some undefined reading error (retry a few times and then return)
    for(int i=0; i<MAX_READ_RETRIES; i++){
      // if we read a partial value, reset the counter back
      if(found_values) fseek_dispatch(stream, -1*((int)found_values), SEEK_CUR);

      // try to read again
      found_values = fread_dispatch(buffer, size, count, stream);
      if(found_values == count) return found_values;
    }

    // otherwise throw an error
    error("Internal error: fread read %zu values (expected %zu).\n", found_values, count);
  }
  return found_values;
}

size_t unsafe_fread(void *buffer, size_t size, size_t count, file_t *stream){
  // this function doesn't ensure we read as many values as we expected
  // it also doesn't retry on failure
  return fread_dispatch(buffer, size, count, stream);
}

size_t safe_fwrite(void *buffer, size_t size, size_t count, file_t *stream){
  size_t written_values = fwrite_dispatch(buffer, size, count, stream);
  if(written_values != count){
    // same scenarios as in safe_fread
    if(safe_feof(stream))
      error("%s", "Internal error: fread tried to read past end of file.\n");

    for(int i=0; i<MAX_WRITE_RETRIES; i++){
      // if we read a partial value, reset the counter back
      if(written_values) fseek_dispatch(stream, -1*((int)written_values), SEEK_CUR);

      // try to read again
      written_values = fwrite_dispatch(buffer, size, count, stream);
      if(written_values == count) return count;
    }

    // otherwise throw an error
    error("Internal error: fwrite wrote %zu values (expected %zu).\n", written_values, count);
  }
  return written_values;
}


/*
 * Below functions are simple wrappers
 * not necessary but included for consistency with safe_fread, safe_fwrite
 * can be extended later
 */

char safe_getc(file_t *f){
  if(!f) error("Error calling getc on null pointer.");
  return getc_dispatch(f);
}

int safe_fseek(file_t *f, long offset, int origin){
  int status = fseek_dispatch(f, offset, origin);
  if(status != 0) // basic error reporting on failure
    error("Internal error calling fseek on file with format %d (offset %ld, origin %d, code %d).",
          f->file_type, offset, origin, status);
  return status;
}

void safe_rewind(file_t *f){
  int status = fseek_dispatch(f, 0, SEEK_SET);
  if(status != 0) // basic error reporting on failure
    error("Internal error calling rewind on file with format %d (code %d).",
          f->file_type, status);
  return;
}

long safe_ftell(file_t *f){
  long pos = ftell_dispatch(f);
  if(pos == -1) // basic error reporting on failure
    error("Internal error calling ftell on file with format %d.", f->file_type);
  return pos;
}

int safe_feof(file_t *f){
  return feof_dispatch(f);
}
