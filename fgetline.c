/* fgetline reads a line from fp into a dynamically
 * allocated string pointed to by *line. Note: it *is*
 * acceptable to do char *p = NULL and pass &p; but it is
 * *not* acceptable to do this: char s[10]; and pass &s!!
 *
 * If you pass the address of a NULL pointer, you needn't
 * set size before passing &size. Otherwise, it must be
 * set to the number of bytes pointed to by the pointer
 * whose address you pass in.
 *
 * The string will be dynamically resized if need be, so as
 * to be able to store the whole line.
 *
 * It is the caller's responsibility to release any memory
 * allocated by this function.
 *
 * flags:
 *   FGL_REDUCE : after the line is captured, reduce the
 *                amount of memory to the minimum necessary
 *                to hold the string. This feature might
 *                fail silently (pretty darn unlikely, but still...)
 *
 * The function returns 0 on success, 1 on EOF, or a negative
 * number if something else went wrong (almost certainly a
 * lack of memory, but could be a stream error).
 *
 * Assertions:
 * If NDEBUG is not defined, the following program bugs will
 * fire an assertion failure: line is NULL, size is NULL, fp is NULL.
 */

#include <assert.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fgetdata.h"

int fgetline(char **line, size_t *size, size_t maxrecsize, FILE *fp, unsigned int flags)
{
  int Result = 0;
  size_t count = 0;
  size_t newsize = 0;
  int ch = 0;
  char *tmp = NULL;

  assert(fp   != NULL);
  assert(line != NULL);
  assert(size != NULL);

  if(NULL == *line)
  {
    *line = malloc(FGDATA_BUFSIZ);
    if(*line != NULL)
    {
      *size = FGDATA_BUFSIZ;
    }
    else
    {
      Result = -1;
    }
  }
  if(0 == Result && feof(fp))
  {
    Result = 1;
  }

  //while(0 == Result && (ch = fgetc(fp)) != '\n' && ch != 13 && ch != EOF)
  while(0 == Result && (ch = fgetc(fp)) != '\n' && ch != EOF)
  {
    if(count + 2 >= *size)
    {
      /* This realloc strategy has been revised recently, and is subject to further revision */
      newsize = 3 * (count + *size) / 2;
      if(newsize > maxrecsize)
      {
        Result = -2;
      }
      else
      {
        tmp = realloc(*line, newsize);
        if(NULL == tmp)
        {
          Result = -3;
        }
        else
        {
          *line = tmp;
          tmp = NULL;
          *size = newsize;
        }
      }
    }
    if(0 == Result && ch!='\r')
    {
      (*line)[count++] = ch;
    }
  }
  if(0 == Result)
  {
    (*line)[count] = '\0';
  }
  if(ch == EOF)
  {
    Result = 1;
  }

  if(0 == Result)
  {
    if(flags & FGDATA_REDUCE)
    {
      newsize = strlen(*line) + 1;
      tmp = realloc(*line, newsize);
      if(tmp != NULL)
      {
        *line = tmp;
        //*size = newsize-1;
        *size = newsize;
      }
    }
  }

  if(0 == Result && ferror(fp))
  {
    Result = -4;
  }

  return Result;
}
