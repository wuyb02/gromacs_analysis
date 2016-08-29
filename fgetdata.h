/* fgetdata.h - header for the fget* function.
 * fgetline() and fgetword() were written by Richard Heathfield.
 * The version documented here is not the final version.
 *
 * fgetline reads a line from fp into a dynamically
 * allocated string pointed to by *line. Note: it *is*
 * acceptable to do char *p = NULL and pass &p; but it is
 * *not* acceptable to do this: char s[10]; and pass &s!!
 *
 * fgetword does the same thing, only for words rather than
 * entire lines.
 *
 * If you pass the address of a NULL pointer, you needn't
 * set size before passing &size. Otherwise, it must be
 * set to the number of bytes pointed to by the pointer
 * whose address you pass in.
 *
 * The string will be dynamically resized if need be, so as
 * to be able to store the whole line (or word).
 *
 * It is the caller's responsibility to release any memory
 * allocated by these functions.
 *
 * flags:
 *   FGDATA_REDUCE : after the line is captured, reduce the
 *                   amount of memory to the minimum necessary
 *                   to hold the string.
 *
 * The function returns 0 on success, 1 on EOF, or a negative
 * number if something else went wrong (almost certainly a
 * lack of memory, but could be a stream error).
 *
 * Note: the source uses assertions for various kinds of
 * program error. You should take advantage of these when
 * testing your code.
 */

#ifndef FGETDATA_H_
#define FGETDATA_H_ 1
#include <stdio.h>

#define FGDATA_BUFSIZ 255 /* adjust to taste */
#define FGDATA_WRDSIZ sizeof("floccinaucinihilipilification")
#define FGDATA_REDUCE  1

int fgetline(char **line, size_t *size, size_t maxrecsize, FILE *fp, unsigned int flags);
int fgetword(char **word, size_t *size, const char *delimiters, size_t maxrecsize, FILE *fp, unsigned int flags);

#endif
