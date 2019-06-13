/* ------- file: -------------------------- getline.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Dec  9 11:43:18 1999 --

       --------------------------                      ----------RH-- */

/* --- Routines for reading formatted input files. --  -------------- */


#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rh.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */
char *sgets(char *str, int num, char **input);
char *trimwhitespace(char *str);

/* --- Global variables --                             -------------- */

extern char messageStr[];


/* ------- begin -------------------------- getLine.c --------------- */

int getLine(FILE *inputFile, char *commentChar, char *line,
	    bool_t exit_on_EOF)
{
  const char routineName[] = "getLine";

  char *linePtr, dummy[2];

  /* --- Reads (into char array line) inputFile till first line
         that does not start with commentChar and returns length of
         line, or EOF at end of file. Empty lines are also ignored.
         --                                            -------------- */

  while((linePtr = fgets(line, MAX_LINE_SIZE, inputFile)) != NULL) {
    if ((sscanf(line, "%1s", dummy) > 0) && (*line != *commentChar))
      break;
  }
  if (linePtr == NULL) {
    if (exit_on_EOF) {
      sprintf(messageStr,
	      "Reached end of input file before all data was read");
      Error(ERROR_LEVEL_2, routineName, messageStr);
    } else
      return EOF;
  } else
    return strlen(line);
  return 0;
}
/* ------- end ---------------------------- getLine.c --------------- */

int getLineString(char **inputString, char *commentChar, char *line,
	    bool_t exit_on_EOF)
{
  const char routineName[] = "getLineString";

  char *linePtr, dummy[2];
  /* --- Does the same as getLine, but it reads from a string instead
         of a file.
         --                                            -------------- */
  while((linePtr = sgets(line, MAX_LINE_SIZE, inputString)) != NULL) {
    if ((sscanf(line, "%1s", dummy) > 0) &&
        (*(trimwhitespace(line)) != *commentChar)) break;
  }
  if (linePtr == NULL) {
    if (exit_on_EOF) {
      sprintf(messageStr, "Reached end of file before all data was read");
      Error(ERROR_LEVEL_2, routineName, messageStr);
    } else
      return EOF;
  } else
    return strlen(line);
  return 0;
}

/* ------- begin -------------------------- checkNread.c ------------ */

void checkNread(int Nread, int Nrequired, const char *routineName,
		int checkPoint)
{
  /* --- Check whether Nread igreater or equal to Nrequired and issue
         error otherwise --                           --------------- */

  if (Nread < Nrequired) {
    sprintf(messageStr, "Unable to read input file\n"
	    " At checkpoint %d. Needed at least %d, read %d item%s",
	    checkPoint, Nrequired, Nread, (Nread > 1) ? "s" : "");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- checkNread.c ------------ */

/* ------- end ---------------------------- UpperCase.c ------------- */

void UpperCase(char *string)
{
  register int n;

  /* --- Change string content to upper case --        -------------- */

  for (n = 0;  n < (int) strlen(string);  n++)
     string[n] = toupper(string[n]);
}
/* ------- end ---------------------------- UpperCase.c ------------- */

/* ------- begin -------------------------- substring.c ------------- */

char *substring(const char *string, int N0, int Nchar)
{
  static char destination[MAX_LINE_SIZE];
  int length = strlen(string);

  /* --- Extract a substring of length Nchar from source string,
         starting at position N0. --                   -------------- */

  if (N0 >= length) N0 = length;
  if (N0 + Nchar >= length) Nchar = length - N0;

  memcpy(destination, string + N0, Nchar);
  destination[Nchar] = '\0';

  return destination;
}
/* ------- end ---------------------------- substring.c ------------- */

char *readWholeFile(const char *filename)
{ /* Reads file into string */
  const char routineName[] = "readWholeFile";
  char *buffer = 0;
  unsigned long length;
  FILE *fp;

  if ((fp = fopen(filename, "rb")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", filename);
    Error(ERROR_LEVEL_2, routineName, filename);
  }
  if (fseek(fp, 0, SEEK_END) != 0) {
      sprintf(messageStr, "Error reading file %s", filename);
      Error(ERROR_LEVEL_2, routineName, filename);
  }
  length = ftell(fp);
  if (length == -1) {
      sprintf(messageStr, "Error reading file %s", filename);
      Error(ERROR_LEVEL_2, routineName, filename);
  }
  if (fseek(fp, 0, SEEK_SET) != 0) {
      sprintf(messageStr, "Error reading file %s", filename);
      Error(ERROR_LEVEL_2, routineName, filename);
  }
  buffer = (char*)malloc(length+1);
  if (buffer) {
      if (fread(buffer, 1, length, fp) < 1) {
          sprintf(messageStr, "Error reading file %s", filename);
          Error(ERROR_LEVEL_2, routineName, filename);
      }
  }
  if (fclose(fp) != 0) {
      sprintf(messageStr, "Error closing file %s", filename);
      Error(ERROR_LEVEL_2, routineName, filename);
  }
  buffer[length] = '\0';
  return buffer;
}


char *sgets(char *str, int num, char **input)
{ /* Works like fgets, but for a string instead of a file object */
    char *next = *input;
    int  numread = 0;

    while (numread + 1 < num && *next) {
        int isnewline = ( *next == '\n' );
        *str++ = *next++;
        numread++;
        /* newline terminates the line but is included */
        if (isnewline) break;
    }
    if (numread == 0) return NULL;  /* "eof" */
    /* must have hit the null terminator or end of line */
    *str = '\0';  /* null terminate this string */
    /* set up input for next call */
    *input = next;
    return str;
}


char *trimwhitespace(char *str)
{ /* Gets rid of leading white space from string */
  char *end;

  while(isspace((unsigned char) *str)) str++;
  if(*str == 0) return str;  /* Got all spaces? */
  end = str + strlen(str) - 1;
  *(end + 1) = 0;    /* Write new null terminator */
  return str;
}
