%option noyywrap
%option yylineno
%{
  #include <assert.h>
  #include <sys/stat.h>
  #include <sys/types.h>

  static FILE * fdepend = NULL, * myout = NULL;
  static char * fname;
  
  static char * paths[100] = { LIBDIR }, grid[80] = "quadtree";
  static int npath = 1, hasgrid = 0;
  static int incode;   // are we in a code block?
  static int somecode; // any code blocks in this file?

  static char * _stack[100]; int stack = -1;
  #define push(s) { char * f = malloc (strlen (s) + 1);	\
                    strcpy (f, s); _stack[++stack] = f; }
  #define pop()  _stack[stack--];

  static void singleslash (char * path, FILE * fp)
  {
    char * s = path, slash = 0;
    while (*s != '\0') {
      if (*s == '/') {
	if (!slash)
	  fputc (*s, fp);
	slash = 1;
      }
      else {
	slash = 0;
	fputc (*s, fp);
      }
      s++;
    }
  }

  static FILE * openpath (const char * name, const char * mode, char ** path) {
    int i;
    for (i = 0; i <= npath; i++) {
      char * p = malloc (strlen (paths[i]) + strlen (name) + 3);
      strcpy (p, paths[i]); strcat (p, "//"); strcat (p, name);
      FILE * fp = fopen (p, mode);
      if (fp) {
	if (fdepend) {
	  fputc ('\t', fdepend); singleslash (p, fdepend);
	  fputs (" \\\n", fdepend);
	}
	*path = p;
	return fp;
      }
      else
	free (p);
    }
    return NULL;
  }

  static int yyerror(const char * s);
  static int comment(void);
  static void echo() {
    if (myout) {
      if (incode) {
	fputs (yytext, myout);
	somecode = 1;
      }
      else { // only keep newlines
	char * s = yytext;
	while (*s != '\0') {
	  if (*s == '\n')
	    fputc ('\n', myout);
	  s++;
	}
      }
    }
  }
%}

SP   [ \t]
WS   [ \t\v\n\f]
ES   (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
BEGINCODE ^[SP]*[~]{3,}c[^\n]*\n
ENDCODE   ^[SP]*[~]{3,}[^\n]*\n

%%

{BEGINCODE} {
  if (incode) {
    yylineno--;
    return yyerror ("code blocks cannot be nested");
  }
  incode = 1;
  if (myout) fputc ('\n', myout);
}

{ENDCODE} {
  if (!incode) {
    yylineno--;
    return yyerror ("not in a code block");
  }
  incode = 0;
  if (myout) fputc ('\n', myout);
}

^{SP}*#{SP}*include{SP}+\".*\"*{SP}*\n {
  echo();
  char * s = strchr(yytext, '"');
  s++;
  char * e = &s[strlen(s) - 1];
  while (*e == ' ' || *e == '\t' || *e == '"' || *e == '\n') {
    *e = '\0'; e--;
  }
  char * path;
  FILE * fp = openpath (s, "r", &path);
  if (fp != NULL) {
    push (path);
    free (path);
    fclose (fp);
  }
}

^{SP}*#{SP}*define{SP}+GRIDNAME{WS}+ {
    echo();
    hasgrid = 1;
    char * s = fname;
    while (strchr (s, '/')) s = strchr (s, '/') + 1;
    strcpy (grid, s);
    if ((s = strchr (grid, '.'))) *s = '\0';
}

"/*"              { echo(); if (comment()) return 1; }
"//".*            { echo(); /* consume //-comment */ }
.                   echo();
[\n]                echo();
({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+  echo(); /* STRING_LITERAL */

%%

int yyerror (const char * s)
{
  if (0) yyunput (0, NULL); // just prevents 'yyunput unused' compiler warning
  fprintf (stderr, "%s:%d: error: %s\n", fname, yylineno, s);
  return 1;
}

static int getput(void)
{
  int c = input();
  if (myout)
    fputc (c, myout);
  return c;
}

static int comment(void)
{
  int c;
  while ((c = getput()) != 0) {
    if (c == '*') {
      while ((c = getput()) == '*')
	;
      if (c == '/')
	return 0;
      if (c == 0)
	break;
    }
  }
  return yyerror ("unterminated comment");
}

void stripname (char * path)
{
  char * s = &path[strlen(path)];
  while (s != path && *s != '/')
    *s-- = '\0';
  if (s == path)
    strcpy (path, ".");
  else
    *s = '\0';
}

char * stripslash (char * path)
{
  char * strip = malloc (strlen (path) + 1), * s = path, * o = strip;
  int slash = 0;
  do {
    if (*s == '/') {
      if (!slash)
	*o++ = *s;
      slash = 1;
    }
    else {
      *o++ = *s;
      slash = 0;
    }
  } while (*s++ != '\0');
  return strip;
}

static int include (char * file, FILE * fin, FILE * fout)
{
  fname = stripslash (file);
  paths[npath] = malloc (strlen (file) + 1);
  strcpy (paths[npath], file);
  stripname (paths[npath]);
  yyin = fin;
  myout = fout;
  yylineno = 1;
  incode = somecode = 0;
  long header = fout ? ftell (fout) : 0;
  int ret = yylex();
  if (fout && !somecode) { 
    // no literate code block found, assume the entire file is pure code
    fseek (fout, header, SEEK_SET);
    rewind (fin);
    int c;
    while ((c = fgetc (fin)) != EOF)
      fputc (c, fout);
  }
  free (fname);
  return ret;
}

FILE * writepath (char * path, const char * mode);
void cleanup (int status, const char * dir);

static int compdir (char * file, char ** out, int nout, const char * dir)
{
  push (file);
  while (stack >= 0) {
    char * path = pop();
    FILE * fin = fopen (path, "r");
    if (fin == NULL) {
      perror (path);
      cleanup (1, dir);
    }
    FILE * fout = NULL;
    if (dir) {
      char * file = strstr (path, "//");
      if (file) file += 2; else file = path;
      char * out = malloc (strlen (dir) + strlen (file) + 2);
      strcpy (out, dir);
      strcat (out, "/");
      strcat (out, file);
      fout = writepath (out, "w");
      if (fout == NULL) {
	perror (out);
	cleanup (1, dir);
      }
      free (out);

      fputs ("# 1 \"", fout);
      singleslash (path, fout);
      fputs ("\"\n", fout);
    }
    if (include (path, fin, fout))
      cleanup (1, dir);
    fclose (fin);
    if (fout)
      fclose (fout);
    out[nout++] = path;
  }
  return nout;
}

int includes (int argc, char ** argv, char ** out, 
	      char ** grid1, int * default_grid,
	      const char * dir)
{
  int depend = 0, nout = 0;
  char * file = NULL, * output = NULL;
  int i;
  for (i = 1; i < argc; i++) {
    if (!strncmp (argv[i], "-grid=", 6))
      strcpy (grid, &argv[i][6]);
    else if (!strcmp (argv[i], "-MD"))
      depend = 1;
    else if (!strcmp (argv[i], "-o"))
      output = argv[i + 1];
    else if (argv[i][0] != '-' && \
	     !strcmp (&argv[i][strlen(argv[i]) - 2], ".c")) {
      if (file) {
	fprintf (stderr, "usage: include [OPTIONS] FILE.c\n");
	cleanup (1, dir);
      }
      file = argv[i];
    }
  }
  if (depend && file) {
    if (!output) output = file;
    char ndep[80], * s = &output[strlen(output)-1];
    while (*s != '.' && s != output) s--;
    if (output != file || s == output)
      /* generate dep files with suffixes included for -o option */
      strcpy (ndep, output);
    else {
      *s = '\0';
      strcpy (ndep, output);
      *s = '.';
    }
    if (strlen(ndep) < 2 || strcmp (&ndep[strlen(ndep)-2], ".d"))
      strcat (ndep, ".d");
    else
      output[strlen(ndep)-2] = '\0'; // strip trailing ".d";
    fdepend = fopen (ndep, "w");
    if (!fdepend) {
      perror (ndep);
      cleanup (1, dir);
    }
    fprintf (fdepend, "%s:\t\\\n", output);
  }
  if (file) {
    nout = compdir (file, out, 0, dir);
    if (!hasgrid) {
      char * path, gridpath[80] = "grid/";
      strcat (gridpath, grid); strcat (gridpath, ".h");
      FILE * fp = openpath (gridpath, "r", &path);
      if (!fp) {
	fprintf (stderr, "include: invalid grid '%s': ", grid);
	perror ("");
	cleanup (1, dir);
      }
      fclose (fp);
      nout = compdir (path, out, nout, dir);
      hasgrid = 0;
    }
    char * path;    
    FILE * fp = openpath ("common.h", "r", &path);
    assert (fp);
    out[nout++] = path;
    fclose (fp);
  }
  if (fdepend) {
    fputc ('\n', fdepend);
    fclose (fdepend);
  }
  *grid1 = grid;
  *default_grid = !hasgrid;
  return nout;
}