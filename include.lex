%option noyywrap
%{
  #include <assert.h>

  static FILE * fdepend = NULL;
  static char * fname;
  static int line;
  
  static char * paths[100] = { LIBDIR }, grid[80] = "quadtree";
  static int npath = 1, hasgrid = 0;

  static char * _stack[100]; int stack = -1;
  #define push(s) { char * f = malloc (sizeof (char)*(strlen (s) + 1));	\
                    strcpy (f, s); _stack[++stack] = f; }
  #define pop()  _stack[stack--];

  static FILE * openpath (const char * name, const char * mode, char ** path) {
    int i;
    for (i = 0; i <= npath; i++) {
      char * p = malloc (sizeof (char)*(strlen (paths[i]) + strlen (name) + 3));
      strcpy (p, paths[i]); strcat (p, "//"); strcat (p, name);
      FILE * fp = fopen (p, mode);
      if (fp) {
	if (fdepend)
	  fprintf (fdepend, "\t%s \\\n", p);
	*path = p;
	return fp;
      }
      else
	free (p);
    }
    return NULL;
  }

#define YY_INPUT(buf,result,max_size) {				\
    int c = fgetc(yyin);					\
    result = (c == EOF) ? YY_NULL : (buf[0] = c, 1);		\
    if (c == '\n') line++;					\
  }

  static int yyerror(const char * s);
  static int comment(void);
%}

SP  [ \t]

%%

^{SP}*#{SP}*include{SP}+\".*\"*{SP}*\n {
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

^{SP}*#{SP}*define{SP}+foreach[ \t(]    { hasgrid = 1; }

"/*"                                    { if (comment()) return 1; }
"//".*                                  { /* consume //-comment */ }
.                                       ;
[\n]                                    ;

%%

int yyerror (const char * s)
{
  fprintf (stderr, "%s:%d: error: %s\n", fname, line, s);
  return 1;
}

static int comment(void)
{
  int c;
  while ((c = input()) != 0) {
    if (c == '*') {
      while ((c = input()) == '*')
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
  char * strip = malloc (sizeof (char)*(strlen (path) + 1)), * s = path, * o = strip;
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

static int include (char * file, FILE * fin)
{
  fname = stripslash (file);
  paths[npath] = malloc (sizeof (char)*(strlen (file) + 1));
  strcpy (paths[npath], file);
  stripname (paths[npath]);
  yyin = fin;
  line = 1;
  int ret = yylex();
  free (fname);
  return ret;
}

static int compdir (char * file, char ** out)
{
  int nout = 0;
  push (file);
  while (stack >= 0) {
    char * path = pop();
    FILE * fin = fopen (path, "r");
    if (fin == NULL) {
      perror (path);
      exit (1);
    }
    if (include (path, fin))
      exit (1);
    fclose (fin);
    out[nout++] = path;
  }
  return nout;
}

int includes (int argc, char ** argv, char ** out, char ** grid1)
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
    else if (argv[i][0] != '-' && !strcmp (&argv[i][strlen(argv[i]) - 2], ".c")) {
      if (file) {
	fprintf (stderr, "usage: include [OPTIONS] FILE.c\n");
	exit (1);
      }
      file = argv[i];
    }
  }
  if (depend && output && file) {
    char ndep[80], * s = &output[strlen(output)-1];
    while (*s != '.' && s != output) s--;
    if (1/*s == output*/) /* always generate dep files with suffixes included */
      strcpy (ndep, output);
    else {
      *s = '\0';
      strcpy (ndep, output);
      *s = '.';
    }
    strcat (ndep, ".d");
    fdepend = fopen (ndep, "w");
    if (!fdepend) {
      perror (ndep);
      exit (1);
    }
    fprintf (fdepend, "%s:\t\\\n", output);
  }
  if (file) {
    nout = compdir (file, out);
    if (!hasgrid) {
      char * path, gridpath[80] = "grid/";
      strcat (gridpath, grid); strcat (gridpath, ".h");
      FILE * fp = openpath (gridpath, "r", &path);
      if (!fp) {
	fprintf (stderr, "include: invalid grid '%s': ", grid);
	perror ("");
	exit (1);
      }
      out[nout++] = path;
      fclose (fp);
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
  *grid1 = hasgrid ? NULL : grid;
  return nout;
}