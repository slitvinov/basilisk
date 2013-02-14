%{
  #include <unistd.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <sys/wait.h>
  #include <assert.h>

  int line, brack, para, inforeach, foreachbrack, foreachpara, hasgrid = 0;
  char foreachs[80], * fname, grid[80] = "quadtree";
  
  char * _stack[100]; int stack = -1;
  #define push(s) { char * f = malloc (sizeof (char)*(strlen (s) + 1));	\
                    strcpy (f, s); _stack[++stack] = f; }
  #define pop()  _stack[stack--];
  
  struct { char * v; int b; } _varstack[100]; int varstack = -1;
  void varpush (char * s, int b) {
    if (s[0] != '\0') {
      char * f = malloc (sizeof (char)*(strlen (s) + 1));
      strcpy (f, s); 
      _varstack[++varstack].v = f; 
      _varstack[varstack].b = b;
    }
  }
  void varpop () {
    while (varstack >= 0 && _varstack[varstack].b > brack)
      free (_varstack[varstack--].v);
  }
  
  char * paths[100] = { "/home/popinet/local/src/atmosphere" };
  int npath = 1;
  
  FILE * openpath (const char * name, const char * mode, char ** path) {
    int i;
    for (i = 0; i <= npath; i++) {
      char * p = malloc (sizeof (char)*(strlen (paths[i]) + strlen (name) + 3));
      strcpy (p, paths[i]); strcat (p, "//"); strcat (p, name);
      FILE * fp = fopen (p, mode);
      if (fp) {
	*path = p;
	return fp;
      }
      else
	free (p);
    }
    return NULL;
  }

  void endforeach() {
    if (varstack >= 0) {
      fputc ('\n', yyout);
      int i = varstack;
      while (i >= 0) {
	char * v = _varstack[i--].v;
	fprintf (yyout, "#undef %s\n", v);
      }
    }
    fprintf (yyout, " end_%s();\n#line %d\n", foreachs, line);
  }

  int _getput (int c) {
    fputc (c, yyout);
    if (c == '\n') line++;
    return c;
  }
#define getput() _getput(input())
%}

ID [a-zA-Z_]
SP [ \t]

%%
\n                 ECHO; ++line;

"("                ECHO; para++;

")"                ECHO; para--;

"{" {
  ECHO; brack++;
}

"}" {
  ECHO; brack--;
  if (brack < 0) {
    fprintf (stderr, 
	     "%s:%d: error: mismatched bracket\n", fname, line);
    return 1;
  }
  varpop();
  if (inforeach && brack == foreachbrack) {
    inforeach = 0;
    endforeach ();
  }
}

foreach{ID}* {
  ECHO;
  strcpy (foreachs, yytext);
  inforeach = 1; foreachbrack = brack; foreachpara = para;
  int c = getput();
  while (strchr (" \t", c)) c = getput();
  if (c != '(') {
    fprintf (stderr, "%s:%d: error: expecting '('\n", fname, line);
    return 1;
  }
  para++;
  while (para > foreachpara && c != EOF) {
    c = getput();
    if (c == '(') para++;
    if (c == ')') para--;
  }
  if (c != ')') {
    fprintf (stderr, "%s:%d: error: expecting ')'\n", fname, line);
    return 1;
  }
  if (varstack >= 0) {
    fputc ('\n', yyout);
    int i = varstack;
    while (i >= 0) {
      char * v = _varstack[i--].v;
      fprintf (yyout, "#define %s(i,j) val(%s,i,j)\n", v, v);
    }
    fprintf (yyout, "#line %d\n", line);
  }
}

end_foreach{ID}*{SP}*"()" {
  if (strncmp(&yytext[4], foreachs, strlen (foreachs))) {
    fprintf (stderr, 
	     "%s:%d: error: "
	     "%s() loop ended with %s\n", 
	     fname, line, foreachs, yytext);
    return 1;
  }
}

;  {
  ECHO;
  if (inforeach && brack == foreachbrack && para == foreachpara) {
    endforeach();
    inforeach = 0;
  }
}

^{SP}*#{SP}*include{SP}+\".*\"*{SP}*\n  {
  ECHO;
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
  line++;
}

[^a-z0-9]var[ \t]+ {
  ECHO;
  if (yytext[0] == '(') para++;
  if (para <= 1) { /* no nested functions */
    register int c = getput();
    for (;;) {
      char var[80] = "", * v = var;
      while (c != EOF && 
	     ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9'))) {
	*v++ = c;
	c = getput();
      }
      *v = '\0';
      
      if (para == 1) { /* function prototype */
	varpush (var, brack + 1);
	if (c == ')') para--;
	break;
      }

      /* declaration in statement */
      varpush (var, brack);
      while (c != ',' && c != ';' && c != EOF)
	c = getput();
      if (c == EOF) {
	fprintf (stderr, "%s:%d: error: end-of-file in declaration\n", fname, line);
	return 1;	
      }
      if (c == ';')
	break;
      c = getput();
      while (c == ' ' || c == '\t')
	c = getput();
    }
  }
}

"/*" {
  ECHO;
  register int c;
  for (;;) {
    while ((c = input()) != '*' && c != EOF)
      fputc (c, yyout);    /* eat up text of comment */
    if (c == '*') {
      fputc (c, yyout);
      while ( (c = input()) == '*' )
	fputc (c, yyout);
      fputc (c, yyout);
      if ( c == '/' )
	break;    /* found the end */
    }
    if (c == EOF) {
      fprintf (stderr, 
	       "%s:%d: error: "
	       "EOF in comment\n", 
	       fname, line);
      return 1;
    }
  }
}

"//" {
  ECHO;
  register int c;
  while ((c = input()) != '\n' && c != EOF)
    fputc (c, yyout);    /* eat up text of comment */
  if (c == '\n') line++;
  fputc (c, yyout);
}

"#" {
  ECHO;                     
  register int oldc = 0, c;
  for (;;) {
    while ((c = input()) != '\n' && c != EOF) {
      fputc (c, yyout);
      oldc = c;    /* eat up text of preproc */
    }
    if (c == '\n') line++;
    fputc (c, yyout);
    if (c == EOF || oldc != '\\')
      break;
  }
  fprintf (yyout, "\n#line %d\n", line);
}

\"[^\"]*\" ECHO;

%%

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

char dir[] = ".endforXXXXXX";

int endfor (char * file, FILE * fin, FILE * fout)
{
  fname = stripslash (file);
  fprintf (fout, "# 1 \"%s\"\n", fname);
  paths[npath] = malloc (sizeof (char)*(strlen (file) + 1));
  strcpy (paths[npath], file);
  stripname (paths[npath]);
  yyin = fin;
  yyout = fout;
  line = 1, brack = 0, para = 0;
  inforeach = 0, foreachbrack = 0, foreachpara = 0;
  int ret = yylex();
  if (!ret && brack > 0) {
    fprintf (stderr, "%s:%d: error: unclosed bracket\n", fname, line);
    ret = 1;
  }
  free (fname);
  return ret;
}

FILE * writepath (char * path, const char * mode)
{
  char * s = path;
  while ((s = strchr (s, '/'))) {
    *s = '\0';
    if (access (path, R_OK|W_OK|X_OK)) {
      if (strlen(path) > 5 && !strcmp(&path[strlen(path)-5], "/grid"))
	hasgrid = 1;
      if (mkdir (path, 0700))
	return NULL;
    }
    *s++ = '/';
  }
  return fopen (path, mode);
}

void cleanup (int status)
{
  char command[1000] = "";
  strcat (command, "cat ");
  strcat (command, dir);
  strcat (command, "/grid.h; rm -r -f ");
  strcat (command, dir);
  system (command);
  exit (status);
}

void compdir (char * file)
{
  int first = 1, hadgrid = 1;
  push (file);
  for (;;) {
    while (stack >= 0) {
      char * path = pop();
      FILE * fin = fopen (path, "r");
      if (fin == NULL) {
	perror (path);
	cleanup (1);
      }
      char * file = strstr (path, "//"); 
      if (file) file += 2; else file = path;
      char * out = malloc (sizeof (char)*(strlen (dir) + strlen (file) + 2));
      strcpy (out, dir);
      strcat (out, "/");
      strcat (out, file);
      FILE * fout = writepath (out, "w");
      if (fout == NULL) {
	perror (out);
	cleanup (1);
      }
      if (first) {
	fprintf (fout, "#include \"grid.h\"\n");
	first = 0;
      }
      if (endfor (path, fin, fout))
	cleanup (1);
      fclose (fout);
      free (out);
      fclose (fin);
      free (path);
    }
    if (hasgrid)
      break;
    else {
      char * path, gridpath[80] = "grid/";
      strcat (gridpath, grid); strcat (gridpath, ".h");
      FILE * fp = openpath (gridpath, "r", &path);
      assert (fp);
      push (path);
      free (path);
      fclose (fp);
      hadgrid = 0;
    }
  }
  
  /* global variables */
  char * out = malloc (sizeof (char)*(strlen (dir) + strlen ("grid.h") + 2));
  strcpy (out, dir); strcat (out, "/grid.h");
  FILE * fout = fopen (out, "w");
  if (varstack < 0)
    fprintf (fout, "struct _Data {};\n");
  else {
    fprintf (fout, "struct _Data {\n  double");
    int i;
    for (i = 0; i < varstack; i++) {
      assert (_varstack[i].b == 0);
      fprintf (fout, " %s,", _varstack[i].v);
    }
    assert (_varstack[i].b == 0);
    fprintf (fout, " %s;", _varstack[i].v);
    fprintf (fout, "\n};\n");
  }
  if (!hadgrid)
    fprintf (fout, "#include \"grid/%s.h\"\n", grid);
  free (out);
  fclose (fout);
}

int main (int argc, char ** argv)
{
  char * cc = getenv ("CC"), command[1000];
  if (cc == NULL) cc = "cc";
  strcpy (command, cc);
  char * file = NULL;
  int i;
  for (i = 1; i < argc; i++) {
    if (!strncmp (argv[i], "-grid=", 6))
      strcpy (grid, &argv[i][6]);
    else if (argv[i][0] != '-' && !strcmp (&argv[i][strlen(argv[i]) - 2], ".c")) {
      if (file) {
	fprintf (stderr, "usage: endfor -grid=[GRID] [OPTIONS] FILE.c\n");
	return 1;
      }
      file = argv[i];
    }
    else {
      strcat (command, " ");
      strcat (command, argv[i]);
    }
  }
  if (!mkdtemp (dir)) {
    perror (dir);
    return 1;
  }
  if (file) {
    compdir (file);
    strcat (command, " ");
    strcat (command, dir);
    strcat (command, "/");
    strcat (command, file);
  }
  /* compilation */
  int status = system (command);
  if (status == -1 ||
      (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || WTERMSIG (status) == SIGQUIT)))
    cleanup (1);
  cleanup (WEXITSTATUS (status));
  return 0;
}
