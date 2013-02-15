%option noyywrap
%{
  #include <unistd.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <sys/wait.h>
  #include <assert.h>

  FILE * fdepend = NULL;

  int line, scope, para, inforeach, foreachscope, foreachpara, hasgrid = 0;
  char foreachs[80], * fname, grid[80] = "quadtree";
  
  char * _stack[100]; int stack = -1;
  #define push(s) { char * f = malloc (sizeof (char)*(strlen (s) + 1));	\
                    strcpy (f, s); _stack[++stack] = f; }
  #define pop()  _stack[stack--];
  
  struct { char * v; int b; } _varstack[100]; int varstack = -1;
  void varpush (const char * s, int b) {
    if (s[0] != '\0') {
      char * f = malloc (sizeof (char)*(strlen (s) + 1));
      strcpy (f, s); 
      _varstack[++varstack].v = f; 
      _varstack[varstack].b = b;
    }
  }
  void varpop () {
    while (varstack >= 0 && _varstack[varstack].b > scope)
      free (_varstack[varstack--].v);
  }

  struct { char * v, * fname; int line; } _newvar[100]; int newvar = 0;
  int newvarpush (const char * s) {
    int i = 0;
    for (i = 0; i < newvar; i++)
      if (!strcmp (s, _newvar[i].v)) {
	fprintf (stderr, "%s:%d: error: redefinition of '%s'\n", fname, line, s);
	fprintf (stderr, "%s:%d: note: previous definition of '%s' was here\n", 
		 _newvar[i].fname, _newvar[i].line, s);
	return 1;
      }
    char * f = malloc (sizeof (char)*(strlen (s) + 1));
    strcpy (f, s);
    _newvar[newvar].v = f;
    f = malloc (sizeof (char)*(strlen (fname) + 1));
    strcpy (f, fname);
    _newvar[newvar].fname = f;
    _newvar[newvar].line = line;
    newvar++;
    varpush (s, 0);
    return 0;
  }
  
  char * paths[100] = { LIBDIR };
  int npath = 1;
  
  FILE * openpath (const char * name, const char * mode, char ** path) {
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

  void endforeach() {
    if (varstack >= 0) {
      fputc ('\n', yyout);
      int i = varstack;
      while (i >= 0) {
	char * v = _varstack[i--].v;
	fprintf (yyout, "#undef %s\n", v);
      }
    }
    fprintf (yyout, " end_%s();\n#line %d\n", foreachs, line - 1);
  }

#define YY_INPUT(buf,result,max_size) {				\
    int c = fgetc(yyin);					\
    result = (c == EOF) ? YY_NULL : (buf[0] = c, 1);		\
    if (c == '\n') line++;					\
  }

  int yyerror(const char * s);
  int getput(void); 
  int comment(void);
%}

ID  [a-zA-Z_0-9]
SP  [ \t]
ES  (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
WS  [ \t\v\n\f]

%%
"("                ECHO; para++;

")" {
  ECHO; para--;
  if (para < 0)
    return yyerror ("mismatched ')'");
}

"{" {
  ECHO; scope++;
}

"}" {
  ECHO; scope--;
  if (scope < 0)
    return yyerror ("mismatched '}'");
  varpop();
  if (inforeach && scope == foreachscope) {
    inforeach = 0;
    endforeach ();
  }
}

foreach{ID}* {
  ECHO;
  strcpy (foreachs, yytext);
  inforeach = 1; foreachscope = scope; foreachpara = para;
  int c = getput();
  while (strchr (" \t", c)) c = getput();
  if (c != '(')
    return yyerror ("expecting '('");
  para++;
  while (para > foreachpara && c != EOF) {
    c = getput();
    if (c == '(') para++;
    if (c == ')') para--;
  }
  if (c != ')')
    return yyerror ("expecting ')");
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
  if (inforeach && scope == foreachscope && para == foreachpara) {
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
}

[^{ID}]new{WS}+var{WS}+({ID}|,|{WS})+ {
  /* new var ID, ID, ID ... */
  fputc (yytext[0], yyout);
  if (yytext[0] == '(') para++;
  char * s = strstr(yytext,"var"); s += 4;
  while (strchr(" \t\v\n\f", *s)) s++;
  while (*s != '\0') {
    char * var = s;
    while (!strchr(" \t\v\n\f,\0", *s)) s++;
    if (*s != '\0') *s++ = '\0';
    while (*s != '\0' && strchr(" \t\v\n\f,", *s))
      s++;
    if (newvarpush (var))
      return 1;
  }
  fprintf (yyout, "\n#line %d\n", line);
}

[^{ID}]var{SP}+{ID}+ {
  /* var ID*/
  ECHO;
  if (yytext[0] == '(') para++;
  if (para == 1) { /* ignore nested functions */
    char * var = &yytext[5];
    while (strchr (" \t", *var)) var++;
    varpush (var, scope + 1);
  }
}

"#" {
  ECHO;                     
  register int oldc = 0, c;
  for (;;) {
    while ((c = getput()) != '\n' && c != EOF)
      oldc = c;    /* eat up text of preproc */
    if (c == EOF || oldc != '\\')
      break;
  }
  fprintf (yyout, "\n#line %d\n", line);
}

"/*"                                    { ECHO; if (comment()) return 1; }
"//"[^\n]*                              { ECHO; /* consume //-comment */ }
({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+	{ ECHO; /* STRING_LITERAL */ }

%%

int yyerror (const char * s)
{
  fprintf (stderr, "%s:%d: error: %s\n", fname, line, s);
  return 1;
}

int comment(void)
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

int getput(void)
{
  int c = input();
  fputc (c, yyout);
  return c;
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
  line = 1, scope = 0, para = 0;
  inforeach = 0, foreachscope = 0, foreachpara = 0;
  int ret = yylex();
  if (!ret) {
    if (scope > 0)
      ret = yyerror ("mismatched '{'");
    else if (para > 0)
      ret = yyerror ("mismatched '('");
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
  char command[80] = "rm -r -f ";
  strcat (command, dir);
  int s = system (command);
  exit (status);
}

void compdir (char * file)
{
  int first = 1, hadgrid = 1;
  char * path;
  FILE * fp = openpath ("common.h", "r", &path);
  assert (fp);
  push (path);
  free (path);
  fclose (fp);
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
  
  /* new variables */
  char * out = malloc (sizeof (char)*(strlen (dir) + strlen ("grid.h") + 2));
  strcpy (out, dir); strcat (out, "/grid.h");
  FILE * fout = fopen (out, "w");
  fputs ("#include \"common.h\"\n", fout);
  if (newvar == 0)
    fprintf (fout, "struct _Data {};\n");
  else {
    fprintf (fout, "struct _Data {\n  double");
    int i;
    for (i = 0; i < newvar - 1; i++)
      fprintf (fout, " %s,", _newvar[i].v);
    fprintf (fout, " %s;", _newvar[i].v);
    fputs ("\n};\nvar", fout);
    for (i = 0; i < newvar - 1; i++)
      fprintf (fout, "\t%s = offsetof(Data,%s),\n", _newvar[i].v, _newvar[i].v);
    fprintf (fout, "\t%s = offsetof(Data,%s);\n", _newvar[i].v, _newvar[i].v);
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
  int depend = 0;
  char * file = NULL, * output = NULL;
  int i;
  for (i = 1; i < argc; i++) {
    if (!strncmp (argv[i], "-grid=", 6))
      strcpy (grid, &argv[i][6]);
    else if (!strcmp (argv[i], "-MD"))
      depend = 1;
    else if (!strcmp (argv[i], "-o")) {
      output = argv[i + 1];
      strcat (command, " ");
      strcat (command, argv[i]);
    }
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
      cleanup (1);
    }
    fprintf (fdepend, "%s:\t\\\n", output);
  }
  if (file) {
    compdir (file);
    strcat (command, " ");
    strcat (command, dir);
    strcat (command, "/");
    strcat (command, file);
  }
  if (fdepend) {
    fputc ('\n', fdepend);
    fclose (fdepend);
  }
  /* compilation */
  int status = system (command);
  if (status == -1 ||
      (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || WTERMSIG (status) == SIGQUIT)))
    cleanup (1);
  cleanup (WEXITSTATUS (status));
  return 0;
}
