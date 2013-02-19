%option noyywrap
%{
  #include <unistd.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <sys/wait.h>
  #include <assert.h>

  int debug = 0;

  int nvar = 0;
  int line, scope, para, inforeach, foreachscope, foreachpara, invardecl;
  char foreachs[80], * fname;
  
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
  
  void endforeach() {
    if (varstack >= 0) {
      fputc ('\n', yyout);
      int i = varstack;
      while (i >= 0) {
	char * v = _varstack[i--].v;
	fprintf (yyout, "#undef %s\n", v);
      }
      fprintf (yyout, " end_%s();\n#line %d\n", foreachs, line - 1);
    }
    else
      fprintf (yyout, " end_%s();", foreachs);
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
  invardecl = 0;
}

[^{ID}]{WS}*var{WS}+{ID}+ {
  ECHO;
  if (yytext[0] == '(') para++;
  char * var = &strstr(yytext,"var")[4];
  while (strchr (" \t\v\n\f", *var)) var++;
  if (para == 0) { /* declaration */
    if (debug)
      fprintf (stderr, "%s:%d: declaration: %s\n", fname, line, var);
    varpush (var, scope);
    invardecl = 1;
  }
  else if (para == 1) { /* function prototype (no nested functions) */
    if (debug)
      fprintf (stderr, "%s:%d: proto: %s\n", fname, line, var);
    varpush (var, scope + 1);
  }
}

,{WS}*{ID}+ {
  ECHO;
  if (invardecl) {
    char * var = &yytext[1];
    while (strchr (" \t\v\n\f", *var)) var++;
    if (debug)
      fprintf (stderr, "%s:%d: declaration: %s\n", fname, line, var);
    varpush (var, scope);
  }
}

[^{ID}]new{WS}+var[^{ID}] {
  /* new var */
  if (yytext[0] == '(') para++;
  fputc (yytext[0], yyout);
  fprintf (yyout, "%d", nvar++);
  unput (yytext[yyleng - 1]);
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
"//".*                                  { ECHO; /* consume //-comment */ }
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

void stripname (char * path);
char * stripslash (char * path);
int includes (int argc, char ** argv, char ** out, char ** grid);

char dir[] = ".qccXXXXXX";

int endfor (char * file, FILE * fin, FILE * fout)
{
  fname = stripslash (file);
  fprintf (fout, "# 1 \"%s\"\n", fname);
  yyin = fin;
  yyout = fout;
  line = 1, scope = 0, para = 0;
  inforeach = 0, foreachscope = 0, foreachpara = 0;
  invardecl = 0;
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
    if (access (path, R_OK|W_OK|X_OK) && mkdir (path, 0700))
      return NULL;
    *s++ = '/';
  }
  return fopen (path, mode);
}

void cleanup (int status)
{
  if (!debug) {
    char command[80] = "rm -r -f ";
    strcat (command, dir);
    int s = system (command);
  }
  exit (status);
}

void compdir (char * file, char ** in, int nin, char * grid)
{
  int i;
  for (i = nin - 1; i >= 0; i--) {
    char * path = in[i];
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
    if (i == 0)
      fprintf (fout, "#include \"grid.h\"\n");
    if (endfor (path, fin, fout))
      cleanup (1);
    fclose (fout);
    free (out);
    fclose (fin);
    free (path);
  }

  /* new variables */
  char * out = malloc (sizeof (char)*(strlen (dir) + strlen ("grid.h") + 2));
  strcpy (out, dir); strcat (out, "/grid.h");
  FILE * fout = fopen (out, "w");
  fprintf (fout,
	   "#include \"common.h\"\n"
	   "int nvar = %d, datasize = %d*sizeof (double);\n", nvar, nvar);
  if (grid)
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
      ;
    else if (!strcmp (argv[i], "-MD"))
      ;
    else if (!strcmp (argv[i], "-debug"))
      debug = 1;
    else if (argv[i][0] != '-' && !strcmp (&argv[i][strlen(argv[i]) - 2], ".c")) {
      if (file) {
	fprintf (stderr, "usage: qcc -grid=[GRID] [OPTIONS] FILE.c\n");
	return 1;
      }
      file = argv[i];
    }
    else {
      strcat (command, " ");
      strcat (command, argv[i]);
    }
  }
  int status;
  if (debug) {
    status = system ("rm -r -f .qcc");
    strcpy (dir, ".qcc");
    status = mkdir (dir, 0700);
  }
  else
    status = (mkdtemp (dir) == NULL);
  if (status) {
    perror (dir);
    return 1;
  }
  if (file) {
    char * out[100], * grid = NULL;
    int nout = includes (argc, argv, out, &grid);
    compdir (file, out, nout, grid);
    strcat (command, " ");
    strcat (command, dir);
    strcat (command, "/");
    strcat (command, file);
  }
  /* compilation */
  status = system (command);
  if (status == -1 ||
      (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || WTERMSIG (status) == SIGQUIT)))
    cleanup (1);
  cleanup (WEXITSTATUS (status));
  return 0;
}
