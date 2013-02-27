%option noyywrap
%{
  #include <unistd.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <sys/wait.h>
  #include <assert.h>

  enum { variable, vector };

  int debug = 0;
  char dir[] = ".qccXXXXXX";

  int nvar = 0, nevents = 0;
  int line;
  int scope, para, inforeach, foreachscope, foreachpara;
  int invardecl, vartype;
  int inval, invalpara;
  int brack, inarray;

  int inevent, eventscope, eventpara;
  char eventarray[100];
  int nexpr[100];

  int foreachdim, foreachdimpara, foreachdimline;
  char foreachdimname[80];
  FILE * foreachdimfp;

  char foreachs[80], * fname;
  
  struct { char * v; int type, args, scope; } _varstack[100]; int varstack = -1;
  void varpush (const char * s, int type, int scope) {
    if (s[0] != '\0') {
      char * f = malloc (sizeof (char)*(strlen (s) + 1));
      strcpy (f, s);
      char * q = f;
      int na = 0;
      while ((q = strchr (q, '['))) {
	*q++ = '\0'; na++;
      }
      _varstack[++varstack].v = f;
      _varstack[varstack].scope = scope;
      _varstack[varstack].args = na;
      _varstack[varstack].type = type;
    }
  }
  void varpop () {
    while (varstack >= 0 && _varstack[varstack].scope > scope)
      free (_varstack[varstack--].v);
  }
  
  void endforeach () {
    inforeach = 0;
    fprintf (yyout, " end_%s();", foreachs);
  }

  int identifier (int c) {
    return ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9'));
  }

  void writefile (FILE * fp, char x, char y) {
    rewind (fp);
    char s[] = "123";
    int c, i = 0;
    while ((c = fgetc(fp)) != EOF) {
      if (i < 3) {
	s[i++] = c; s[i] = '\0';
      }
      else {
	if (!identifier(s[0]) && !identifier(s[2])) {
	  if      (s[1] == 'x') s[1] = x;
	  else if (s[1] == 'y') s[1] = y;
	}
	fputc (s[0], yyout);
	s[0] = s[1]; s[1] = s[2]; s[2] = c;
      }
    }
    if (i > 0)
      fputs (s, yyout);
  }

  void endforeachdim () {
    foreachdim = 0;
    fclose (yyout);
    yyout = foreachdimfp;
    FILE * fp = fopen (foreachdimname, "r");
    writefile (fp, 'x', 'y');
    fprintf (yyout, 
	     "\n#undef val1\n"
	     "#define val1(a,i,j) val(a,j,i)\n"
	     "#line %d\n", foreachdimline);
    writefile (fp, 'y', 'x');
    fclose (fp);
  }

  void endevent() {
    fprintf (yyout, "\n  return 0;\n}\n#line %d\n", line);
    inevent = 0;
    nevents++;
  }

#define YY_INPUT(buf,result,max_size)			      \
  {							      \
    int c = fgetc(yyin);				      \
    result = (c == EOF) ? YY_NULL : (buf[0] = c, 1);	      \
    if (c == '\n') line++;				      \
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
  para--; if (para == invalpara) inval = 0;
  if (para < 0)
    return yyerror ("mismatched ')'");
  if (inevent > 0 && inevent < 4 && scope == eventscope && eventpara == para + 1) {
    if (!eventarray[nevents])
      fprintf (yyout, ");\n"
	       "  *ip = i; *tp = t;\n"
	       "  return ret;\n"
	       "}\n");
    fprintf (yyout, 
	     "static int event_action%d (int i, double t) {\n"
	     "  #line %d\n",
	     nevents, line);
    assert (nevents < 100);
    nexpr[nevents] = inevent;
    inevent = 4;
  }
  else
    ECHO;
}

"{" {
  ECHO; scope++;
}

"}" {
  ECHO; scope--;
  if (scope < 0)
    return yyerror ("mismatched '}'");
  varpop();
  if (foreachdim && scope == foreachdim)
    endforeachdim ();
  if (inforeach && scope == foreachscope)
    endforeach (line);
  else if (inevent && scope == eventscope)
    endevent ();
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
    if (c == ')') { para--; if (para == invalpara) inval = 0; }
  }
  if (c != ')')
    return yyerror ("expecting ')");
  if (varstack >= 0)
    fprintf (yyout, 
	     "\n#undef val1\n"
	     "#define val1(a,i,j) val(a,i,j)\n"
	     "#line %d\n", line);
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
  if (foreachdim && scope == foreachdim && para == foreachdimpara) {
    ECHO;
    endforeachdim ();
  }
  if (inforeach && scope == foreachscope && para == foreachpara) {
    ECHO;
    endforeach (line - 1);
  }
  else if (inevent > 0 && inevent < 3 && para == eventpara) {
    if (eventarray[nevents])
      return yyerror ("cannot mix arrays and expressions");
    fprintf (yyout, ");\n"
	     "  *ip = i; *tp = t;\n"
	     "  return ret;\n"
	     "}\n"
	     "static int event_expr%d%d (int * ip, double * tp) {\n"
	     "  int i = *ip; double t = *tp;\n"
	     "  #line %d\n"
	     "  int ret = (", nevents, inevent++, line);
  }
  else if (inevent == 4 && scope == eventscope && para == eventpara - 1) {
    ECHO;
    endevent ();
  }
  else
    ECHO;
  invardecl = 0;
}

[^{ID}]{WS}*(scalar|vector){WS}+[a-zA-Z0-9\[\]]+ {
  ECHO;
  if (yytext[0] == '(') para++;
  char * var = strstr(yytext,"scalar");
  vartype = variable;
  if (!var) {
    var = strstr(yytext,"vector");
    vartype = vector;
  }
  var = &var[7];
  while (strchr (" \t\v\n\f", *var)) var++;
  if (para == 0) { /* declaration */
    if (debug)
      fprintf (stderr, "%s:%d: declaration: %s\n", fname, line, var);
    varpush (var, vartype, scope);
    invardecl = scope + 1;
  }
  else if (para == 1) { /* function prototype (no nested functions) */
    if (debug)
      fprintf (stderr, "%s:%d: proto: %s\n", fname, line, var);
    varpush (var, vartype, scope + 1);
  }
}

,{WS}*{ID}+ {
  if (invardecl == scope + 1) {
    char * var = &yytext[1];
    while (strchr (" \t\v\n\f", *var)) var++;
    if (debug)
      fprintf (stderr, "%s:%d: declaration: %s\n", fname, line, var);
    varpush (var, vartype, scope);
  }
  else
    REJECT;
  ECHO;
}

new{WS}+scalar {
  fprintf (yyout, "%d", nvar++);
}

new{WS}+vector {
  fprintf (yyout, "{%d,%d}", nvar++, nvar++);
}

[^{ID}]val{WS}*[(]    {
  if (yytext[0] == '(') para++;
  inval = 1; invalpara = para++;
  ECHO;
}

[a-zA-Z_0-9\.]+{WS}*\[{WS}*. {
  /* v[... */
  int found = 0;
  if (inforeach && !inval) {
    char * s = yytext;
    while (!strchr(" \t\v\n\f[.", *s)) s++;
    int i, len = s - yytext;
    for (i = varstack; i >= 0 && !found; i--)
      if (strlen(_varstack[i].v) == len && !strncmp(yytext, _varstack[i].v, len) &&
	  (*s != '.' || _varstack[i].type == vector)) {
	found = 1;
	s = yytext;
	while (!strchr(" \t\v\n\f[", *s)) s++;
	*s = '\0';
	if (yytext[yyleng-1] == ']')
	  /* v[] */
	  fprintf (yyout, "val1(%s,0,0)", yytext);
	else {
	  fprintf (yyout, "val1(%s", yytext);
	  if (_varstack[i].args > 0) {
	    /* v[...][... */
	    fputc ('[', yyout);
	    fputc (yytext[yyleng-1], yyout);
	    inarray = brack++;
	    int j = _varstack[i].args;
	    while (j) {
	      int c = getput();
	      if (c == EOF)
		return yyerror ("unexpected EOF");
	      if (c == '[') brack++;
	      if (c == ']') {
		brack--;
		if (brack == inarray)
		  j--;
	      }
	    }
	    int c = input();
	    if (c != '[') {
	      fprintf (stderr, "%s:%d: error: expecting '[' not '", fname, line);
	      fputc (c, stderr);
	      fputs ("'\n", stderr);
	      return 1;
	    }
	    fputc (',', yyout);
	    c = input();
	    if (c == ']')
	      /* v[...][] */
	      fputs ("0,0)", yyout);
	    else {
	      fputc (c, yyout);
	      inarray = ++brack;
	    }
	  }
	  else {
	    /* v[...] */
	    fputc (',', yyout);
	    fputc (yytext[yyleng-1], yyout);
	    inarray = ++brack;
	  }
	}
      }
  }
  if (!found)
    REJECT;
}

[\n]       ECHO;
"["        ECHO; brack++;
"]"        {
  if (inarray == brack) {
    fputc (')', yyout);
    inarray = 0;
  }
  else
    ECHO;
  brack--;
  if (brack < 0)
    return yyerror ("mismatched ']'");  
}

[^{ID}]event{WS}*[(] {
  fputc (yytext[0], yyout);
  /* event (... */
  fprintf (yyout, 
	   "static int event_expr%d%d (int * ip, double * tp) {\n"
	   "  int i = *ip; double t = *tp;\n"
	   "  #line %d\n"
	   "  int ret = (", nevents, inevent++, line);
  eventscope = scope; eventpara = ++para;
  eventarray[nevents] = 0;
}

[it]{WS}*={WS}*[{][^}]*[}] {
  if (inevent == 1) {
    eventarray[nevents] = yytext[0];
    yytext[yyleng-1] = '\0';
    fprintf (yyout, "1);\n"
	     "  *ip = i; *tp = t;\n"
	     "  return ret;\n"
	     "}\n"
	     "static %s event_array%d[] = %s,-1};\n",
	     yytext[0] == 'i' ? "int" : "double", nevents, strchr (yytext, '{'));
  }
  else
    REJECT;
}

foreach_dimension{WS}*[(]{WS}*[)] {
  if (!inforeach)
    ECHO;
  else {
    foreachdimline = line;
    foreachdim = scope; foreachdimpara = para;
    foreachdimfp = yyout;
    strcpy (foreachdimname, dir);
    strcat (foreachdimname, "/dimension.h");
    yyout = fopen (foreachdimname, "w");
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

int endfor (char * file, FILE * fin, FILE * fout)
{
  fname = stripslash (file);
  fprintf (fout, "# 1 \"%s\"\n", fname);
  yyin = fin;
  yyout = fout;
  line = 1, scope = para = 0;
  inforeach = foreachscope = foreachpara = 0;
  invardecl = 0;
  inval = invalpara = 0;
  brack = inarray = 0;
  inevent = 0;
  foreachdim = 0;
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
    int s = system (command); s = s;
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

  char * out = malloc (sizeof (char)*(strlen (dir) + strlen ("grid.h") + 2));
  strcpy (out, dir); strcat (out, "/grid.h");
  FILE * fout = fopen (out, "w");
  /* new variables */
  fprintf (fout,
	   "#include \"common.h\"\n"
	   "int nvar = %d, datasize = %d*sizeof (double);\n", nvar, nvar);
  /* events */
  int j;
  for (i = 0; i < nevents; i++) {
    fprintf (fout, "static int event_action%d (int i, double t);\n", i);
    for (j = 0; j < nexpr[i]; j++)
      fprintf (fout,
	       "static int event_expr%d%d (int * ip, double * tp);\n",
	       i, j);
    if (eventarray[i])
      fprintf (fout, "static %s event_array%d[];\n", eventarray[i] == 'i' ? "int" : "double", i);
  }
  fputs ("Event Events[] = {\n", fout);
  for (i = 0; i < nevents; i++) {
    fprintf (fout, "  { false, %d, event_action%d, {", nexpr[i], i);
    for (j = 0; j < nexpr[i] - 1; j++)
      fprintf (fout, "event_expr%d%d, ", i, j);
    fprintf (fout, "event_expr%d%d}, ", i, j);
    if (eventarray[i] == 'i')
      fprintf (fout, "event_array%d, ", i);
    else
      fprintf (fout, "NULL, ");
    if (eventarray[i] == 't')
      fprintf (fout, "event_array%d},\n", i);
    else
      fprintf (fout, "NULL},\n");
  }
  fputs ("  { true }\n};\n", fout);
  /* grid */
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
