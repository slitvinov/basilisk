%option noyywrap
%{
  #include <unistd.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <sys/wait.h>
  #include <assert.h>

  enum { scalar, vector };

  typedef int Scalar;
  typedef struct { Scalar x, y; } Vector;
  typedef struct { Vector x, y; } Tensor;

  int debug = 0, fpe = 0;
  char dir[] = ".qccXXXXXX";

  int nvar = 0, nevents = 0;
  int nscalars = 0, nvectors = 0, ntensors = 0;
  Scalar * scalars = NULL;
  Vector * vectors = NULL;
  Tensor * tensors = NULL;
  int line;
  int scope, para, inforeach, foreachscope, foreachpara, 
    inforeach_boundary, inforeach_face;
  int invardecl, vartype;
  int inval, invalpara;
  int brack, inarray;

  int inevent, eventscope, eventpara;
  char eventarray[100];
  int nexpr[100];

  int foreachdim, foreachdimpara, foreachdimline;
  FILE * foreachdimfp;

  char foreachs[80], * fname;
  FILE * foreachfp;
  char reduction[10][4], reductvar[10][80];
  int nreduct;

  int inboundary;
  char * boundaryfunc;
  int nboundary = 0;

  int infunction, functionscope, functionpara;
  FILE * dopen (const char * fname, const char * mode);

  int foreach_face_line;
  enum { face_x, face_y, face_xy };
  int foreach_face_xy;

  typedef struct { char * v; int type, args, scope; } var_t;
  var_t _varstack[100]; int varstack = -1;
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
  
  int identifier (int c) {
    return ((c >= 'a' && c <= 'z') || 
	    (c >= 'A' && c <= 'Z') || 
	    (c >= '0' && c <= '9'));
  }

  var_t * varlookup (char * s, int len) {
    int i;
    for (i = varstack; i >= 0; i--)
      if (strlen(_varstack[i].v) == len && 
	  !strncmp(s, _varstack[i].v, len) &&
	  (s[len-1] != '.' || _varstack[i].type == vector))
	return &_varstack[i];
    return NULL;
  }

  void writefile (FILE * fp, char x, char y, int line1,
		  const char * condition) {
    fputs ("\n"
	   "#undef val\n"
	   "#undef fine\n"
	   "#undef coarse\n"
	   "#undef neighbor\n",
	   yyout);
    if (x == 'x')
      fputs ("#define val(a,k,l) data(k,l)[a]\n"
	     "#define fine(a,k,l) _fine(a,k,l)\n"
	     "#define coarse(a,k,l) _coarse(a,k,l)\n"
	     "#define neighbor(k,l) _neighbor(k,l)\n",
	     yyout);
    else
      fputs ("#define val(a,k,l) data(l,k)[a]\n"
	     "#define fine(a,k,l) _fine(a,l,k)\n"
	     "#define coarse(a,k,l) _coarse(a,l,k)\n"
	     "#define neighbor(k,l) _neighbor(l,k)\n",
	     yyout);
    if (condition)
      fprintf (yyout, "if (%s) {\n", condition);
    fprintf (yyout, "#line %d\n", line1);

    rewind (fp);
    char s[] = "123";
    int c, i = 0;
    while ((c = fgetc(fp)) != EOF) {
      if (i < 3) {
	s[i++] = c; s[i] = '\0';
      }
      else {
	if (s[0] == '.' && !identifier(s[2])) {
	  if      (s[1] == 'x') s[1] = x;
	  else if (s[1] == 'y') s[1] = y;
	}
	fputc (s[0], yyout);
	s[0] = s[1]; s[1] = s[2]; s[2] = c;
      }
    }
    if (i > 0)
      fputs (s, yyout);
    if (condition)
      fputs (" } ", yyout);
  }

  void endforeachdim () {
    foreachdim = 0;
    fclose (yyout);
    yyout = foreachdimfp;
    fputc ('{', yyout);
    FILE * fp = dopen ("dimension.h", "r");
    writefile (fp, 'y', 'x', foreachdimline, NULL);
    writefile (fp, 'x', 'y', foreachdimline, NULL);
    fclose (fp);
    fputc ('}', yyout);
    fprintf (yyout, "\n#line %d\n", line);
  }

  void endforeach () {
    if (inforeach_face) {
      fclose (yyout);
      yyout = foreachfp;
      fputs ("foreach()", yyout);
      FILE * fp = dopen ("foreach_face.h", "r");
      if (foreach_face_xy == face_xy) {
	fputs (" { int jg = -1; VARIABLES; ", yyout);
	writefile (fp, 'y', 'x', foreach_face_line, "is_face_x()");
	fputs (" } { int ig = -1; VARIABLES; ", yyout);
	writefile (fp, 'x', 'y', foreach_face_line, "is_face_x()");
	fputs (" } ", yyout);
      }
      else if (foreach_face_xy == face_x) {
	fputs (" { int ig = -1; VARIABLES; ", yyout);
	writefile (fp, 'x', 'y', foreach_face_line, "is_face_x()");
	fputs (" } ", yyout);
      }
      else {
	fputs (" { int jg = -1; VARIABLES; ", yyout);
	writefile (fp, 'x', 'y', foreach_face_line, "is_face_y()");
	fputs (" } ", yyout);
      }
      fputs (" end_foreach()\n", yyout);
      if (foreach_face_xy != face_x) {
	fputs ("boundary_ghost (top, ({", yyout);
	if (foreach_face_xy == face_xy)
	  writefile (fp, 'y', 'x', foreach_face_line, NULL);
	else
	  writefile (fp, 'x', 'y', foreach_face_line, NULL);
	fputs ("}));\n", yyout);
      }
      if (foreach_face_xy != face_y) {
	fputs ("boundary_ghost (right, ({", yyout);
	writefile (fp, 'x', 'y', foreach_face_line, NULL);
	fputs ("}));\n", yyout);
      }
      fprintf (yyout, "#line %d\n", line);
      fclose (fp);
    }
    else if (nreduct > 0) {
      fputs ("\n#undef _OMPEND\n#define _OMPEND ", yyout);
      int i;
      for (i = 0; i < nreduct; i++)
	fprintf (yyout, "OMP(omp critical) if (_%s %s %s) %s = _%s; ",
		 reductvar[i], strcmp(reduction[i], "min") ? ">" : "<",
		 reductvar[i], reductvar[i], reductvar[i]);
      fprintf (yyout,
	       "\nend_%s();\n"
	       "#undef _OMPSTART\n"
	       "#undef _OMPEND\n"
	       "#define _OMPSTART\n"
	       "#define _OMPEND\n"
	       "#line %d\n",
	       foreachs, line);
    }
    else
      fprintf (yyout, " end_%s();", foreachs);
    inforeach = inforeach_boundary = inforeach_face = 0;
  }

  void endevent() {
    fprintf (yyout, "  return 0; } ");
    inevent = 0;
    nevents++;
  }

  void infunction_declarations() {
    if (debug)
      fprintf (stderr, "%s:%d: function declarations\n", fname, line);
    fputs (" int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED (jg);"
	   " POINT_VARIABLES; ", yyout);
  }

#define nonspace(s) { while (strchr(" \t\v\n\f", *s)) s++; }

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
  if (inforeach == 1 && scope == foreachscope && para == foreachpara) {
    ECHO;
    if (inforeach_face) {
      fclose (yyout);
      yyout = dopen ("foreach_face.h", "w");
      foreach_face_line = line;
      foreach_face_xy = face_xy;
      FILE * fp = dopen ("foreach.h", "r");
      int c;
      while ((c = fgetc (fp)) != EOF)
	if (c == 'x')
	  foreach_face_xy = face_x;
	else if (c == 'y')
	  foreach_face_xy = face_y;
      fclose (fp);
    }
    else {
      fclose (yyout);
      yyout = foreachfp;
      if (nreduct > 0) {
	fprintf (yyout, "\n#undef _OMPSTART\n#define _OMPSTART ");
	int i;
	for (i = 0; i < nreduct; i++)
	  fprintf (yyout, "double _%s = %s; ", reductvar[i], reductvar[i]);
	fprintf (yyout, "\n#line %d\n", line);
      }
      FILE * fp = dopen ("foreach.h", "r");
      int c;
      while ((c = fgetc (fp)) != EOF)
	fputc (c, yyout);
      fclose (fp);
    }
    inforeach = 2;
  }
  else if (inevent > 0 && inevent < 4 && 
	   scope == eventscope && eventpara == para + 1) {
    if (!eventarray[nevents])
      fprintf (yyout, "); "
	       "  *ip = i; *tp = t; "
	       "  return ret; "
	       "} ");
    fprintf (yyout, 
	     "static int event_%d (int i, double t) { ",
	     nevents);
    assert (nevents < 100);
    nexpr[nevents] = inevent;
    inevent = 4;
  }
  else
    ECHO;
}

"{" {
  ECHO;
  if (infunction && functionpara == 1 && scope == functionscope)
    infunction_declarations();
  scope++;
}

"}" {
  ECHO; scope--;
  if (scope < 0)
    return yyerror ("mismatched '}'");
  varpop();
  if (foreachdim && scope == foreachdim)
    endforeachdim ();
  if (infunction && scope == functionscope)
    infunction = 0;
  if (inforeach && scope == foreachscope)
    endforeach ();
  else if (inevent && scope == eventscope)
    endevent ();
}

foreach{ID}* {
  strcpy (foreachs, yytext);
  inforeach = 1; foreachscope = scope; foreachpara = para;
  nreduct = 0;
  foreachfp = yyout;
  yyout = dopen ("foreach.h", "w");
  inforeach_boundary = (!strncmp(foreachs, "foreach_boundary", 16));
  inforeach_face = (!strcmp(foreachs, "foreach_face"));
  ECHO;
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
  if (infunction && scope == functionscope) {
    if (scope > 0) {
      fputs ("; ", yyout);
      infunction_declarations();
    }
    infunction = 0;
  }
  if (inboundary) {
    ECHO;
    fputs (" } ", yyout);
    if (scope == 0)
      /* file scope */
      fprintf (yyout, "static void _boundary%d (void) { %s } ",
	       nboundary++, boundaryfunc);
    else
      /* function scope */
      fputs (boundaryfunc, yyout);
    free (boundaryfunc);
    inboundary = inforeach_boundary = inforeach = 0;
  }
  else if (inforeach && scope == foreachscope && para == foreachpara) {
    ECHO;
    endforeach ();
  }
  else if (inevent > 0 && inevent < 3 && para == eventpara) {
    if (eventarray[nevents])
      return yyerror ("cannot mix arrays and expressions");
    fprintf (yyout, "); "
	     "  *ip = i; *tp = t; "
	     "  return ret; "
	     "} "
	     "static int event_expr%d%d (int * ip, double * tp) { "
	     "  int i = *ip; double t = *tp; "
	     "  int ret = (", nevents, inevent++);
  }
  else if (inevent == 4 && scope == eventscope && para == eventpara - 1) {
    ECHO;
    endevent ();
  }
  else
    ECHO;
  invardecl = 0;
}

[^{ID}](scalar|vector|tensor){WS}+[a-zA-Z0-9\[\]]+ {
  ECHO;
  if (yytext[0] == '(') para++;
  char * var = strstr(yytext,"scalar");
  vartype = scalar;
  if (!var) {
    var = strstr(yytext,"vector");
    if (!var)
      var = strstr(yytext,"tensor");
    vartype = vector;
  }
  var = &var[7];
  nonspace (var);
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
    nonspace (var);
    if (debug)
      fprintf (stderr, "%s:%d: declaration: %s\n", fname, line, var);
    varpush (var, vartype, scope);
  }
  else
    REJECT;
  ECHO;
}

new{WS}+scalar {
  if (scope > 0)
    fprintf (yyout, "new_scalar (%d)", nvar);
  else {
    fprintf (yyout, "%d", nvar);
    nscalars++;
    scalars = realloc (scalars, sizeof (Scalar)*nscalars);
    scalars[nscalars-1] = nvar;
  }
  nvar++;
}

new{WS}+vector {
  if (scope > 0)
    fprintf (yyout, "new_vector ((vector){%d,%d})", nvar, nvar + 1);
  else {
    fprintf (yyout, "{%d,%d}", nvar, nvar + 1);
    nvectors++;
    vectors = realloc (vectors, sizeof (Vector)*nvectors);
    vectors[nvectors-1] = (Vector){nvar,nvar + 1};
  }
  nvar += 2;
}

new{WS}+tensor {
  if (scope > 0)
    fprintf (yyout, "new_tensor ((tensor){{%d,%d},{%d,%d}})",
	     nvar, nvar + 1, nvar + 2, nvar + 3);
  else {
    fprintf (yyout, "{{%d,%d},{%d,%d}}",
	     nvar, nvar + 1, nvar + 2, nvar + 3);
    ntensors++;
    tensors = realloc (tensors, sizeof (Tensor)*ntensors);
    tensors[ntensors-1] = (Tensor){{nvar,nvar + 1},{nvar + 2, nvar + 3}};
  }  
  nvar += 4;
}

[^{ID}]val{WS}*[(]    {
  if (yytext[0] == '(') para++;
  inval = 1; invalpara = para++;
  ECHO;
}

[a-zA-Z_0-9\.]+{WS}*\[{WS}*. {
  /* v[... */
  var_t * var = NULL;
  if ((inforeach || infunction) && !inval) {
    char * s = yytext;
    while (!strchr(" \t\v\n\f[.", *s)) s++;
    if ((var = varlookup (yytext, s - yytext))) {
      s = yytext;
      while (!strchr(" \t\v\n\f[", *s)) s++;
      *s = '\0';
      if (yytext[yyleng-1] == ']')
	/* v[] */
	fprintf (yyout, "val(%s,0,0)", yytext);
      else {
	fprintf (yyout, "val(%s", yytext);
	if (var->args > 0) {
	  /* v[...][... */
	  fputc ('[', yyout);
	  fputc (yytext[yyleng-1], yyout);
	  inarray = brack++;
	  int j = var->args;
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
	  unput (yytext[yyleng-1]);
	  inarray = ++brack;
	}
      }
    }
  }
  if (!var)
    REJECT;
}

[a-zA-Z_0-9\.]+{WS}*\[{WS}*(left|right|top|bottom){WS}*\]{WS}*= {
  /* v[top] = */
  char * s = yytext;
  while (!strchr(" \t\v\n\f[.", *s)) s++;
  if (varlookup (yytext, s - yytext)) {
    s = yytext;
    while (!strchr(" \t\v\n\f[", *s)) s++;
    *s++ = '\0';
    nonspace (s);
    char * b = s;
    while (!strchr(" \t\v\n\f]", *s)) s++;
    *s++ = '\0';
    char * func = malloc ((strlen (yytext) + 1)*sizeof (char));
    strcpy (func, yytext);
    char * s1 = func;
    while (*s1 != '\0') {
      if (*s1 == '.')
	*s1 = '_';
      s1++;
    }
    boundaryfunc = malloc ((strlen ("_boundary[][] = _;") + 
			    2*strlen (b) + 
			    strlen(yytext) + 
			    strlen (func) + 1)*sizeof (char));
    sprintf (boundaryfunc, "_boundary[%s][%s] = _%s%s;", b, yytext, func, b);
    fprintf (yyout, 
	     "double _%s%s (Point point, scalar s) {"
	     " int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED (jg);"
	     " POINT_VARIABLES; "
	     " return ", 
	     func, b);
    free (func);
    inboundary = inforeach_boundary = 1;
    inforeach = 2;
  }
  else
    REJECT;
}

ghost {
  if (inforeach_boundary)
    fputs ("ig,jg", yyout);
  else
    ECHO;
}

[^{ID}]Point{WS}+point[^{ID}] {
  /* Point point */
  if (yytext[0] == '(') para++;
  infunction = 1; functionscope = scope; functionpara = para;
  if (yytext[yyleng-1] == ')') para--;
  ECHO;
}

for{WS}*[(]{WS}*(scalar|vector){WS}+{ID}+{WS}+in{WS}+{ID}+{WS}*[)] {
  /* for (scalar .. in .. ) */
  char * s = strchr (&yytext[3], 'r'); s++;
  int vartype = s[-2] == 'o' ? vector : scalar;
  nonspace (s);
  char * id = s;
  while (!strchr (" \t\v\n\f", *s)) s++;
  *s++ = '\0';
  s = strchr (s, 'n'); s++;
  nonspace (s);
  char * list = s;
  while (!strchr (" \t\v\n\f)", *s)) s++;
  *s++ = '\0';
  static int i = 0;
  if (vartype == scalar)
    fprintf (yyout,
	     "for (scalar %s = *%s, *_i%d = %s; %s >= 0; %s = *++_i%d)",
	     id, list, i, list, id, id, i);
  else
    fprintf (yyout,
	     "for (vector %s = *%s, *_i%d = %s; %s.x >= 0; %s = *++_i%d)",
	     id, list, i, list, id, id, i);
  i++;
  varpush (id, vartype, scope);
}

for{WS}*[(][^)]+{WS}+in{WS}+[^)]+[)] {
  /* for (a,b in c,d) */
  char * id[10], * list[10];
  int nid = 0, nlist = 0, inin = 0;
  char * s = strchr (yytext, '('); s++;
  s = strtok (s, " \t\v\n\f,");
  while (s) {
    if (!strcmp (s, "in"))
      inin = 1;
    else if (inin)
      list[nlist++] = s;
    else
      id[nid++] = s;
    s = strtok (NULL, " \t\v\n\f,)");
  }
  if (nlist != nid)
    return yyerror ("lists must be the same size");
  int i;
  static int index = 0;
  char * sc = NULL, * vc = NULL;
  for (i = 0; i < nid; i++) {
    var_t * var = varlookup (id[i], strlen(id[i]));
    if (!var) {
      fprintf (stderr, "%s:%d: error: '%s' is not a scalar or vector\n", 
	       fname, line, id[i]);
      return 1;
    }
    fprintf (yyout, "%s * _i%d = %s; ", 
	     var->type == vector ? "vector" : "scalar",
	     index + i, list[i]);
    if (!sc && var->type == scalar) sc = id[i];
    if (!vc && var->type == vector) vc = id[i];
  }
  fprintf (yyout, "for (%s = *%s", id[0], list[0]);
  for (i = 1; i < nid; i++)
    fprintf (yyout, ", %s = *%s", id[i], list[i]);
  if (sc)
    fprintf (yyout, "; %s >= 0; ", sc);
  else
    fprintf (yyout, "; %s.x >= 0; ", vc);
  fprintf (yyout, "%s = *++_i%d", id[0], index);
  for (i = 1; i < nid; i++)
    fprintf (yyout, ", %s = *++_i%d", id[i], index + i);
  fputc (')', yyout);
  index += nid;
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
	   "event_expr%d%d (int * ip, double * tp) {"
	   "  int i = *ip; double t = *tp;"
	   "  int ret = (", nevents, inevent++);
  eventscope = scope; eventpara = ++para;
  eventarray[nevents] = 0;
}

[it]{WS}*={WS}*[{][^}]*[}] {
  if (inevent == 1) {
    eventarray[nevents] = yytext[0];
    yytext[yyleng-1] = '\0';
    fprintf (yyout, "1); "
	     "  *ip = i; *tp = t; "
	     "  return ret; "
	     "} "
	     "static %s event_array%d[] = %s,-1}; ",
	     yytext[0] == 'i' ? "int" : "double", 
	     nevents, strchr (yytext, '{'));
  }
  else
    REJECT;
}

foreach_dimension{WS}*[(]{WS}*[)] {
  foreachdimline = line;
  foreachdim = scope; foreachdimpara = para;
  foreachdimfp = yyout;
  yyout = dopen ("dimension.h", "w");
}

reduction[(](min|max):{ID}+[)] {
  char * s = strchr (yytext, '('), * s1 = strchr (yytext, ':');
  *s1 = '\0'; s1++;
  assert (nreduct < 10);
  strcpy (reduction[nreduct], ++s);
  yytext[yyleng-1] = '\0';
  strcpy (reductvar[nreduct++], s1);
}

{ID}+ {
  if (inforeach) {
    int i;
    for (i = 0; i < nreduct; i++)
      if (!strcmp (yytext, reductvar[i])) {
	fputc ('_', yyout);
	break;
      }
  }
  ECHO;
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
int includes (int argc, char ** argv, char ** out, 
	      char ** grid, int * default_grid);

int endfor (char * file, FILE * fin, FILE * fout)
{
  fname = stripslash (file);
  fprintf (fout, "# 1 \"%s\"\n", fname);
  yyin = fin;
  yyout = fout;
  line = 1, scope = para = 0;
  inforeach = foreachscope = foreachpara = 
    inforeach_boundary = inforeach_face = 0;
  invardecl = 0;
  inval = invalpara = 0;
  brack = inarray = 0;
  inevent = 0;
  foreachdim = 0;
  inboundary = 0;
  infunction = 0;
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

FILE * dopen (const char * fname, const char * mode)
{
  char * out = malloc (sizeof (char)*(strlen (dir) + strlen (fname) + 2));
  strcpy (out, dir); strcat (out, "/"); strcat (out, fname);
  FILE * fout = fopen (out, mode);
  free (out);
  return fout;
}

void compdir (char * file, char ** in, int nin, char * grid, int default_grid)
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
    if (i == 0) {
      if (fpe)
	fputs ("#define _GNU_SOURCE 1\n"
	       "#include <fenv.h>\n", fout);
      fputs ("#include \"grid.h\"\n", fout);
    }
    if (endfor (path, fin, fout))
      cleanup (1);
    fclose (fout);
    free (out);
    fclose (fin);
    free (path);
  }

  FILE * fout = dopen ("grid.h", "w");
  /* new variables */
  fprintf (fout,
	   "#include \"common.h\"\n"
	   "int nvar = %d, datasize = %d*sizeof (double);\n"
	   "scalar all[] = {",
	   nvar, nvar);
  for (i = 0; i < nvar; i++)
    fprintf (fout, "%d,", i);
  fputs ("-1};\n", fout);
  /* events */
  int j;
  for (i = 0; i < nevents; i++) {
    fprintf (fout, "static int event_%d (int i, double t);\n", i);
    for (j = 0; j < nexpr[i]; j++)
      fprintf (fout,
	       "static int event_expr%d%d (int * ip, double * tp);\n",
	       i, j);
    if (eventarray[i])
      fprintf (fout, "static %s event_array%d[];\n", 
	       eventarray[i] == 'i' ? "int" : "double", i);
  }
  fputs ("Event Events[] = {\n", fout);
  for (i = 0; i < nevents; i++) {
    fprintf (fout, "  { false, %d, event_%d, {", nexpr[i], i);
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
  /* boundaries */
  for (i = 0; i < nboundary; i++)
    fprintf (fout, "static void _boundary%d (void);\n", i);
  /* methods */
  fprintf (fout,
	   "scalar %s_new_scalar (scalar);\n"
	   "#define new_scalar %s_new_scalar\n"
	   "vector %s_new_vector (vector);\n"
	   "#define new_vector %s_new_vector\n"
	   "tensor %s_new_tensor (tensor);\n"
	   "#define new_tensor %s_new_tensor\n",
	   grid, grid, grid, grid, grid, grid);
  fputs ("static void init_solver (void) {\n"
	 "  for (int b = 0; b < nboundary; b++)\n"
	 "    boundary[b] = calloc (nvar, sizeof (BoundaryFunc));\n",
	 fout);
  /* refinement functions */
  fputs ("  refine = calloc (nvar, sizeof (RefineFunc));\n", fout);
  if (fpe)
    fputs ("  feenableexcept (FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);\n", fout);
  for (i = 0; i < nscalars; i++)
    fprintf (fout, "  new_scalar (%d);\n", scalars[i]);
  for (i = 0; i < nvectors; i++)
    fprintf (fout, "  new_vector ((vector){%d,%d});\n", 
	     vectors[i].x, vectors[i].y);
  for (i = 0; i < ntensors; i++)
    fprintf (fout, "  new_tensor ((tensor){{%d,%d},{%d,%d}});\n", 
	     tensors[i].x.x, tensors[i].x.y,
	     tensors[i].y.x, tensors[i].y.y);
  for (i = 0; i < nboundary; i++)
    fprintf (fout, "  _boundary%d();\n", i);
  fputs ("  init_events();\n}\n", fout);
  /* grid */
  if (default_grid)
    fprintf (fout, "#include \"grid/%s.h\"\n", grid);
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
    else if (!strcmp (argv[i], "-fpe"))
      fpe = 1;
    else if (argv[i][0] != '-' && 
	     !strcmp (&argv[i][strlen(argv[i]) - 2], ".c")) {
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
    int default_grid, nout = includes (argc, argv, out, &grid, &default_grid);
    compdir (file, out, nout, grid, default_grid);
    strcat (command, " ");
    strcat (command, dir);
    strcat (command, "/");
    strcat (command, file);
  }
  /* compilation */
  if (debug)
    fprintf (stderr, "command: %s\n", command);
  status = system (command);
  if (status == -1 ||
      (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				WTERMSIG (status) == SIGQUIT)))
    cleanup (1);
  cleanup (WEXITSTATUS (status));
  return 0;
}
