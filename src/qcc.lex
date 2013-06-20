%option noyywrap
%{
  #include <unistd.h>
  #include <string.h>
  #include <ctype.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <sys/wait.h>
  #include <assert.h>

  enum { scalar, vector, tensor };

  typedef struct { int i; char * name; } Scalar;
  typedef struct { int x, y; char * name; } Vector;
  typedef struct { Vector x, y; char * name; } Tensor;

  int debug = 0, fpe = 0, catch = 0;
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
  int inreturn;

  #define EVMAX 100
  int inevent, eventscope, eventpara;
  char eventarray[EVMAX], * eventfile[EVMAX], * eventid[EVMAX];
  char * eventarray_elems[EVMAX];
  int nexpr[EVMAX], eventline[EVMAX];

  int foreachdim, foreachdimpara, foreachdimline;
  FILE * foreachdimfp;

  #define REDUCTMAX 10
  char foreachs[80], * fname;
  FILE * foreachfp;
  char reduction[REDUCTMAX][4], reductvar[REDUCTMAX][80];
  int nreduct;

  int inboundary;
  char * boundaryfunc, boundaryvar[80];
  int nboundary = 0;

  int infunction, infunctiondecl, functionscope, functionpara;
  FILE * dopen (const char * fname, const char * mode);

  int foreach_face_line;
  enum { face_x, face_y, face_xy };
  int foreach_face_xy;

  char ** args = NULL, ** argss = NULL;
  int nargs = 0, inarg;

  typedef struct { 
    char * v; 
    int type, args, scope, automatic; 
    int i[4];
  } var_t;
  var_t _varstack[100]; int varstack = -1;
  void varpush (const char * s, int type, int scope) {
    if (s[0] != '\0') {
      char * f = malloc (strlen (s) + 1);
      strcpy (f, s);
      char * q = f;
      int na = 0;
      while ((q = strchr (q, '['))) {
	*q++ = '\0'; na++;
      }
      _varstack[++varstack] = (var_t) { f, type, na, scope, 0, {-1} };
    }
  }

  char * makelist (const char * input, int type);

  void delete_automatic (int scope) {
    char * list = NULL;
    for (int i = varstack; i >= 0 && _varstack[i].scope > scope; i--) {
      var_t var = _varstack[i];
      if (var.automatic) {
	if (list == NULL) {
	  list = malloc (strlen (var.v) + 3);
	  strcpy (list, "{");
	}
	else {
	  list = realloc (list, strlen (list) + strlen (var.v) + 3);
	  strcat (list, ",");
	}
	strcat (list, var.v);
      }
    }
    if (list) {
      strcat (list, "}");
      char * slist = makelist (list, scalar);
      if (debug)
	fprintf (stderr, "%s:%d: deleting %s\n", fname, line, list);
      fprintf (yyout, " delete (%s); ", slist);
      free (slist);
      free (list);
    }
  }

  void varpop () {
    delete_automatic (scope);
    while (varstack >= 0 && _varstack[varstack].scope > scope)
      free (_varstack[varstack--].v);
  }

  var_t * vartop() {
    return varstack >= 0 ? &_varstack[varstack] : NULL;
  }

  int identifier (int c) {
    return ((c >= 'a' && c <= 'z') || 
	    (c >= 'A' && c <= 'Z') || 
	    (c >= '0' && c <= '9'));
  }

  var_t * varlookup (char * s, int len) {
    for (int i = varstack; i >= 0; i--)
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
	fputs ("foreach_boundary_ghost (top) {\n", yyout);
	if (foreach_face_xy == face_xy)
	  writefile (fp, 'y', 'x', foreach_face_line, NULL);
	else
	  writefile (fp, 'x', 'y', foreach_face_line, NULL);
	fputs ("} end_foreach_boundary_ghost();\n", yyout);
	fputs ("#ifdef foreach_boundary_ghost_halo\n"
	       "foreach_boundary_ghost_halo (top) {\n", yyout);
	if (foreach_face_xy == face_xy)
	  writefile (fp, 'y', 'x', foreach_face_line, NULL);
	else
	  writefile (fp, 'x', 'y', foreach_face_line, NULL);
	fputs ("} end_foreach_boundary_ghost_halo();\n"
	       "#endif\n", yyout);
      }
      if (foreach_face_xy != face_y) {
	fputs ("foreach_boundary_ghost (right) {\n", yyout);
	writefile (fp, 'x', 'y', foreach_face_line, NULL);
	fputs ("} end_foreach_boundary_ghost();\n", yyout);
	fputs ("#ifdef foreach_boundary_ghost_halo\n"
	       "foreach_boundary_ghost_halo (right) {\n", yyout);
	writefile (fp, 'x', 'y', foreach_face_line, NULL);
	fputs ("} end_foreach_boundary_ghost_halo();\n"
	       "#endif\n", yyout);
      }
      fprintf (yyout, "#line %d\n", line);
      fclose (fp);
    }
    if (nreduct > 0)
      fprintf (yyout,
	       "\nend_%s();\n"
	       "#undef _OMPSTART\n"
	       "#undef _OMPEND\n"
	       "#define _OMPSTART\n"
	       "#define _OMPEND\n"
	       "#line %d\n",
	       foreachs, line);
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
    if (!infunctiondecl) {
      if (debug)
	fprintf (stderr, "%s:%d: function declarations\n", fname, line);
      fputs (" int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED (jg);"
	     " POINT_VARIABLES; ", yyout);
      infunctiondecl = 1;
    }
  }

  char * boundaryrep (char * name) {
    if (!inboundary || strcmp(name, boundaryvar))
      return name;
    static char s[] = "_s";
    return s;
  }

  char * makelist (const char * input, int type) {
    char * text = malloc (strlen(input) - 1);
    strncpy (text, &input[1], strlen (input) - 1);
    text[strlen (input) - 2] = '\0';
    char * s = strtok (text, " \t\v\n\f,");
    int listtype = 100;
    while (s) {
      char * dot = strchr (s, '.');
      var_t * var = varlookup (s, strlen(s) - (dot ? strlen(dot) : 0));
      if (!var) { // not a scalar or vector list, give up
	free (text);
	return NULL;
      }
      int vtype = var->type;
      while (dot) {
	vtype--;
	dot = strchr (dot+1, '.');
      }
      if (vtype < listtype)
	listtype = vtype;
      s = strtok (NULL, " \t\v\n\f,)");
    }
    if (listtype < type) { // incompatible list types
      free (text);
      return NULL;
    }      
    if (type >= 0)
      listtype = type;
    char * list = calloc (strlen("((scalar []){") + 1, sizeof (char));
    switch (listtype) {
    case scalar: strcpy (list, "((scalar []){"); break;
    case vector: strcpy (list, "((vector []){"); break;
    case tensor: strcpy (list, "((tensor []){"); break;
    default: assert (0); break;
    }
    strncpy (text, &input[1], strlen (input) - 1);
    text[strlen (input) - 2] = '\0';
    s = strtok (text, " \t\v\n\f,");
    while (s) {
      char * dot = strchr (s, '.');
      var_t * var = varlookup (s, strlen(s) - (dot ? strlen(dot) : 0));
      int vtype = var->type;
      while (dot) {
	vtype--;
	dot = strchr (dot+1, '.');
      }
      char member[80] = "";
      if (scope > 0 || var->i[0] < 0) { // dynamic allocation
	switch (listtype) {
	case scalar: {
	  switch (vtype) {
	  case 0: sprintf (member, "%s,", s); break;
	  case 1: sprintf (member, "%s.x,%s.y,", s, s); break;
	  case 2: sprintf (member, "%s.x.x,%s.x.y,%s.y.x,%s.y.y,",
			   s, s, s, s); break;
	  default: assert (0);
	  }
	  break;
	}
	case vector: {
	  switch (vtype) {
	  case 1: sprintf (member, "{%s.x,%s.y},", s, s); break;
	  case 2: sprintf (member, "{%s.x.x,%s.x.y},{%s.y.x,%s.y.y},",
			   s, s, s, s); break;
	  default: assert (0);
	  }
	  break;	  
	}
	case tensor: {
	  switch (vtype) {
	  case 2: sprintf (member, "{{%s.x.x,%s.x.y},{%s.y.x,%s.y.y}},",
			   s, s, s, s); break;
	  default: assert (0);
	  }
	  break;	  
	}
	default: assert(0);
	}
      }
      else { // static allocation
	int n = vtype - listtype;
	char coord[20];
	for (int i = 0; i < (1 << n); i++) {
	  switch (listtype) {
	  case scalar:
	    sprintf (coord, "%d,", var->i[i]); break;
	  case vector:
	    sprintf (coord, "{%d,%d},", var->i[2*i], var->i[2*i+1]); break;
	  case tensor:
	    sprintf (coord, "{{%d,%d},{%d,%d}},",
		     var->i[4*i], var->i[4*i+1], var->i[4*i+2], var->i[4*i+3]);
	    break;
	  default: assert (0);
	  }
	  strcat (member, coord);
	}
      }
      list = realloc (list, (strlen(list) + strlen(member) + 1)*sizeof(char));
      strcat (list, member);
      s = strtok (NULL, " \t\v\n\f,)");
    }
    free (text);
    char end[20];
    switch (listtype) {
    case scalar: strcpy (end, "-1})"); break;
    case vector: strcpy (end, "{-1,-1}})"); break;
    case tensor: strcpy (end, "{{-1,-1},{-1,-1}}})"); break;
    default: assert (0);
    }
    list = realloc (list, (strlen(list) + strlen(end) + 1)*sizeof(char));
    strcat (list, end);
    return list;
  }

  void new_field (var_t * var) {
    if (scope > 0)
      fprintf (yyout, "= new_%s(\"%s\")",
	       var->type == scalar ? "scalar" : 
	       var->type == vector ? "vector" : 
	       var->type == tensor ? "tensor" :
	       "internal_error",
	       var->v);
    else if (var->type == scalar) {
      fprintf (yyout, "= %d", nvar);
      var->i[0] = nvar;
      nscalars++;
      scalars = realloc (scalars, sizeof (Scalar)*nscalars);
      scalars[nscalars-1].i = nvar++;
      scalars[nscalars-1].name = strdup (var->v);
    }
    else if (var->type == vector) {
      fprintf (yyout, "= {%d,%d}", nvar, nvar + 1);    
      for (int i = 0; i < 2; i++)
	var->i[i] = nvar + i;
      nvectors++;
      vectors = realloc (vectors, sizeof (Vector)*nvectors);
      vectors[nvectors-1] = (Vector){nvar,nvar+1};
      vectors[nvectors-1].name = strdup (var->v);
      nvar += 2;
    }
    else if (var->type == tensor) {
      fprintf (yyout, "= {{%d,%d},{%d,%d}}",
	       nvar, nvar + 1, nvar + 2, nvar + 3);
      for (int i = 0; i < 4; i++)
	var->i[i] = nvar + i;
      ntensors++;
      tensors = realloc (tensors, sizeof (Tensor)*ntensors);
      tensors[ntensors-1] = (Tensor){{nvar,nvar+1},{nvar+2, nvar+3}};
      tensors[ntensors-1].name = strdup (var->v);
      nvar += 4;
    }
    else
      assert (0);
  }

  void declaration (char * var) {
    if (!strcmp (&var[strlen(var)-2], "[]")) {
      // automatic
      var[strlen(var)-2] = '\0';
      varpush (var, vartype, scope);
      var_t * v = varlookup (var, strlen(var));
      v->automatic = 1;
      fputs (yytext, yyout);
      new_field (v);
    }
    else {
      varpush (var, vartype, scope);
      fputs (yytext, yyout);
    }
    if (debug)
      fprintf (stderr, "%s:%d: declaration: %s\n", fname, line, var);
  }
  
#define nonspace(s) { while (strchr(" \t\v\n\f", *s)) s++; }
#define space(s) { while (!strchr(" \t\v\n\f", *s)) s++; }

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
SCALAR [a-zA-Z_0-9]+[.xyz]*

%%

"("                ECHO; para++;

")" {
  para--; if (para == invalpara) inval = 0;
  if (para < 0)
    return yyerror ("mismatched ')'");
  if (inforeach == 1 && scope == foreachscope && para == foreachpara) {
    ECHO;
    fclose (yyout);
    yyout = foreachfp;
    if (nreduct > 0) {
      fputs ("\n#undef _OMPSTART\n#define _OMPSTART ", yyout);
      for (int i = 0; i < nreduct; i++)
	fprintf (yyout, "double _%s = %s; ", reductvar[i], reductvar[i]);
      fputs ("\n#undef _OMPEND\n#define _OMPEND ", yyout);
      for (int i = 0; i < nreduct; i++)
	fprintf (yyout, "OMP(omp critical) if (_%s %s %s) %s = _%s; ",
		 reductvar[i], strcmp(reduction[i], "min") ? ">" : "<",
		 reductvar[i], reductvar[i], reductvar[i]);
      fprintf (yyout, "\n#line %d\n", line);
    }
    if (inforeach_face) {
      yyout = dopen ("foreach_face.h", "w");
      foreach_face_line = line;
      foreach_face_xy = face_xy;
      if (nreduct == 0) { // foreach_face (x)
	FILE * fp = dopen ("foreach.h", "r");
	int c;
	while ((c = fgetc (fp)) != EOF)
	  if (c == 'x')
	    foreach_face_xy = face_x;
	  else if (c == 'y')
	    foreach_face_xy = face_y;
	fclose (fp);
      }
    }
    else {
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
	     "static int %s (int i, double t) { ",
	     eventid[nevents]);
    assert (nevents < EVMAX);
    nexpr[nevents] = inevent;
    inevent = 4;
  }
  else if (inarg == para + 1) {
    fputs ("})", yyout);
    inarg = 0;
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
  scope--;
  if (scope < 0)
    return yyerror ("mismatched '}'");
  varpop();
  ECHO;
  if (foreachdim && scope == foreachdim)
    endforeachdim ();
  if (infunction && scope <= functionscope) {
    infunction = 0;
    if (debug)
      fprintf (stderr, "%s:%d: outfunction\n", fname, line);
  }
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
  int insthg = 0;
  if (foreachdim && scope == foreachdim && para == foreachdimpara) {
    ECHO; insthg = 1;
    endforeachdim ();
  }
  if (infunction && scope == functionscope) {
    if (scope > 0 && !infunctiondecl) {
      fputs ("; ", yyout); insthg = 1;
      infunction_declarations();
    }
    else if (!infunctiondecl) {
      infunction = 0;
      if (debug)
	fprintf (stderr, "%s:%d: outfunction\n", fname, line);
    }
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
	     "static int %s_expr%d (int * ip, double * tp) { "
	     "  int i = *ip; double t = *tp; "
	     "  int ret = (", eventid[nevents], inevent++);
  }
  else if (inevent == 4 && scope == eventscope && para == eventpara - 1) {
    ECHO;
    endevent ();
  }
  else if (inreturn) {
    fputs ("; }", yyout);
    inreturn = 0;
  }
  else if (!insthg)
    ECHO;
  invardecl = 0;
}

{ID}+{WS}*[=]{WS}*new{WS}+(scalar|vector|tensor) |
[=]{WS}*new{WS}+(scalar|vector|tensor) {
  char * type = strchr (yytext, '=');
  type = strstr (type, "new"); space(type); nonspace(type); 
  int newvartype = (!strcmp (type, "scalar") ? scalar :
		    !strcmp (type, "vector") ? vector :
		    !strcmp (type, "tensor") ? tensor :
		    -1);
  var_t * var;
  if (yytext[0] == '=') {
    var = vartop();
    assert (var && var->scope == scope);
    if (var->type != newvartype)
      return yyerror ("type mismatch in new");
  }
  else {
    char * s = yytext; space (s); *s = '\0';
    var = varlookup (yytext, strlen(yytext));
    if (!var) {
      fprintf (stderr, "%s:%d: undeclared %s '%s'", fname, line, type, yytext);
      return 1;
    }
    if (var->type != newvartype)
      return yyerror ("type mismatch in new");
    fprintf (yyout, "%s ", yytext);
  }
  new_field (var);
  if (debug)
    fprintf (stderr, "%s:%d: new %s: %s\n", fname, line, type, var->v);
}

(scalar|vector|tensor){WS}+[a-zA-Z0-9_\[\]]+ {
  char * var = strstr(yytext,"scalar");
  vartype = scalar;
  if (!var) {
    var = strstr(yytext,"vector");
    vartype = vector;
    if (!var) {
      var = strstr(yytext,"tensor");
      vartype = tensor;
    }
  }
  var = &var[7];
  nonspace (var);
  if (*var == '[') {
    // scalar [..
    for (; *var != '\0'; var++)
      if (*var == '[') brack++;
      else if (*var == ']') brack--;
    ECHO;
  }
  else if (para == 0) { /* declaration */
    declaration (var);
    invardecl = scope + 1;
  }
  else if (para == 1) { /* function prototype (no nested functions) */
    ECHO;
    if (debug)
      fprintf (stderr, "%s:%d: proto: %s\n", fname, line, var);
    varpush (var, vartype, scope + 1);
    invardecl = 0;
  }
  else
    ECHO;
}

,{WS}*[a-zA-Z0-9_\[\]]+ {
  if (invardecl == scope + 1) {
    char * var = &yytext[1];
    nonspace (var);
    declaration (var);
  }
  else
    REJECT;
}

return{WS} {
  // returning from a function: delete automatic fields before returning
  // note that this assumes that the function scope is always 1 
  // (i.e. no nested functions allowed).
  inreturn = 1;
  fputs ("{ ", yyout);
  delete_automatic (0);
  ECHO;
}

val{WS}*[(]    {
  inval = 1; invalpara = para++;
  ECHO;
}

{SCALAR}{WS}*\[{WS}*. {
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
	fprintf (yyout, "val(%s,0,0)", boundaryrep(yytext));
      else {
	fprintf (yyout, "val(%s", boundaryrep(yytext));
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

{SCALAR}{WS}*\[{WS}*(left|right|top|bottom){WS}*\]{WS}*= {
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
    char * func = malloc (strlen (yytext) + 1);
    strcpy (func, yytext);
    char * s1 = func;
    while (*s1 != '\0') {
      if (*s1 == '.')
	*s1 = '_';
      s1++;
    }
    boundaryfunc = malloc (strlen ("_method[].boundary[] = _;") + 
			   2*strlen (b) + 
			   strlen(yytext) + 
			   strlen (func) + 1);
    sprintf (boundaryfunc, "_method[%s].boundary[%s] = _%s%s;", 
	     yytext, b, func, b);
    fprintf (yyout, 
	     "double _%s%s (Point point, scalar _s) {"
	     " int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED (jg);"
	     " POINT_VARIABLES; "
	     " return ", 
	     func, b);
    free (func);
    inboundary = inforeach_boundary = 1;
    strcpy (boundaryvar, yytext);
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

Point{WS}+point[^{ID}] {
  /* Point point */
  if (debug)
    fprintf (stderr, "%s:%d: infunction\n", fname, line);
  infunction = 1; infunctiondecl = 0;
  functionscope = scope; functionpara = para;
  if (yytext[yyleng-1] == ')') para--;
  ECHO;
}

for{WS}*[(]{WS}*(scalar|vector){WS}+{ID}+{WS}+in{WS}+ {
  /* for (scalar .. in .. */
  char * s = strchr (&yytext[3], 'r'); s++;
  int vartype = s[-2] == 'o' ? vector : scalar;
  nonspace (s);
  char * id = s;
  while (!strchr (" \t\v\n\f", *s)) s++;
  *s = '\0';
  char * list = malloc (sizeof(char)); list[0] = '\0';
  int nl = 1, c;
  char last[] = ")";
  while ((c = input()) != EOF) {
    if (c == ')')
      break;
    else {
      list = realloc (list, (nl + 1)*sizeof(char));
      list[nl - 1] = c;
      list[nl++] = '\0';
    }
  }
  if (c != ')')
    return yyerror ("expecting ')'");
  static int i = 0;
  if (list[0] == '{') {
    char * oldlist = list;
    list = makelist (oldlist, vartype);
    free (oldlist);
    if (list == NULL)
      return yyerror ("invalid list");
  }
  else
    fprintf (yyout, "if (%s) ", list);
  if (debug)
    fprintf (stderr, "%s:%d: %s\n", fname, line, list);
  fprintf (yyout,
	   "for (%s %s = *%s, *_i%d = %s; %s%s >= 0; %s = *++_i%d%s",
	   vartype == scalar ? "scalar" : "vector", id,
	   list, i, list, 
	   id, vartype == scalar ? "" : ".x", 
	   id, i, last);
  free (list);
  i++;
  varpush (id, vartype, scope);
}

for{WS}*[(][^)]+,[^)]+{WS}+in{WS}+[^)]+,[^)]+[)] {
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
  static int index = 0;
  char * sc = NULL, * vc = NULL;
  for (int i = 0; i < nid; i++) {
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
  for (int i = 1; i < nid; i++)
    fprintf (yyout, ", %s = *%s", id[i], list[i]);
  if (sc)
    fprintf (yyout, "; %s >= 0; ", sc);
  else
    fprintf (yyout, "; %s.x >= 0; ", vc);
  fprintf (yyout, "%s = *++_i%d", id[0], index);
  for (int i = 1; i < nid; i++)
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

event{WS}+{ID}+{WS}*[(] {
  /* event (... */
  char * id = yytext;
  space (id); nonspace (id);
  char * s = id; while (!strchr(" \t\v\n\f(", *s)) s++; *s = '\0';
  fprintf (yyout, 
	   "static int %s_expr%d (int * ip, double * tp) {"
	   "  int i = *ip; double t = *tp;"
	   "  int ret = (", id, inevent++);
  eventscope = scope; eventpara = ++para;
  eventarray[nevents] = 0;
  eventfile[nevents] = strdup (fname);
  eventid[nevents] = strdup (id);
  eventline[nevents] = line;
}

[it]{WS}*={WS}*[{][^}]*[}] {
  if (inevent == 1) {
    eventarray[nevents] = yytext[0];
    yytext[yyleng-1] = '\0';
    eventarray_elems[nevents] = strdup (strchr (yytext, '{'));
    fputs ("1); "
	   "  *ip = i; *tp = t; "
	   "  return ret; "
	   "} ", yyout);
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

reduction{WS}*[(](min|max):{ID}+[)] {
  char * s = strchr (yytext, '('), * s1 = strchr (yytext, ':');
  *s1 = '\0'; s1++;
  assert (nreduct < REDUCTMAX);
  strcpy (reduction[nreduct], ++s);
  yytext[yyleng-1] = '\0';
  strcpy (reductvar[nreduct++], s1);
}

[{]({SP}*[a-zA-Z_0-9.]+{SP}*,{SP}*)*[a-zA-Z_0-9.]+{SP}*[}] {
  // {h, zb, ...}
  char * list = makelist (yytext, -1);
  if (list == NULL)
    ECHO;
  else {
    if (debug)
      fprintf (stderr, "%s:%d: %s\n", fname, line, list);
    fputs (list, yyout);
    free (list);
  }
}

{ID}+ {
  if (inforeach) {
    for (int i = 0; i < nreduct; i++)
      if (!strcmp (yytext, reductvar[i])) {
	fputc ('_', yyout);
	break;
      }
  }
  ECHO;
}

{SCALAR}[.][a-wA-Z_0-9]+ {
  // scalar methods
  char * dot1 = strchr (yytext, '.');
  var_t * var = varlookup (yytext, strlen(yytext) - strlen(dot1));
  if (var) {
    char * dot = dot1;
    while (dot1) {
      dot = dot1;
      dot1 = strchr (dot1 + 1, '.');
    }
    dot[0] = '\0'; dot++;
    fprintf (yyout, "_method[%s].%s", yytext, dot);
    if (debug)
      fprintf (stderr, "%s:%d: _method[%s].%s\n", fname, line, yytext, dot);
  }
  else
    ECHO;
}

{ID}+{WS}+{ID}+{WS}*[(]{WS}*struct{WS}+{ID}+{WS}+{ID}+{WS}*[)] {
  // function declaration with struct argument
  ECHO;
  char * s = yytext;
  space (s); nonspace (s);
  args = realloc (args, sizeof (char *)*++nargs);
  argss = realloc (argss, sizeof (char *)*nargs);
  char * s1 = s; while (!strchr(" \t\v\n\f(", *s1)) s1++;
  *s1++ = '\0'; s1 = strstr (s1, "struct"); space (s1); nonspace (s1);
  char * s2 = s1; space (s2); *s2 = '\0';
  args[nargs-1] = strdup (s);
  argss[nargs-1] = strdup (s1);
  if (debug)
    fprintf (stderr, "%s:%d: function '%s' with struct '%s'\n", 
	     fname, line, s, s1);
}

{ID}+{WS}*[(]{WS}*{ID}+{WS}*[)] {
  // function call without 'args' assignment
  char * s = yytext; space (s);
  int len = s - yytext;
  for (int i = 0; i < nargs && !inarg; i++)
    if (strlen(args[i]) == len && !strncmp (args[i], yytext, len)) {
      ECHO;
      inarg = 1;
    }
  if (!inarg)
    REJECT;
  inarg = 0;
}

{ID}+{WS}*[(] {
  // function call with 'args' assignment
  char * s = yytext; space (s);
  int len = s - yytext;
  for (int i = 0; i < nargs && !inarg; i++)
    if (strlen(args[i]) == len && !strncmp (args[i], yytext, len)) {
      ECHO; para++;
      inarg = para;
      fprintf (yyout, "(struct %s){", argss[i]);
    }
  if (!inarg)
    REJECT;
}

{ID}+{WS}*=[^=] {
  // arguments of function call with 'args' assignment
  if (inarg && para == inarg) {
    fputc('.', yyout);
    ECHO;
  }
  else
    REJECT;
}

"# "[0-9]+" "({SP}?\"([^\"\\\n]|{ES})*\")+ {
  // line numbers
  ECHO;
  char * ln = yytext;
  space (ln); nonspace(ln);
  char * name = ln; space (name); *name++ = '\0'; nonspace(name);
  line = atoi(ln) - 1;
  free (fname);
  name[strlen(name)-1] = '\0'; name++;
  fname = strdup (name);
}

^{SP}*@{SP}*def{SP}+ {
  // @def ... @
  fputs ("#define ", yyout);
  register int c;
  while ((c = input()) != EOF && c != '@') {
    if (c == '\n')
      fputc (' ', yyout);
    else
      fputc (c, yyout);
  }
  fprintf (yyout, "\n#line %d\n", line);
}

^{SP}*@{SP}*{ID}+ {
  yytext = strchr(yytext, '@'); yytext++;
  fprintf (yyout, "#%s", yytext);
  register int oldc = 0, c;
  for (;;) {
    while ((c = getput()) != '\n' && c != EOF)
      oldc = c;    /* eat up text of preproc */
    if (c == EOF || oldc != '\\')
      break;
  }
}

^{SP}*@.*" Pragma(" {
  yytext = strchr(yytext, '@'); yytext++;
  char * s = strstr (yytext, "Pragma("); *s++ = '\0';
  fprintf (yyout, "#%s", yytext);
  fputs ("_Pragma(", yyout);
  register int oldc = 0, c;
  for (;;) {
    while ((c = getput()) != '\n' && c != EOF)
      oldc = c;    /* eat up text of preproc */
    if (c == EOF || oldc != '\\')
      break;
  }  
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
	      char ** grid, int * default_grid,
	      const char * dir);

int endfor (FILE * fin, FILE * fout)
{
  yyin = fin;
  yyout = fout;
  line = 0, scope = para = 0;
  inforeach = foreachscope = foreachpara = 
    inforeach_boundary = inforeach_face = 0;
  invardecl = 0;
  inval = invalpara = 0;
  brack = inarray = 0;
  inevent = inreturn = 0;
  foreachdim = 0;
  inboundary = 0;
  infunction = 0;
  inarg = 0;
  int ret = yylex();
  if (!ret) {
    if (scope > 0)
      ret = yyerror ("mismatched '{'");
    else if (para > 0)
      ret = yyerror ("mismatched '('");
  }
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

void cleanup (int status, const char * dir)
{
  if (!debug && dir) {
    char command[80] = "rm -r -f ";
    strcat (command, dir);
    int s = system (command); s = s;
  }
  exit (status);
}

FILE * dopen (const char * fname, const char * mode)
{
  char * out = malloc (strlen (dir) + strlen (fname) + 2);
  strcpy (out, dir); strcat (out, "/"); strcat (out, fname);
  FILE * fout = fopen (out, mode);
  free (out);
  return fout;
}

void compdir (FILE * fin, FILE * fout, char * grid)
{
  if (endfor (fin, fout))
    cleanup (1, dir);
  fclose (fout);

  fout = dopen ("grid.h", "w");
  /* new variables */
  fprintf (fout,
	   "int datasize = %d*sizeof (double);\n",
	   nvar);
  /* events */
  for (int i = 0; i < nevents; i++) {
    char * id = eventid[i];
    fprintf (fout, "static int %s (int i, double t);\n", id);
    for (int j = 0; j < nexpr[i]; j++)
      fprintf (fout,
	       "static int %s_expr%d (int * ip, double * tp);\n",
	       id, j);
    if (eventarray[i])
      fprintf (fout, "static %s %s_array[] = %s,-1};\n", 
	       eventarray[i] == 'i' ? "int" : "double", id,
	       eventarray_elems[i]);
  }
  fputs ("Event Events[] = {\n", fout);
  for (int i = 0; i < nevents; i++) {
    char * id = eventid[i];
    fprintf (fout, "  { 0, %d, %s, {", nexpr[i], id);
    int j;
    for (j = 0; j < nexpr[i] - 1; j++)
      fprintf (fout, "%s_expr%d, ", id, j);
    fprintf (fout, "%s_expr%d}, ", id, j);
    if (eventarray[i] == 'i')
      fprintf (fout, "%s_array, ", id);
    else
      fprintf (fout, "((void *)0), ");
    if (eventarray[i] == 't')
      fprintf (fout, "%s_array,\n", id);
    else
      fprintf (fout, "((void *)0),\n");
    fprintf (fout, "    \"%s\", %d},\n", eventfile[i], eventline[i]);
  }
  fputs ("  { 1 }\n};\n", fout);
  /* boundaries */
  for (int i = 0; i < nboundary; i++)
    fprintf (fout, "static void _boundary%d (void);\n", i);
  /* methods */
  fprintf (fout, "void %s_methods(void);\n", grid);
  fputs ("static void init_solver (void) {\n", fout);
  /* scalar methods */
  fputs ("  _method = calloc (datasize/sizeof(double), sizeof (Methods));\n", 
	 fout);
  /* list of all scalars */
  fprintf (fout, 
	   "  all = malloc (sizeof (scalar)*%d);\n"
	   "  for (int i = 0; i < %d; i++)\n"
	   "    all[i] = i;\n"
	   "  all[%d] = -1;\n", nvar + 1, nvar, nvar);
  if (fpe)
    fputs ("  set_fpe();\n", fout);
  if (catch)
    fputs ("  catch_fpe();\n", fout);
  fprintf (fout, "  %s_methods();\n", grid);
  for (int i = 0; i < nscalars; i++)
    fprintf (fout, "  init_scalar (%d, \"%s\");\n", 
	     scalars[i].i, scalars[i].name);
  for (int i = 0; i < nvectors; i++)
    fprintf (fout, "  init_vector ((vector){%d,%d}, \"%s\");\n", 
	     vectors[i].x, vectors[i].y, vectors[i].name);
  for (int i = 0; i < ntensors; i++)
    fprintf (fout, "  init_tensor ((tensor){{%d,%d},{%d,%d}}, \"%s\");\n", 
	     tensors[i].x.x, tensors[i].x.y,
	     tensors[i].y.x, tensors[i].y.y, tensors[i].name);
  for (int i = 0; i < nboundary; i++)
    fprintf (fout, "  _boundary%d();\n", i);
  fputs ("  init_events();\n}\n", fout);
  fclose (fout);
}

int main (int argc, char ** argv)
{
  char * cc = getenv ("CC99"), command[1000], command1[1000] = "";
  if (cc == NULL)
    strcpy (command, CC99);
  else
    strcpy (command, cc);
  strcat (command, " -I");
  strcat (command, LIBDIR);
  char * file = NULL;
  int i, dep = 0, tags = 0;
  for (i = 1; i < argc; i++) {
    if (!strncmp (argv[i], "-grid=", 6))
      ;
    else if (!strcmp (argv[i], "-MD"))
      dep = 1;
    else if (!strcmp (argv[i], "-tags"))
      tags = 1;
    else if (!strcmp (argv[i], "-debug"))
      debug = 1;
    else if (!strcmp (argv[i], "-fpe"))
      fpe = 1;
    else if (!strcmp (argv[i], "-catch"))
      catch = 1;
    else if (!strcmp (argv[i], "-o")) {
      strcat (command1, " ");
      strcat (command1, argv[i++]);
      if (i < argc) {
	strcat (command1, " ");
	strcat (command1, argv[i]);
      }
    }
    else if (argv[i][0] != '-' && 
	     (tags || !strcmp (&argv[i][strlen(argv[i]) - 2], ".c"))) {
      if (file) {
	fprintf (stderr, "usage: qcc -grid=[GRID] [OPTIONS] FILE.c\n");
	return 1;
      }
      file = argv[i];
    }
    else if (!file) { 
      strcat (command, " ");
      strcat (command, argv[i]);
    }
    else {
      strcat (command1, " ");
      strcat (command1, argv[i]);
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
    int default_grid;
    includes (argc, argv, out, &grid, &default_grid, dep || tags ? NULL : dir);
    if (!dep && !tags) {
      char * basename = strdup (file), * ext = basename;
      while (*ext != '\0' && *ext != '.') ext++;
      char * cpp = malloc (strlen(basename) + strlen("-cpp") + strlen(ext) + 1);
      if (*ext == '.') {
	*ext = '\0';
	strcpy (cpp, basename);
	strcat (cpp, "-cpp");
	*ext = '.';
	strcat (cpp, ext);
      }
      else {
	strcpy (cpp, basename);
	strcat (cpp, "-cpp");
      }
      FILE * fin = dopen (file, "r");
      if (!fin) {
	perror (file);
	cleanup (1, dir);
      }
      FILE * fout = dopen (cpp, "w");
      if (fpe)
	fputs ("#include <string.h>\n"
	       "#include <fenv.h>\n", 
	       fout);
      if (catch)
	fputs ("#define _CATCH last_point = point;\n", fout);
      else
	fputs ("#define _CATCH\n", fout);
      fputs ("#include \"common.h\"\n", fout);
      /* catch */
      if (catch)
	fputs ("void catch_fpe (void);\n", fout);
      /* undefined value */
      if (fpe)
	/* Initialises unused memory with "signaling NaNs".  
	 * This is probably not very portable, tested with
	 * gcc (Debian 4.4.5-8) 4.4.5 on Linux 2.6.32-5-amd64.
	 * This blog was useful:
	 *   http://codingcastles.blogspot.co.nz/2008/12/nans-in-c.html 
	 */
	fputs ("double undefined;\n"
	       "static void set_fpe (void) {\n"
	       "  long lnan = 0x7ff0000000000001;\n"
	       "  assert (sizeof (long) == sizeof (double));\n"
	       "  memcpy (&undefined, &lnan, sizeof (long));\n"
	       "  feenableexcept (FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);\n"
	       "}\n",
	       fout);
      else
	fputs ("#define undefined DBL_MAX\n", fout);
      fputs ("@include \"grid.h\"\n", fout);
      /* grid */
      if (default_grid)
	fprintf (fout, "#include \"grid/%s.h\"\n", grid);
      int c;
      while ((c = fgetc (fin)) != EOF)
	fputc (c, fout);
      fclose (fout);
      fclose (fin);
      fout = dopen (file, "w");
      if (!fout) {
	perror (file);
	cleanup (1, dir);
      }

      char preproc[1000];
      strcpy (preproc, "cd ");
      strcat (preproc, dir);
      strcat (preproc, " && ");
      strcat (preproc, command);
      strcat (preproc, " -E ");
      strcat (preproc, cpp);

      fin = popen (preproc, "r");
      if (!fin) {
	perror (preproc);
	cleanup (1, dir);
      }

      compdir (fin, fout, grid);
      int status = pclose (fin);
      if (status == -1 ||
	  (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				    WTERMSIG (status) == SIGQUIT)))
	cleanup (1, dir);

      strcat (command, " ");
      strcat (command, dir);
      strcat (command, "/");
      strcat (command, file); 
     strcat (command, command1);
    }
  }
  else if (dep || tags) {
    fprintf (stderr, "usage: qcc -grid=[GRID] [OPTIONS] FILE.c\n");
    cleanup (1, dir);
  }
  /* compilation */
  if (!dep && !tags) {
    if (debug)
      fprintf (stderr, "command: %s\n", command);
    status = system (command);
    if (status == -1 ||
	(WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				  WTERMSIG (status) == SIGQUIT)))
      cleanup (1, dir);
    cleanup (WEXITSTATUS (status), dir);
  }
  cleanup (0, dir);
}
