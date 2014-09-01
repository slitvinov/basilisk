%option noyywrap
%{
  #include <unistd.h>
  #include <string.h>
  #include <ctype.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <sys/wait.h>
  #include <assert.h>

  enum { scalar, vector, tensor, struct_type = -1 };

  typedef struct { int i; char * name; } Scalar;
  typedef struct { int x, y, face; char * name; } Vector;
  typedef struct { Vector x, y; char * name; } Tensor;

  int debug = 0, catch = 0, nolineno = 0, events = 0;
  char dir[] = ".qccXXXXXX";

  int nvar = 0, nconst = 0, nevents = 0;
  int line;
  int scope, para, inforeach, foreachscope, foreachpara, 
    inforeach_boundary, inforeach_face, nmaybeconst = 0;
  int invardecl, vartype, varsymmetric, varface, varvertex, varmaybeconst;
  char * varconst;
  int inval, invalpara;
  int brack, inarray;
  int inreturn;

  #define EVMAX 100
  int inevent, eventscope, eventpara;
  char eventarray[EVMAX], * eventfile[EVMAX], * eventid[EVMAX];
  char * eventfunc[EVMAX];
  char * eventarray_elems[EVMAX];
  int nexpr[EVMAX], eventline[EVMAX], eventparent[EVMAX], eventchild[EVMAX];
  int eventlast[EVMAX];

  int foreachdim, foreachdimpara, foreachdimline;
  FILE * foreachdimfp;

  int foreach_child, foreach_child_scope, foreach_child_para;

  int inattr, attrscope, attrline;
  FILE * attrfp;

  #define REDUCTMAX 10
  char foreachs[80], * fname;
  FILE * foreachfp;
  char reduction[REDUCTMAX][4], reductvar[REDUCTMAX][80];
  int nreduct;

  int inboundary;
  char boundaryvar[80], boundarydir[80];
  FILE * boundaryfp = NULL, * boundaryheader = NULL;
  int boundarycomponent, nboundary = 0, nsetboundary = 0;
  int boundaryindex[80];

  int infunction, infunctiondecl, functionscope, functionpara, inmain;
  char * return_type = NULL;
  FILE * dopen (const char * fname, const char * mode);

  int infunctionproto;

  int foreach_line;
  enum { face_x, face_y, face_xy };
  int foreach_face_xy;

  char ** args = NULL, ** argss = NULL;
  int nargs = 0, inarg;

  typedef struct { 
    char * v, * constant;
    int type, args, scope, automatic, symmetric, face, vertex, maybeconst;
    int i[4];
    char * conditional;
  } var_t;
  var_t _varstack[100]; int varstack = -1;
  var_t * varpush (const char * s, int type, int scope, int maybeconst) {
    var_t * v = NULL;
    if (s[0] != '\0') {
      char * f = malloc (strlen (s) + 1);
      strcpy (f, s);
      char * q = f;
      int na = 0;
      while ((q = strchr (q, '['))) {
	*q++ = '\0'; na++;
      }
      _varstack[++varstack] = (var_t) { f, NULL, type, na, scope, 
					0, 0, 0, 0, maybeconst, {-1} };
      v = &(_varstack[varstack]);
    }
    return v;
  }
  var_t ** foreachconst = NULL;

  char * makelist (const char * input, int type);

  char * automatic_list (int scope, int conditional) {
    char * list = NULL;
    int i;
    for (i = varstack; i >= 0 && _varstack[i].scope > scope; i--) {
      var_t var = _varstack[i];
      if (var.automatic && !var.constant && (conditional || !var.conditional)) {
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
    return list;
  }

  void delete_automatic (int scope) {
    char * list = automatic_list (scope, 0);
    if (list) {
      strcat (list, "}");
      char * slist = makelist (list, scalar);
      if (debug)
	fprintf (stderr, "%s:%d: deleting %s\n", fname, line, list);
      fprintf (yyout, " delete (%s); ", slist);
      free (slist);
      free (list);
    }
    int i;
    for (i = varstack; i >= 0 && _varstack[i].scope > scope; i--) {
      var_t var = _varstack[i];
      if (var.automatic && var.conditional) {
	if (debug)
	  fprintf (stderr, "%s:%d: deleting conditional %s\n", 
		   fname, line, var.v);
	char list[80];
	sprintf (list, "{%s}", var.v);
	char * slist = makelist (list, scalar);
	fprintf (yyout, " { if (!%s) delete (%s); } ", 
		 var.conditional, slist);
	free (slist);
      }
    }
  }

  void varpop () {
    delete_automatic (scope);
    while (varstack >= 0 && _varstack[varstack].scope > scope) {
      if (debug)
	fprintf (stderr, "%s:%d: '%s' out of scope\n", fname, line, 
		 _varstack[varstack].v);
      free (_varstack[varstack].v);
      free (_varstack[varstack].conditional);
      free (_varstack[varstack--].constant);
    }
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
    int i;
    for (i = varstack; i >= 0; i--)
      if (strlen(_varstack[i].v) == len && 
	  !strncmp(s, _varstack[i].v, len)) {
	if (_varstack[i].type < 0) // not a field
	  return NULL;
	if (s[len-1] == '.' && _varstack[i].type != vector)
	  return NULL;
	return &_varstack[i];
      }
    return NULL;
  }

  int get_vartype (char * s) {
    int len = strlen (s), i;
    for (i = varstack; i >= 0; i--)
      if (strlen(_varstack[i].v) == len && 
	  !strncmp(s, _varstack[i].v, len)) {
	if (s[len-1] == '.' && _varstack[i].type != vector)
	  return -100;
	return _varstack[i].type;
      }
    return -100;
  }

  int rotate (FILE * fin, FILE * fout, int n);

  void writefile (FILE * fp, char x, char y, int line1,
		  const char * condition) {
    if (condition)
      fprintf (yyout, " if (%s) {", condition);
    fprintf (yyout, "\n#line %d\n", line1);
    rotate (fp, yyout, x == 'y');
    if (condition)
      fputs (" } ", yyout);
  }

  void endforeachdim () {
    fclose (yyout);
    yyout = foreachdimfp;
    if (foreachdim > 1)
      fputc ('{', yyout);
    FILE * fp = dopen ("_dimension.h", "r");
    writefile (fp, 'y', 'x', foreachdimline, NULL);
    writefile (fp, 'x', 'y', foreachdimline, NULL);
    fclose (fp);
    if (foreachdim > 1)
      fputc ('}', yyout);
    foreachdim = 0;
  }

  void foreachbody() {
    if (inforeach_face) {
      // foreach_face()
      fputs ("foreach_face_generic()", yyout);
      FILE * fp = dopen ("_foreach_body.h", "r");
      if (foreach_face_xy == face_xy) {
	fputs (" { int jg = -1; VARIABLES; ", yyout);
	writefile (fp, 'y', 'x', foreach_line, "is_face_y()");
	fputs (" } { int ig = -1; VARIABLES; ", yyout);
	writefile (fp, 'x', 'y', foreach_line, "is_face_x()");
	fputs (" } ", yyout);
      }
      else if (foreach_face_xy == face_x) {
	fputs (" { int ig = -1; VARIABLES; ", yyout);
	writefile (fp, 'x', 'y', foreach_line, "is_face_x()");
	fputs (" } ", yyout);
      }
      else {
	fputs (" { int jg = -1; VARIABLES; ", yyout);
	writefile (fp, 'x', 'y', foreach_line, "is_face_y()");
	fputs (" } ", yyout);
      }
      fputs (" end_foreach_face_generic()\n", yyout);
      fprintf (yyout, "#line %d\n", line);
      fclose (fp);
    }
    else {
      // foreach()
      FILE * fp = dopen ("_foreach.h", "r");
      int c;
      while ((c = fgetc (fp)) != EOF)
	fputc (c, yyout);
      fclose (fp);
      fp = dopen ("_foreach_body.h", "r");
      while ((c = fgetc (fp)) != EOF)
	fputc (c, yyout);
      fclose (fp);
    }
    fprintf (yyout, " end_%s();", foreachs);
  }

  char * macro_from_var (const char * name) {
    char * macro = strdup (name), * s = macro;
    while (*s != '\0') {
      if (*s == '.')
	*s = '_';
      s++;
    }
    return macro;
  }

  // combinations for each constant field
  void maybeconst_combinations (int line, void (* body)(void)) {
    if (nmaybeconst) {
      int n = 1 << nmaybeconst, bits;
      for (bits = 0; bits < n; bits++) {
	fputs ("\nif (", yyout);
	int i;
	for (i = 0; i < nmaybeconst; i++) {
	  fprintf (yyout, "%sis_constant(%s%s)",
		   (bits & (1 << i)) ? "" : "!",
		   foreachconst[i]->v,
		   foreachconst[i]->type == vector ? ".x" : "");
	  if (i < nmaybeconst - 1)
	    fputs (" && ", yyout);
	}
	fputs (") {\n", yyout);
	for (i = 0; i < nmaybeconst; i++)
	  if (foreachconst[i]->type == scalar) {
	    if (bits & (1 << i))
	      fprintf (yyout,
		       "const double _const_%s = _constant[%s-_NVARMAX];\n"
		       "NOT_UNUSED(_const_%s);\n"
		       "#undef val_%s\n"
		       "#define val_%s(a,i,j) _const_%s\n",
		       foreachconst[i]->v, foreachconst[i]->v, 
		       foreachconst[i]->v, foreachconst[i]->v, 
		       foreachconst[i]->v, foreachconst[i]->v);
	    else
	      fprintf (yyout, 
		       "#undef val_%s\n"
		       "#define val_%s(a,i,j) val(a,i,j)\n",
		       foreachconst[i]->v, foreachconst[i]->v);
	  }
	  else if (foreachconst[i]->type == vector) {
	    if (bits & (1 << i))
	      fprintf (yyout, "const struct { double x, y; } _const_%s = "
		       "{_constant[%s.x -_NVARMAX],"
		       " _constant[%s.y -_NVARMAX]};\n"
		       "NOT_UNUSED(_const_%s);\n",
		       foreachconst[i]->v, foreachconst[i]->v, 
		       foreachconst[i]->v, foreachconst[i]->v);
	    int c;
	    for (c = 'x'; c <= 'y'; c++)
	      if (bits & (1 << i))
		fprintf (yyout,
			 "#undef val_%s_%c\n"
			 "#define val_%s_%c(a,i,j) _const_%s.%c\n",
			 foreachconst[i]->v, c,
			 foreachconst[i]->v, c, foreachconst[i]->v, c);
	      else
		fprintf (yyout, 
			 "#undef val_%s_%c\n"
			 "#define val_%s_%c(a,i,j) val(a,i,j)\n",
			 foreachconst[i]->v, c, foreachconst[i]->v, c);
	  }
	fprintf (yyout, "#line %d\n", line);
	body();
	fputs (" }", yyout);
      }
    }
    else
      body();
  }

  void endforeach () {
    fclose (yyout);
    yyout = foreachfp;
    maybeconst_combinations (foreach_line, foreachbody);
    if (nreduct > 0)
      fprintf (yyout,
	       "\n"
	       "#undef _OMPSTART\n"
	       "#undef _OMPEND\n"
	       "#define _OMPSTART\n"
	       "#define _OMPEND\n"
	       "#line %d\n",
	       line + 1);
    fputs (" }", yyout);
    inforeach = inforeach_boundary = inforeach_face = 0;
  }

  void endevent() {
    fprintf (yyout, "  return 0; } ");
    inevent = 0;
    nevents++;
  }

  void endattr() {
    inattr = 0;
    fclose (yyout); yyout = attrfp;
    FILE * fp = dopen ("_attribute.h", "r");
    FILE * out = dopen ("_attributes.h", "a");
    fprintf (out, "\n#line %d \"%s\"\n", attrline, fname);
    int c;
    while ((c = fgetc (fp)) != EOF) {
      if (c == '\n')
	fputc (c, yyout);
      fputc (c, out);
    }
    fclose (fp);
    fclose (out);
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

  void boundary_staggering (const char * dir, int component, FILE * fp) {
    if (!strcmp (dir, "left")) {
      if (component == 'y')
	fputs (" int ig = -1, jg = -1;", fp);
      else
	fputs (" int ig = -1, jg = 0;", fp);
    }
    else if (!strcmp (dir, "right")) {
      if (component == 'y')
	fputs (" int ig = 1, jg = -1;", fp);
      else
	fputs (" int ig = 1, jg = 0;", fp);
    }
    else if (!strcmp (dir, "top")) {
      if (component == 'x')
	fputs (" int ig = -1, jg = 1;", fp);
      else
	fputs (" int ig = 0, jg = 1;", fp);
    }
    else if (!strcmp (dir, "bottom")) {
      if (component == 'x')
	fputs (" int ig = -1, jg = -1;", fp);
      else
	fputs (" int ig = 0, jg = -1;", fp);      
    }
    else
      assert (0);
    fputs (" NOT_UNUSED(ig); NOT_UNUSED (jg);", fp);
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
	int n = vtype - listtype, i;
	char * constant = var->constant ? "_NVARMAX + " : "";
	char coord[80];
	for (i = 0; i < (1 << n); i++) {
	  switch (listtype) {
	  case scalar:
	    sprintf (coord, "%s%d,", constant, var->i[i]); break;
	  case vector:
	    sprintf (coord, "{%s%d,%s%d},", 
		     constant, var->i[2*i], constant, var->i[2*i+1]); break;
	  case tensor:
	    sprintf (coord, "{{%s%d,%s%d},{%s%d,%s%d}},",
		     constant, var->i[4*i], constant, var->i[4*i+1], 
		     constant, var->i[4*i+2], constant, var->i[4*i+3]);
	    break;
	  default: assert (0);
	  }
	  strcat (member, coord);
	}
      }
      list = realloc (list, (strlen(list) + strlen(member) + 1)*sizeof(char));
      list = strcat (list, member);
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
    if (scope > 0) {
      // automatic variables
      fprintf (yyout, " new_%s%s%s(\"%s\"",
	       var->constant ? "const_" : "",
	       var->symmetric ? "symmetric_" : 
	       var->face ? "face_" : 
	       var->vertex ? "vertex_" : 
	       "",
	       var->type == scalar ? "scalar" : 
	       var->type == vector ? "vector" : 
	       var->type == tensor ? "tensor" :
	       "internal_error",
	       var->v);
      if (var->constant) {
	if (var->type == scalar) {
	  fprintf (yyout, ", %d, %s", nconst, var->constant);
	  nconst++;
	}
	else {
	  fprintf (yyout, ", %d, (double [])%s", 
		   nconst, var->constant);
	  nconst += 2;
	}
      }
      fputc (')', yyout);
    }
    else {
      // global constants
      if (var->constant) {
	if (var->type == scalar) {
	  fprintf (yyout, " _NVARMAX + %d", nconst);
	  var->i[0] = nconst++;
	}
	else if (var->type == vector) {
	  fprintf (yyout, " {_NVARMAX + %d,_NVARMAX + %d}", 
		   nconst, nconst + 1);
	  int i;
	  for (i = 0; i < 2; i++)
	    var->i[i] = nconst++;
	}
	else
	  assert (0);
      }
      // global variables
      else if (var->type == scalar) {
	fprintf (yyout, " %d", nvar);
	var->i[0] = nvar++;
      }
      else if (var->type == vector) {
	fprintf (yyout, " {%d,%d}", nvar, nvar + 1);    
	int i;
	for (i = 0; i < 2; i++)
	  var->i[i] = nvar++;
      }
      else if (var->type == tensor) {
	fprintf (yyout, " {{%d,%d},{%d,%d}}",
		 nvar, nvar + 1, nvar + 2, nvar + 3);
	int i;
	for (i = 0; i < 4; i++)
	  var->i[i] = nvar++;
      }
      else
	assert (0);
    }
    if (!var->constant && var->maybeconst) {
      if (debug)
	fprintf (stderr, "%s:%d: '%s' cannot be a constant anymore\n", 
		 fname, line, var->v);
      var->maybeconst = 0;
    }
  }

  void declaration (char * var, char * text) {
    var_t * v;
    if (!strcmp (&var[strlen(var)-2], "[]")) {
      // automatic
      var[strlen(var)-2] = '\0';
      v = varpush (var, vartype, scope, 0);
      v->automatic = 1;
      v->symmetric = varsymmetric;
      v->face = varface;    
      v->vertex = varvertex;    
      if (varconst)
	v->constant = strdup (varconst);
      fputs (text, yyout);
      fputc ('=', yyout);
      new_field (v);
    }
    else {
      v = varpush (var, vartype, scope, varmaybeconst);
      v->symmetric = varsymmetric;
      v->face = varface;    
      v->vertex = varvertex;
      fputs (text, yyout);
    }
    if (debug)
      fprintf (stderr, "%s:%d: declaration: %s\n", fname, line, var);
  }

  static int homogeneize (FILE * in, FILE * fp)
  {
    /* replace neumann/dirichlet(...) with neumann/dirichlet(0) */
    char * cond[2] = {"dirichlet", "neumann"};
    char str[1024];
    while (fgets (str, 1024, in)) {
      char * s = str;
      while (*s != '\0') {
	int i;
	for (i = 0; i < 2; i++)
	  if (!strncmp (s, cond[i], strlen (cond[i]))) {
	    fputs (cond[i], fp);
	    s += strlen (cond[i]);
	    while (strchr (" \t\v\n\f", *s))
	      s++;
	    if (*s == '(') {
	      fputs ("_homogeneous()", fp);
	      s++;
	      int para = 1;
	      while (para > 0 && *s != '\0') {
		if (*s == '(') para++;
		if (*s == ')') para--;
		s++;
	      }
	      if (para > 0)
		return 1;
	    }
	    break;
	  }
	fputc (*s++, fp);
      }
    }
    return 0;
  }
  
  void boundary_body (void) {
    FILE * fp = dopen ("_inboundary.h", "r");
    int c;
    while ((c = fgetc (fp)) != EOF)
      fputc (c, yyout);
    fclose (fp);
  }

  int yyerror(const char * s);

  void boundary_homogeneous_body (void) {
    FILE * fp = dopen ("_inboundary.h", "r");
    if (homogeneize (fp, yyout))
      yyerror ("expecting ')'");
    fclose (fp);
  }

#define nonspace(s) { while (strchr(" \t\v\n\f", *s)) s++; }
#define space(s) { while (!strchr(" \t\v\n\f", *s)) s++; }

#define YY_INPUT(buf,result,max_size)			      \
  {							      \
    int c = fgetc(yyin);				      \
    if (c == '\n') { line++; }				      \
    result = (c == EOF) ? YY_NULL : (buf[0] = c, 1);	      \
  }

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
      int i;
      for (i = 0; i < nreduct; i++)
	if (strcmp (reduction[i], "+"))
	  fprintf (yyout, "double _%s = %s; ", reductvar[i], reductvar[i]);
      fputs ("\n#undef _OMPEND\n#define _OMPEND ", yyout);
      for (i = 0; i < nreduct; i++) {
	if (strcmp (reduction[i], "+"))
	  fprintf (yyout,
		   "OMP(omp critical) if (_%s %s %s) %s = _%s; ",
		   reductvar[i], strcmp(reduction[i], "min") ? ">" : "<",
		   reductvar[i], reductvar[i], reductvar[i]);
	fprintf (yyout,
		 "mpi_all_reduce (%s, MPI_DOUBLE, %s); ",
		 reductvar[i], 
		 !strcmp(reduction[i], "min") ? "MPI_MIN" : 
		 !strcmp(reduction[i], "max") ? "MPI_MAX" : 
		 "MPI_SUM");
      }
      fprintf (yyout, "\n#line %d\n", foreach_line);
    }
    yyout = dopen ("_foreach_body.h", "w");
    if (inforeach_face) {
      foreach_face_xy = face_xy;
      if (nreduct == 0) { // foreach_face (x)
	FILE * fp = dopen ("_foreach.h", "r");
	int c;
	while ((c = fgetc (fp)) != EOF)
	  if (c == 'x')
	    foreach_face_xy = face_x;
	  else if (c == 'y')
	    foreach_face_xy = face_y;
	fclose (fp);
      }
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
	     "static int %s (const int i, const double t, Event * _ev) { ",
	     eventfunc[nevents]);
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
  if (foreachdim != 1 || scope != 0 || infunctionproto)
    ECHO;
  if (infunction && functionpara == 1 && scope == functionscope)
    infunction_declarations();
  if (inmain == 1 && scope == 0) {
    fputs (" init_solver(); ", yyout);
    inmain = 2;
  }
  scope++;
}

"}" {
  scope--;
  if (scope < 0)
    return yyerror ("mismatched '}'");
  if (inmain == 2 && scope == 0) {
    fputs (" free_solver(); ", yyout);
    inmain = 0;
  }    
  varpop();
  if (foreachdim && scope == foreachdim - 1) {
    if (scope != 0 || infunctionproto)
      ECHO;
    endforeachdim ();
  }
  else if (!inattr || scope != attrscope)
    ECHO;
  if (foreach_child && foreach_child_scope == scope) {
    fputs (" end_foreach_child(); }", yyout);
    foreach_child = 0;
  }
  if (infunction && scope <= functionscope) {
    infunction = 0;
    if (debug)
      fprintf (stderr, "%s:%d: outfunction\n", fname, line);
  }
  if (inforeach && scope == foreachscope)
    endforeach ();
  else if (inevent && scope == eventscope)
    endevent();
  if (inattr && scope == attrscope)
    endattr ();
  if (scope == 0)
    infunctionproto = 0;
}


(foreach_child|foreach_child_direction) {
  fputs (" { ", yyout);
  ECHO;
  foreach_child = 1; foreach_child_scope = scope; foreach_child_para = para;
}

foreach{ID}* {
  fputs (" { ", yyout);
  strcpy (foreachs, yytext);
  inforeach = 1; foreachscope = scope; foreachpara = para;
  free (foreachconst); 
  foreachconst = NULL;
  nmaybeconst = 0;
  nreduct = 0;
  foreach_line = line;
  foreachfp = yyout;
  yyout = dopen ("_foreach.h", "w");
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

attribute{WS}+"{" {
  inattr = 1; attrscope = scope++; attrline = line;
  attrfp = yyout;
  yyout = dopen ("_attribute.h", "w");
}

; {
  int insthg = 0;
  if (scope == 0)
    infunctionproto = 0;
  if (foreachdim && scope == foreachdim - 1 && para == foreachdimpara) {
    ECHO; insthg = 1;
    endforeachdim ();
  }
  if (foreach_child && scope == foreach_child_scope && para == foreach_child_para) {
    ECHO; insthg = 1;
    fputs (" end_foreach_child(); }", yyout); foreach_child = 0;
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
    fclose (yyout);

    yyout = boundaryheader;
    maybeconst_combinations (line, boundary_body);
    fputs (" return 0.; } ", yyout);
    fprintf (yyout, 
	     "static double _boundary%d_homogeneous (Point point, scalar _s) {",
	     nboundary);
    boundary_staggering (boundarydir, boundarycomponent, yyout);
    fputs (" POINT_VARIABLES; ", yyout);
    maybeconst_combinations (line, boundary_homogeneous_body);
    fputs (" return 0.; }\n", yyout);

    yyout = boundaryfp;
    if (scope == 0) {
      /* file scope */
      fprintf (yyout, 
	       "static void _set_boundary%d (void) { ",
	       nboundary);
      boundaryindex[nsetboundary++] = nboundary;
    }
    /* function/file scope */
    fprintf (yyout,
	     "_attribute[%s].boundary[%s] = _boundary%d; "
	     "_attribute[%s].boundary_homogeneous[%s] = "
	     "_boundary%d_homogeneous;",
	     boundaryvar, boundarydir, nboundary,
	     boundaryvar, boundarydir, nboundary);
    if (scope == 0)
      /* file scope */
      fputs (" } ", yyout);
    nboundary++;
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
	     "static int %s_expr%d (int * ip, double * tp, Event * _ev) { "
	     "  int i = *ip; double t = *tp; "
	     "  int ret = (", eventfunc[nevents], inevent++);
  }
  else if (inevent == 4 && scope == eventscope && para == eventpara - 1) {
    ECHO;
    endevent ();
  }
  else if (inreturn) {
    fputs ("; ", yyout);
    delete_automatic (0);
    if (inmain == 2)
      fputs (" free_solver(); ", yyout);
    fputs (" return _ret; }", yyout);
    inreturn = 0;
  }
  else if (!insthg)
    ECHO;
  varpop();
  invardecl = varmaybeconst = 0;
}

{ID}+{WS}*[=]{WS}*new{WS}+(symmetric|face|vertex){0,1}*{WS}*(scalar|vector|tensor) |
[=]{WS}*new{WS}+(symmetric|face|vertex){0,1}{WS}*(scalar|vector|tensor) {
  char * type = strchr (yytext, '=');
  type = strstr (type, "new"); space(type); nonspace(type);
  char * symmetric = strstr (type, "symmetric");
  if (symmetric) {
    space(type); nonspace(type);
  }
  char * face = strstr (type, "face");
  if (face) {
    space(type); nonspace(type);
  }
  char * vertex = strstr (type, "vertex");
  if (vertex) {
    space(type); nonspace(type);
  }
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
  var->symmetric = (symmetric != NULL);
  var->face   = (face != NULL);
  var->vertex = (vertex != NULL);
  fputc ('=', yyout);
  new_field (var);
  if (debug)
    fprintf (stderr, "%s:%d: new %s%s: %s\n", 
	     fname, line, 
	     var->symmetric ? "symmetric " : 
	     var->face ? "face " :
	     var->vertex ? "vertex " :
	     "", 
	     type, var->v);
}

{ID}+{WS}*[=]{WS}*automatic{WS}*[(][^)]*[)] |
[=]{WS}*automatic{WS}*[(][^)]*[)] {
  var_t * var;
  if (yytext[0] == '=') {
    var = vartop();
    assert (var && var->scope == scope);
  }
  else {
    char * s = yytext; space (s); *s++ = '\0';
    var = varlookup (yytext, strlen(yytext));
    if (!var) {
      fprintf (stderr, "%s:%d: '%s' undeclared", fname, line, yytext);
      return 1;
    }
    fprintf (yyout, "%s ", yytext);
    yytext = s;
  }
  var->automatic = 1;
  char * arg = strchr (yytext, '(');
  if (arg[1] == ')')
    fputc ('=', yyout);
  else {
    var->conditional = malloc (strlen(arg) + 5);
    sprintf (var->conditional, "%s%s", arg, 
	     var->type == vector ? ".x" : var->type == tensor ? ".x.x" : "");
    fprintf (yyout, "= %s ? %s :", var->conditional, arg);
  }
  new_field (var);
}

\({WS}*const{WS}*\)    varmaybeconst = 1;

symmetric{WS}+tensor{WS}+[a-zA-Z0-9_\[\]]+ |
face{WS}+vector{WS}+[a-zA-Z0-9_\[\]]+ |
vertex{WS}+scalar{WS}+[a-zA-Z0-9_\[\]]+ |
(scalar|vector|tensor){WS}+[a-zA-Z0-9_\[\]]+ {
  varsymmetric = (strstr(yytext, "symmetric") != NULL);
  varface = (strstr(yytext, "face") != NULL);
  varvertex = (strstr(yytext, "vertex") != NULL);
  varconst = NULL;
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
  char * text = var;
  var = &var[7];
  nonspace (var);
  if (*var == '[') {
    // scalar [..
    for (; *var != '\0'; var++)
      if (*var == '[') brack++;
      else if (*var == ']') brack--;
    fputs (text, yyout);
  }
  else if (para == 0) { /* declaration */
    declaration (var, text);
    invardecl = scope + 1;
  }
  else if (para == 1) { /* function prototype (no nested functions) */
    fputs (text, yyout);
    if (debug)
      fprintf (stderr, "%s:%d: proto: %s\n", fname, line, var);
    varpush (var, vartype, scope + 1, varmaybeconst);
    invardecl = varmaybeconst = 0;
  }
  else
    fputs (text, yyout);
}

const{WS}+(symmetric{WS}+|face{WS}+|vertex{WS}+|{WS}*)(scalar|vector|tensor){WS}+{ID}+{WS}*=[^;]+ {
  ECHO;
  char * s = strchr (yytext, '='); s--;
  while (strchr (" \t\v\n\f", *s)) s--; s++; *s = '\0';
  fprintf (stderr, "%s:%d: warning: did you mean '%s[]'?\n", 
	   fname, line, yytext);
}

const{WS}+(symmetric{WS}+|face{WS}+|vertex{WS}+|{WS}*)(scalar|vector|tensor){WS}+{ID}+\[{WS}*\]{WS}*=[^;]+ {
  // const scalar a[] = 1.;
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
  char * text = var;
  var = &var[7];
  nonspace (var);
  char * cst = strchr (var, ']'); *++cst = '\0';
  cst = strchr (cst + 1, '=') + 1;
  if (para != 0)
    return yyerror ("constant fields can only appear in declarations");
  varconst = cst;
  varsymmetric = varface = varvertex = 0;
  declaration (var, text);
}

,{WS}*[a-zA-Z0-9_\[\]]+ {
  if (invardecl == scope + 1) {
    char * var = &yytext[1];
    nonspace (var);
    declaration (var, yytext);
  }
  else
    REJECT;
}

return{WS} {
  char * list = automatic_list (0, 1);
  if (list) {
    // returning from a function: delete automatic fields before returning
    // note that this assumes that the function scope is always 1 
    // (i.e. no nested functions allowed).
    free (list);
    inreturn = 1;
    fprintf (yyout, "{ %s _ret = ", return_type);
  }
  else
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
      if (var->constant)
	fprintf (yyout, "_val_constant(%s", boundaryrep(yytext));
      else if (var->maybeconst) {
	int found = 0, i;
	for (i = 0; i < nmaybeconst && !found; i++)
	  found = !strcmp (foreachconst[i]->v, var->v);
	if (!found) {
	  foreachconst = realloc (foreachconst, 
				  ++nmaybeconst*sizeof (var_t *));
	  foreachconst[nmaybeconst - 1] = var;
	  if (debug)
	    fprintf (stderr, "%s:%d: '%s' may be const\n", 
		     fname, line, var->v);
	}
	char * macro = macro_from_var (yytext);
	fprintf (yyout, "val_%s(%s", macro, yytext);
	free (macro);
      }
      else
	fprintf (yyout, "val(%s", boundaryrep(yytext));
      if (yytext[yyleng-1] == ']')
	/* v[] */
	fputs (",0,0)", yyout);
      else {
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
	    fprintf (stderr, "%s:%d: error: expecting '[' not '", 
		     fname, line);
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
  var_t * var;
  if ((var = varlookup (yytext, s - yytext))) {
    s = yytext;
    while (!strchr(" \t\v\n\f[", *s)) s++;
    *s++ = '\0';
    nonspace (s);
    strcpy (boundarydir, s);
    s = boundarydir;
    while (!strchr(" \t\v\n\f]", *s)) s++;
    *s++ = '\0';
    char * s1 = yytext;
    boundarycomponent = 0;
    while (*s1 != '\0') {
      if (*s1 == '.') {
	s1++;
	boundarycomponent = *s1;
      }
      else
	s1++;
    }
    if (!var->face)
      boundarycomponent = 0;
    boundaryfp = yyout;

    yyout = boundaryheader;
    fprintf (yyout,
	     "#line %d \"%s\"\n"
	     "static double _boundary%d (Point point, scalar _s) {",
	     line, fname, nboundary);
    boundary_staggering (boundarydir, boundarycomponent, yyout);
    fputs (" POINT_VARIABLES; ", yyout);
    nmaybeconst = 0;
    inboundary = inforeach_boundary = 1;
    strcpy (boundaryvar, yytext);
    inforeach = 2;

    yyout = dopen ("_inboundary.h", "w");
    fputs ("return ", yyout);
  }
  else
    REJECT;
}

dirichlet{WS}*[(] {
  para++;
  if (inboundary &&
      ((boundarycomponent == 'x' && (!strcmp(boundarydir, "left") ||
				     !strcmp(boundarydir, "right"))) ||
       (boundarycomponent == 'y' && (!strcmp(boundarydir, "top") ||
				     !strcmp(boundarydir, "bottom")))))
    fputs (&yytext[9], yyout);
  else
    ECHO;
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

for{WS}*[(]{WS}*(scalar|vector|tensor){WS}+{ID}+{WS}+in{WS}+ {
  /* for (scalar .. in .. */
  char * s = strchr (&yytext[3], 'r'); s++;
  int vartype = s[-3] == 's' ? tensor : s[-3] == 't' ? vector : scalar;
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
	   "for (%s %s = *%s, *_i%d = %s; *((scalar *)&%s) >= 0; %s = *++_i%d%s",
	   vartype == scalar ? "scalar" : 
	   vartype == vector ? "vector" : 
                               "tensor", 
	   id,
	   list, i, list, 
	   id, 
	   id, i, last);
  free (list);
  i++;
  varpush (id, vartype, scope + 1, 0);
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
  static int index = 0, i;
  char * sc = NULL, * vc = NULL, * tc = NULL;
  for (i = 0; i < nid; i++) {
    var_t * var = varlookup (id[i], strlen(id[i]));
    if (!var) {
      fprintf (stderr, 
	       "%s:%d: error: '%s' is not a scalar, vector or tensor\n", 
	       fname, line, id[i]);
      return 1;
    }
    fprintf (yyout, "%s * _i%d = %s; ", 
	     var->type == vector ? "vector" : 
	     var->type == scalar ? "scalar" :
	                           "tensor",
	     index + i, list[i]);
    if (!sc && var->type == scalar) sc = id[i];
    if (!vc && var->type == vector) vc = id[i];
    if (!tc && var->type == tensor) tc = id[i];
  }
  fprintf (yyout, "if (%s) for (%s = *%s", list[0], id[0], list[0]);
  for (i = 1; i < nid; i++)
    fprintf (yyout, ", %s = *%s", id[i], list[i]);
  fprintf (yyout, "; *((scalar *)&%s) >= 0; ", id[0]);
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

event{WS}+{ID}+{WS}*[(] {
  /* event (... */
  char ids[80], * id = ids;
  strncpy (id, yytext, 80);
  space (id); nonspace (id);
  char * s = id; while (!strchr(" \t\v\n\f(", *s)) s++; *s = '\0';
  eventparent[nevents] = eventchild[nevents] = -1;
  eventlast[nevents] = 0;
  eventid[nevents] = strdup (id);
  // check whether an event with the same name already exists
  // if it does, append a number (up to 9)
  char i = '0', found = 1;
  int len = strlen(id), lastfound = -1;
  id[len+2] = '\0';
  while (i <= '9' && found) {
    found = 0;
    int j;
    for (j = 0; j < nevents && !found; j++)
      if (!strcmp (id, eventfunc[j])) {
	lastfound = j;
	found = 1;
      }
    if (found) {
      id[len] = '_';
      id[len+1] = i++;
    }
  }
  if (lastfound >= 0) {
    eventparent[lastfound] = nevents;
    eventchild[nevents] = lastfound;
  }
  fprintf (yyout, 
	   "static int %s_expr%d (int * ip, double * tp, Event * _ev) {"
	   "  int i = *ip; double t = *tp;"
	   "  int ret = (", id, inevent++);
  eventscope = scope; eventpara = ++para;
  eventarray[nevents] = 0;
  eventfile[nevents] = strdup (fname);
  eventline[nevents] = line;
  eventfunc[nevents] = strdup (id);
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

,{WS}*last {
  if (inevent == 1)
    eventlast[nevents] = 1;
  else
    ECHO;
}

end {
  if (inevent == 1)
    fputs ("1234567890", yyout);
  else
    ECHO;
}

{ID}+{WS}+{ID}+{WS}*\( {
  if (scope == 0 && para == 0) {
    // function prototype
    para++;
    ECHO;
    char * s1 = yytext; space (s1); *s1++ = '\0';
    nonspace (s1);
    char * s2 = s1; while (!strchr (" \t\v\n\f(", *s2)) s2++; *s2 = '\0';
    free (return_type);
    return_type = strdup (yytext);
    infunctionproto = 1;
    if (debug)
      fprintf (stderr, "%s:%d: function '%s' returns '%s'\n", 
	       fname, line, s1, yytext);
    if (!strcmp (s1, "main") && !strcmp (yytext, "int"))
      inmain = 1;
  }
  else
    REJECT;
}

foreach_dimension{WS}*[(]{WS}*[)] {
  foreachdim = scope + 1; foreachdimpara = para;
  foreachdimfp = yyout;
  yyout = dopen ("_dimension.h", "w");
  foreachdimline = line;
}

reduction{WS}*[(](min|max|\+):{ID}+[)] {
  if (strchr(yytext, '+'))
    ECHO;
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
  var_t * var;
  if (inforeach) {
    int i;
    for (i = 0; i < nreduct; i++)
      if (strcmp ("+", reduction[i]) && !strcmp (yytext, reductvar[i])) {
	fputc ('_', yyout);
	break;
      }
    ECHO;
  }
  else if (scope == 0 && para == 0 &&
	   (var = varlookup (yytext, strlen(yytext))) &&
	   var->constant) {
    // replace global scalar/vector constants
    if (var->type == scalar)
      fprintf (yyout, "(_NVARMAX + %d)", var->i[0]);
    else if (var->type == vector)
      fprintf (yyout, "{_NVARMAX + %d,_NVARMAX + %d}", var->i[0], var->i[1]);
    else
      assert (0);
  }
  else
    ECHO;
}

{SCALAR}[.][a-wA-Z_0-9]+ {
  // scalar attributes
  char * dot1 = strchr (yytext, '.');
  var_t * var = varlookup (yytext, strlen(yytext) - strlen(dot1));
  if (var) {
    char * dot = dot1;
    while (dot1) {
      dot = dot1;
      dot1 = strchr (dot1 + 1, '.');
    }
    dot[0] = '\0'; dot++;
    fprintf (yyout, "_attribute[%s].%s", yytext, dot);
    if (debug)
      fprintf (stderr, "%s:%d: _attribute[%s].%s\n", fname, line, yytext, dot);
  }
  else
    ECHO;
}

{ID}+{WS}+{ID}+{WS}*[(]{WS}*struct{WS}+{ID}+{WS}+{ID}+{WS}*[)] {
  // function declaration with struct argument
  ECHO;
  char * s = yytext;
  space (s); *s++ = '\0'; nonspace (s);
  free (return_type);
  return_type = strdup (yytext);
  args = realloc (args, sizeof (char *)*++nargs);
  argss = realloc (argss, sizeof (char *)*nargs);
  char * s1 = s; while (!strchr(" \t\v\n\f(", *s1)) s1++;
  *s1++ = '\0'; s1 = strstr (s1, "struct"); space (s1); nonspace (s1);
  char * s2 = s1; space (s2); *s2++ = '\0';
  char * s3 = s2; space (s3); if (*s3 == '\0') s3[-1] = '\0';
  args[nargs-1] = strdup (s);
  argss[nargs-1] = strdup (s1);
  varpush (s2, struct_type, scope + 1, 0);
  if (debug)
    fprintf (stderr, 
	     "%s:%d: function '%s' with struct '%s' '%s' returns '%s'\n", 
	     fname, line, s, s1, s2, return_type);
}

{ID}+{WS}*[(]{WS}*{ID}+{WS}*[)] {
  // function call with a single 'args' assignment
  if (inarg)
    REJECT;
  char * s = yytext; space (s);
  int len = s - yytext, i;
  for (i = 0; i < nargs && !inarg; i++) {
    if (strlen(args[i]) == len && !strncmp (args[i], yytext, len)) {
      char * s = strchr (yytext, '('); s++; nonspace (s);
      char * s1 = s; space (s1); *s1 = '\0';
      s1 = strchr (s, ')'); if (s1) *s1 = '\0';
      if (debug)
	fprintf (stderr, "%s:%d: single argument '%s'\n", fname, line, s);
      if (get_vartype (s) == struct_type)
	fprintf (yyout, "%s (%s)", args[i], s);
      else
	fprintf (yyout, "%s ((struct %s){%s})", args[i], argss[i], s);
      inarg = 1;
    }
  }
  if (!inarg)
    REJECT;
  inarg = 0;
}

{ID}+{WS}*[(] {
  // function call with multiple 'args' assignment
  if (inarg)
    REJECT;
  char * s = yytext;
  while (!strchr(" \t\v\n\f(", *s)) s++;
  int len = s - yytext, i;
  for (i = 0; i < nargs && !inarg; i++)
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

#{SP}+[0-9]+{SP}+\"[^\"]+\" {
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

static{WS}+FILE{WS}*[*]{WS}*{ID}+{WS}*= {
  // static FILE * fp = ...
  ECHO;
  char * s = yytext;
  while (*s != '*') s++;
  s++;
  nonspace (s);
  char * id = s;
  space (s);
  *s = '\0';
  fprintf (yyout, "NULL; if (!%s) %s = ", id, id);
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
  line = 1, scope = para = 0;
  inforeach = foreachscope = foreachpara = 
    inforeach_boundary = inforeach_face = 0;
  invardecl = 0;
  inval = invalpara = 0;
  brack = inarray = 0;
  inevent = inreturn = inattr = 0;
  foreachdim = 0;
  foreach_child = 0;
  inboundary = nboundary = nsetboundary = 0;
  infunction = 0;
  infunctionproto = inmain = 0;
  inarg = 0;
  boundaryheader = dopen ("_boundary.h", "w");
  int ret = yylex();
  if (!ret) {
    if (scope > 0)
      ret = yyerror ("mismatched '{'");
    else if (para > 0)
      ret = yyerror ("mismatched '('");
  }
  fclose (boundaryheader);
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

void write_event (int i, FILE * fout)
{
  char * func = eventfunc[i];
  fprintf (fout, "  event_register ((Event){ 0, %d, %s, {", nexpr[i], func);
  int j;
  for (j = 0; j < nexpr[i] - 1; j++)
    fprintf (fout, "%s_expr%d, ", func, j);
  fprintf (fout, "%s_expr%d}, ", func, j);
  if (eventarray[i] == 'i')
    fprintf (fout, "%s_array, ", func);
  else
    fprintf (fout, "((void *)0), ");
  if (eventarray[i] == 't')
    fprintf (fout, "%s_array,\n", func);
  else
    fprintf (fout, "((void *)0),\n");
  fprintf (fout, "    \"%s\", %d, \"%s\"});\n", eventfile[i], 
	   nolineno ? 0 : eventline[i], eventid[i]);
}

static const char * typestr (var_t var) {
  return (var.type == scalar ? "scalar" :
	  var.type == vector ? "vector" :
	  var.type == tensor ? "tensor" : "");
}

void compdir (FILE * fin, FILE * fout, FILE * swigfp, 
	      char * swigname,
	      char * grid)
{
  if (endfor (fin, fout))
    cleanup (1, dir);
  fclose (fout);

  fout = dopen ("_boundarydecl.h", "w");
  int i;
  for (i = 0; i < nboundary; i++)
    fprintf (fout, 
       "static double _boundary%d (Point point, scalar _s);\n"
       "static double _boundary%d_homogeneous (Point point, scalar _s);\n", 
	     i, i);
  fclose (fout);

  fout = dopen ("_grid.h", "w");
  /* new variables */
  fprintf (fout,
	   "int datasize = %d*sizeof (double);\n",
	   nvar);
  /* attributes */
  FILE * fp = dopen ("_attributes.h", "a");
  fputs ("\n"
	 "} _Attributes;\n"
	 "_Attributes * _attribute;\n",
	 fp);
  fclose (fp);
  /* event functions */
  for (i = 0; i < nevents; i++) {
    char * id = eventfunc[i];
    fprintf (fout, 
	     "static int %s (const int i, const double t, Event * _ev);\n", id);
    int j;
    for (j = 0; j < nexpr[i]; j++)
      fprintf (fout,
	       "static int %s_expr%d (int * ip, double * tp, Event * _ev);\n",
	       id, j);
    if (eventarray[i])
      fprintf (fout, "static %s %s_array[] = %s,-1};\n", 
	       eventarray[i] == 'i' ? "int" : "double", id,
	       eventarray_elems[i]);
  }
  /* boundaries */
  for (i = 0; i < nsetboundary; i++)
    fprintf (fout, "static void _set_boundary%d (void);\n", 
	     boundaryindex[i]);
  /* init_solver() */
  fputs ("void init_solver (void) {\n", fout);
  /* MPI */
  fputs ("#if _MPI\n"
	 "  void mpi_init();\n"
	 "  mpi_init();\n"
	 "#endif\n",
	 fout);
  /* events */
  fputs ("  Events = malloc (sizeof (Event));\n"
	 "  Events[0].last = 1;\n", fout);
  int last;
  for (last = 0; last <= 1; last++)
    for (i = 0; i < nevents; i++)
      if (eventchild[i] < 0 && eventlast[i] == last) {
	int j = i;
	while (j >= 0) {
	  write_event (j, fout);
	  j = eventparent[j];
	}
      }
  fputs ("  { 1 }\n};\n", fout);
  /* boundaries */
  for (int i = 0; i < nsetboundary; i++)
    fprintf (fout, "static void _set_boundary%d (void);\n", 
	     boundaryindex[i]);
  /* methods */
  fputs ("void init_solver (void) {\n", fout);
  /* scalar methods */
  fputs ("  _method = calloc (datasize/sizeof(double), sizeof (Methods));\n", 
	 fout);
  /* list of all scalars */
  fprintf (fout, 
	   "  all = malloc (sizeof (scalar)*%d);\n"
	   "  for (int i = 0; i < %d; i++)\n"
	   "    all[i] = i;\n"
	   "  all[%d] = -1;\n", nvar + 1, nvar, nvar);
  fputs ("#if _GNU_SOURCE\n"
	 "  set_fpe();\n"
	 "#endif\n", fout);
  if (catch)
    fputs ("  catch_fpe();\n", fout);
  fprintf (fout, "  %s_methods();\n", grid);
  for (i = varstack; i >= 0; i--) {
    var_t var = _varstack[i];
    if (var.i[0] >= 0) {
      if (var.constant) {
	// global constants
	if (var.type == scalar)
	  fprintf (fout, 
		   "  init_const_scalar (_NVARMAX+%d, \"%s\", %s);\n",
		   var.i[0], var.v, var.constant);
	else if (var.type == vector)
	  fprintf (fout, 
		   "  init_const_vector ((vector){_NVARMAX+%d,_NVARMAX+%d},"
		   " \"%s\", (double [])%s);\n",
		   var.i[0], var.i[1], var.v, var.constant);
	else
	  assert (0);
      }
      // global variables
      else if (var.type == scalar)
	fprintf (fout, "  init_scalar (%d, \"%s\");\n",
		 var.i[0], var.v);
      else if (var.type == vector)
	fprintf (fout, "  init_%svector ((vector){%d,%d}, \"%s\");\n",
		 var.face ? "face_" : 
		 var.vertex ? "vertex_" : 
		 "",
		 var.i[0], var.i[1], var.v);
      else if (var.type == tensor)
	fprintf (fout, "  init_tensor ((tensor){{%d,%d},{%d,%d}}, \"%s\");\n", 
		 var.i[0], var.i[1], var.i[2], var.i[3], var.v);
      else
	assert (0);
    }
  }
  for (i = 0; i < nsetboundary; i++)
    fprintf (fout, "  _set_boundary%d();\n", boundaryindex[i]);
  fputs ("  init_grid (N);\n"
	 "}\n", fout);
  fclose (fout);
  
  /* SWIG interface */
  if (swigfp) {
    int i;
    fputs ("\n%{\n", swigfp);
    for (i = varstack; i >= 0; i--) {
      var_t var = _varstack[i];
      if (var.i[0] >= 0 && !var.constant)
	fprintf (swigfp, "  extern %s %s;\n", typestr (var), var.v);
    }
    fputs ("%}\n\n", swigfp);
    for (i = varstack; i >= 0; i--) {
      var_t var = _varstack[i];
      if (var.i[0] >= 0 && !var.constant)
	fprintf (swigfp, "extern %s %s;\n", typestr (var), var.v);
    }
    fputs ("\n%pythoncode %{\n", swigfp);
    for (i = varstack; i >= 0; i--) {
      var_t var = _varstack[i];
      if (var.i[0] >= 0 && !var.constant)
	fprintf (swigfp, "%s = %s(_%s.cvar.%s)\n",
		 var.v, typestr (var), swigname, var.v);
    }
    fputs ("%}\n", swigfp);
    fclose (swigfp);
  }
}

int main (int argc, char ** argv)
{
  char * cc = getenv ("CC99"), command[1000], command1[1000] = "";
  if (cc == NULL)
    strcpy (command, CC99);
  else
    strcpy (command, cc);
  char * file = NULL;
  int i, dep = 0, tags = 0, source = 0, swig = 0;
  for (i = 1; i < argc; i++) {
    if (!strncmp (argv[i], "-grid=", 6))
      ;
    else if (!strcmp (argv[i], "-MD"))
      dep = 1;
    else if (!strcmp (argv[i], "-tags"))
      tags = 1;
    else if (!strcmp (argv[i], "-python"))
      swig = 1;
    else if (!strcmp (argv[i], "-debug"))
      debug = 1;
    else if (!strcmp (argv[i], "-events"))
      events = 1;
    else if (!strcmp (argv[i], "-catch"))
      catch = 1;
    else if (!strcmp (argv[i], "-source"))
      source = 1;
    else if (catch && !strncmp (argv[i], "-O", 2))
      ;
    else if (!strcmp (argv[i], "-nolineno"))
      nolineno = 1;
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
    FILE * swigfp = NULL;
    char swigname[80] = "";
    if (swig) {
      strcpy (swigname, file);
      char * dot = strchr (swigname, '.');
      *dot = '\0'; strcat (swigname, ".i");
      swigfp = fopen (swigname, "a");
      if (!swigfp) {
	fprintf (stderr, "qcc: could not open '%s': ", swigname);
	return 1;
      }
      *dot = '\0';
    }
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
      FILE * fp = dopen ("_attributes.h", "w");
      fputs ("typedef struct {\n", fp);
      fclose (fp);
      FILE * fout = dopen (cpp, "w");
      if (swig)
	fputs ("@include <Python.h>\n", fout);
      fputs ("@if _XOPEN_SOURCE < 700\n"
	     "  @undef _XOPEN_SOURCE\n"
	     "  @define _XOPEN_SOURCE 700\n"
	     "@endif\n"
	     "@if _GNU_SOURCE\n"
	     "@include <stdint.h>\n"
	     "@include <string.h>\n"
	     "@include <fenv.h>\n"
	     "@endif\n",
	     fout);
      if (catch)
	fputs ("#define TRASH 1\n"
	       "#define _CATCH last_point = point;\n", 
	       fout);
      else
	fputs ("#define _CATCH\n", fout);
      fputs ("#include \"common.h\"\n", fout);
      /* catch */
      if (catch)
	fputs ("void catch_fpe (void);\n", fout);
      /* undefined value */
      /* Initialises unused memory with "signaling NaNs".  
       * This is probably not very portable, tested with
       * gcc (Debian 4.4.5-8) 4.4.5 on Linux 2.6.32-5-amd64.
       * This blog was useful:
       *   http://codingcastles.blogspot.co.nz/2008/12/nans-in-c.html 
       */
      fputs ("@if _GNU_SOURCE\n"
	     "double undefined;\n"
	     "static void set_fpe (void) {\n"
	     "  int64_t lnan = 0x7ff0000000000001;\n"
	     "  assert (sizeof (int64_t) == sizeof (double));\n"
	     "  memcpy (&undefined, &lnan, sizeof (double));\n"
	     "  feenableexcept (FE_DIVBYZERO|FE_INVALID);\n"
	     "}\n"
	     "@else\n"
	     "@  define undefined DBL_MAX\n"
	     "@endif\n", fout);
      /* grid */
      if (default_grid)
	fprintf (fout, "#include \"grid/%s.h\"\n", grid);
      /* declaration of boundary condition fonctions */
      fputs ("@include \"_boundarydecl.h\"\n", fout);
      char s[81];
      while (fgets (s, 81, fin)) {
	// replace '#include <' with '@include <'
	char * s1 = s; while (strchr (" \t", *s1)) s1++;
	if (*s1 == '#') {
	  char *s2 = s1 + 1; while (strchr (" \t", *s2)) s2++;
	  if (!strncmp (s2, "include", 7)) {
	    while (!strchr (" \t", *s2)) s2++;
	    while (strchr (" \t", *s2)) s2++;
	    if (*s2 == '<')
	      *s1 = '@';
	  }
	}
	fputs (s, fout);
      }
      if (swigfp)
	fputs ("@include \"python.h\"\n", fout);
      fclose (fout);
      fclose (fin);
      fout = dopen (file, "w");
      if (!fout) {
	perror (file);
	cleanup (1, dir);
      }

      char preproc[1000], * cppcommand = getenv ("CPP99");
      strcpy (preproc, "cd ");
      strcat (preproc, dir);
      strcat (preproc, " && ");
      if (!cppcommand && strcmp (CPP99, ""))
	cppcommand = CPP99;
      if (!cppcommand) {
	strcat (preproc, command);
	strcat (preproc, " -E");
      }
      else
	strcat (preproc, cppcommand);
      strcat (preproc, " -I. -I");
      strcat (preproc, LIBDIR);
      strcat (preproc, " ");
      if (events) {
	strcat (preproc, " -DDEBUG_EVENTS=1 -DBASILISK=\"\\\"");
	strcat (preproc, BASILISK);
	strcat (preproc, "\\\"\" ");
      }
      strcat (preproc, cpp);
      if (debug)
	fprintf (stderr, "preproc: %s\n", preproc);

      fin = popen (preproc, "r");
      if (!fin) {
	perror (preproc);
	cleanup (1, dir);
      }

      compdir (fin, fout, swigfp, swigname, grid);
      int status = pclose (fin);
      if (status == -1 ||
	  (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				    WTERMSIG (status) == SIGQUIT)))
	cleanup (1, dir);

      fout = dopen ("_tmp", "w");

      fin = dopen (file, "r");
      char line[1024];
      int c;
      // includes _attributes.h
      while (fgets (line, 1024, fin))
	if (!strcmp (line, "#include \"_attributes.h\"\n"))
	  break;
	else
	  fputs (line, fout);
      fp = dopen ("_attributes.h", "r");
      while ((c = fgetc (fp)) != EOF)
	fputc (c, fout);
      fclose (fp);
      // includes _boundarydecl.h
      while (fgets (line, 1024, fin))
	if (!strcmp (line, "#include \"_boundarydecl.h\"\n"))
	  break;
	else
	  fputs (line, fout);
      fp = dopen ("_boundarydecl.h", "r");
      while ((c = fgetc (fp)) != EOF)
	fputc (c, fout);
      fclose (fp);
      // rest of the file
      while (fgets (line, 1024, fin))
	fputs (line, fout);
      fclose (fin);

      fin = dopen ("_boundary.h", "r");
      while ((c = fgetc (fin)) != EOF)
	fputc (c, fout);
      fclose (fin);

      fin = dopen ("_grid.h", "r");
      while ((c = fgetc (fin)) != EOF)
	fputc (c, fout);
      fclose (fin);

      fclose (fout);

      char src[80], dst[80];
      strcpy (src, dir); strcat (src, "/_tmp");
      if (source) {
	strcpy (dst, "_");
      }
      else {
	strcpy (dst, dir); strcat (dst, "/");
      }
      strcat (dst, file);
      rename (src, dst);

      strcat (command, " -I");
      strcat (command, LIBDIR);
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
  else
    strcat (command, command1);
  /* compilation */
  if (!dep && !tags && !source) {
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
  return 0;
}
