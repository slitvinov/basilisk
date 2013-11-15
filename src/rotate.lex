%option noyywrap
%option yylineno
%{
  static int n = 0;

  static int identifier (int c) {
    return ((c >= 'a' && c <= 'z') || 
	    (c >= 'A' && c <= 'Z') || 
	    (c >= '0' && c <= '9') || c == '_');
  }

  static void rotate_string (char * string) {
    char s[] = "123", * j = string;
    int i = 0, len = strlen (string);
    for (int k = 0; k <= len + 1; k++, j++) {
      if (i < 3) {
	s[i++] = k <= len ? *j : '\0'; s[i] = '\0';
      }
      else {
	if ((s[0] == '.' || s[0] == '_') && 
	    s[1] >= 'x' && s[1] <= 'y' &&
	    !identifier(s[2]))
	  s[1] = 'x' + (s[1] + n - 'x') % 2;
	fputc (s[0], yyout);
	s[0] = s[1]; s[1] = s[2]; 
	s[2] = s[1] != '\0' ? *j : '\0';
      }
    }
    fputs (s, yyout);
  }

%}

ID     [a-zA-Z0-9_]
SP     [ \t]
WS     [ \t\v\n\f]
ES     (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))

%%

{ID}*_[xy] |
[.][xy]  rotate_string (yytext);

left   fputs (n % 2 ? "bottom" : "left",   yyout);
right  fputs (n % 2 ? "top"    : "right",  yyout);
top    fputs (n % 2 ? "right"  : "top",    yyout);
bottom fputs (n % 2 ? "left"   : "bottom", yyout);

{ID}+   ECHO;

val_{ID}*{WS}*\( |
(val|fine|coarse|allocated|neighbor){WS}*\( {
  int para = 1, dimension = 2;
  char * index[5];
  int len[5], i = 0, c = input();
  for (int j = 0; j < 5; j++) {
    index[j] = malloc (sizeof(char));
    index[j][0] = '\0';
    len[j] = 1;
  }
  while (para > 0 && c != ')' && c != EOF) {
    if (c == '(') para++;
    else if (c == ')') para--;
    if (c == ',' && para == 1 && i < 5)
      i++;
    else {
      index[i] = realloc (index[i], sizeof(char)*++len[i]);
      index[i][len[i]-2] = c;
      index[i][len[i]-1] = '\0';
    }
    c = input();
  }
  int start = 0;
  if (!strncmp (yytext, "val", 3) || 
      !strncmp (yytext, "fine", 4) || 
      !strncmp (yytext, "coarse", 6)) {
    if (i == dimension) {
      rotate_string (yytext);
      rotate_string (index[0]);
    }
    else {
      ECHO;
      fputs (index[0], yyout);
    }
    fputc (',', yyout);
    start = 1;
  }
  else
    ECHO;
  if (i == dimension - 1 + start)
    for (int j = 0; j < dimension; j++) {
      int k = (j + n) % dimension;
      rotate_string (index[k + start]);
      fputc (j < dimension - 1 ? ',' : ')', yyout);
    }
  else
    //  more than dimension indices: do not rotate
    for (int j = start; j <= i; j++) {
      fputs (index[j], yyout);
      fputc (j < i ? ',' : ')', yyout);
    }
  for (int j = 0; j < 5; j++)
    free (index[j]);
}

"//".*            { ECHO; /* consume //-comment */ }
.                   ECHO;
[\n]                ECHO;
({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+  ECHO; /* STRING_LITERAL */

%%

// rotate dimensions n times
int rotate (FILE * fin, FILE * fout, int nrotate)
{
  if (0) yyunput (0, NULL); // just prevents 'yyunput unused' compiler warning
  rewind (fin);
  yyin = fin;
  yyout = fout;
  yylineno = 1;
  n = nrotate;
  return yylex();
}

#if TEST
int main (int argc, char * argv[])
{
  rotate (stdin, stdout, atoi (argv[1]));
}
#endif
