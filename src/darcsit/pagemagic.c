// Returns successfully if a file starts with /** or %{ or """ or :<<'DOC'

#include <stdio.h>
#include <string.h>

int main()
{
  int c;
  do
    c = getchar();
  while (c != EOF && strchr(" \t\v\n\f", c));
  if (c == '/') {
    c = getchar();
    if (c != '*')
      return 1;
    c = getchar();
    if (c != '*')
      return 1;
  }
  else if (c == '%') {
    c = getchar();
    if (c != '{')
      return 1;    
  }
  else if (c == '"') {
    c = getchar();
    if (c != '"')
      return 1;    
    c = getchar();
    if (c != '"')
      return 1;    
  }
  else if (c == '#') {
    do
      c = getchar();
    while (c != EOF && c != '\n');
    do
      c = getchar();
    while (c != EOF && strchr(" \t\v\n\f", c));
    char * s = ":<<'DOC'";
    while (*s) {
      if (c != *s)
	return 1;
      s++;
      if (*s)
	c = getchar();
    }
  }
  else
    return 1;  
  do
    c = getchar();
  while (c != EOF && strchr(" \t", c));
  return strchr("\v\n\f\r", c) == NULL;
}
