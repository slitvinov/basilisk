// Returns successfully if a file starts with /**

#include <stdio.h>
#include <string.h>

int main()
{
  int c;
  do
    c = getchar();
  while (c != EOF && strchr(" \t\v\n\f", c));
  if (c != '/')
    return 1;
  c = getchar();
  if (c != '*')
    return 1;
  c = getchar();
  if (c != '*')
    return 1;
  do
    c = getchar();
  while (c != EOF && strchr(" \t", c));
  return strchr("\v\n\f\r", c) == NULL;
}
