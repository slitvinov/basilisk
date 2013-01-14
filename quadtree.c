/* Linear quadtree implementation based on 
 * 
 * K. Aizawa, K. Motomura, S. Kimura, R. Kadowaki and J. Fan
 * "Constant time neighbor finding in quadtrees"
 * ISCCSP 2008, Malta, 12-14 March 2008.
 * http://www.lcad.icmc.usp.br/~jbatista/procimg/quadtree_neighbours.pdf
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

int size (int l)
{
  return ((1 << (2*(l + 1))) - 1)/3;
}

int level (int p)
{
  p = 3*p + 1;
  int l = 0;
  while (p > 0) { p /= 4; l++; }
  return l - 1;
}

int code (int p, int r)
{
  int l = level (p);
  return (p - size (l - 1)) << 2*(r - l);
}

int index (int code, int l, int r)
{
  return size(l - 1) + (code >> 2*(r - l));
}

const char* byte_to_binary(int x )
{
  static char b[sizeof(int)*8+1] = {0};
  int y;
  long long z;
  for (z=1LL<<sizeof(int)*8-1,y=0; z>0; z>>=1,y++)
    {
      b[y] = ( ((x & z) == z) ? '1' : '0');
    }
  
  b[y] = 0;
  
  return b;
}

int x (int p, int r)
{
  int n = code (p, r), a = 0, m = 1, b = 1;
  for (int i = 0; i < 2*r - 1; i += 2, m *= 4, b *= 2)
    a += ((m & n) != 0)*b;
  return a;
}

int y (int p, int r)
{
  int n = code (p, r), a = 0, m = 2, b = 1;
  for (int i = 1; i < 2*r; i += 2, m *= 4, b *= 2)
    a += ((m & n) != 0)*b;
  return a;
}

int repeat (int a, int r)
{
  int s = 0;
  for (int i = 0; i < r; i++, a *= 4)
    s += a;
  return s;
}

static int left, right = 1, top = 2, bottom;

#define quad(n, d) ((((n|bottom)+(d&left))&left)|(((n|left)+(d&bottom))&bottom))

int neighbor (int p, int d, int r)
{
  int l = level (p);
  int n = code (p, r);
  d <<= (2*(r - l));
  return index (quad(n, d), l, r);
}

void * quadtree (int r, size_t s)
{
  int n = size (r);
  void * q = malloc (s*n);
  left = repeat (1, r);
  bottom = repeat (2, r);
  for (int p = size (r - 1); p < n; p++) {
    int q = neighbor (p, bottom | left, r);
    printf ("%d %d %d ( %d , %d ) | ( %d , %d ) ",
	    p, index (code (p, r), level (p), r), level (p),
	    x(p, r), y(p, r),
	    x(q, r), y(q, r));
    printf ("\n");
  }
  return q;
}

int main (int argc, char * argv[])
{
  void * q = quadtree (atoi(argv[1]), sizeof (double));
  return 0;
}
