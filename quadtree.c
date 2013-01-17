/* Linear quadtree implementation based on 
 * 
 * K. Aizawa, K. Motomura, S. Kimura, R. Kadowaki and J. Fan
 * "Constant time neighbor finding in quadtrees"
 * ISCCSP 2008, Malta, 12-14 March 2008.
 * http://www.lcad.icmc.usp.br/~jbatista/procimg/quadtree_neighbours.pdf
 *
 * This uses a Z-grid ordering
 */

#define GRID "Linear quadtree"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

/* the maximum level */
static int quad_r;

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

int code (int p, int l)
{
  return (p - size (l - 1)) << 2*(quad_r - l);
}

int index (int code, int l)
{
  return size(l - 1) + (code >> 2*(quad_r - l));
}

int quad_x (int p)
{
  int n = code (p, level(p)), a = 0, m = 1, b = 1;
  for (int i = 0; i < 2*quad_r - 1; i += 2, m *= 4, b *= 2)
    a += ((m & n) != 0)*b;
  return a;
}

int quad_y (int p)
{
  int n = code (p, level(p)), a = 0, m = 2, b = 1;
  for (int i = 1; i < 2*quad_r; i += 2, m *= 4, b *= 2)
    a += ((m & n) != 0)*b;
  return a;
}

int repeat (int a)
{
  int s = 0;
  for (int i = 0; i < quad_r; i++, a *= 4)
    s += a;
  return s;
}

static int left, right = 1, top = 2, bottom;

#define quad(n, d) ((((n|bottom)+(d&left))&left)|(((n|left)+(d&bottom))&bottom))

static int quad_id[3][3];

int quad_neighbor (int p, int i, int j)
{
  int d = quad_id[i+1][j+1];
  if (d == 0) return p;
  int l = level (p);
  int n = code (p, l);
  d <<= (2*(quad_r - l));
  return index (quad(n, d), l);
}

int quad_neighbor_finest (int p, int i, int j)
{
  int d = quad_id[i+1][j+1];
  if (d == 0) return p;
  int s = size (quad_r - 1);
  int n = p - s;
  return s + quad(n,d);
}

void * quadtree (int r, size_t s)
{
  quad_r = r;
  void * q = malloc (s*size (r));
  left = repeat (1);
  bottom = repeat (2);
  quad_id[0][2] = top|left;    quad_id[1][2] = top;    quad_id[2][2] = top|right;
  quad_id[0][1] = left;        quad_id[1][1] = 0;      quad_id[2][1] = right;
  quad_id[0][0] = bottom|left; quad_id[1][0] = bottom; quad_id[2][0] = bottom|right;
  return q;
}

#if 0
int main (int argc, char * argv[])
{
  int r = atoi(argv[1]);
  void * q = quadtree (r, sizeof (double));
  for (int p = size (r - 1); p < size (r); p++) {
    int q = neighbor (p, bottom | left, r);
    printf ("%d %d %d ( %d , %d ) | ( %d , %d ) ",
	    p, index (code (p, r), level (p), r), level (p),
	    x(p, r), y(p, r),
	    x(q, r), y(q, r));
    printf ("\n");
  }
  return 0;
}
#endif
