#define GRID "Multigrid"

#include <stdio.h>

#define GHOSTS 1        // number of ghost layers
#define I (i - GHOSTS)
#define J (j - GHOSTS)

size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*GHOSTS;
  return n*n;
}

size_t _totalsize (int l)
{
  size_t s = 0;
  while (l >= 0)
    s += _size(l--);
  return s;
}

int mlevel (int n)
{
  int r = 0;
  while (n > 1) { n /= 2; r++; }
  return r;
}

#define data(k,l) (((Data *)_m)[(i + k)*(_n + 2*GHOSTS) + (j + l)])

void * coarser_level (void * m, int n)
{
  return (void *) (((char *)m) + sizeof(Data)*(n + 2*GHOSTS)*(n + 2*GHOSTS));
}

void * finer_level (void * m, int n)
{
  return (void *) (((char *)m) - sizeof(Data)*4*(n + GHOSTS)*(n + GHOSTS));
}

#define foreach(m,n) {					\
  void * _m = m; int _n = n;				\
  for (int i = GHOSTS; i < _n + GHOSTS; i++)		\
    for (int j = GHOSTS; j < _n + GHOSTS; j++) {	\
      VARIABLES
#define end_foreach() }}

#define foreach_boundary(m,n,d) {		 \
  void * _m = m; int _n = n;			 \
  for (int _k = 1; _k <= _n; _k++) {		 \
    int i = d > left ? _k : d == right ? _n : 1; \
    int j = d < top  ? _k : d == top   ? _n : 1; \
    VARIABLES
#define end_foreach_boundary() }}

#define foreach_fine_to_coarse(m,n) {				\
  void * _m = m, * _mf = _m;					\
  _m = coarser_level (_m, n); int _n = (n)/2;			\
  for (int level = mlevel(_n);					\
       _n > 0;							\
       _mf = _m, _m = coarser_level (_m, _n), _n /= 2, level--)	\
    for (int i = GHOSTS; i < _n + GHOSTS; i++)			\
      for (int j = GHOSTS; j < _n + GHOSTS; j++)
#define end_foreach_fine_to_coarse() }

#define fine(a,k,l) field(CELL(_mf)[(2*i-GHOSTS+k)*2*(_n + GHOSTS) + (2*j-GHOSTS+l)].d, a, double)

void * init_grid (int n)
{
  int r = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (stderr, "multigrid: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    r++;
  }
  void * m = calloc (_totalsize(r), sizeof (Data));
  return m;
}

#define free_grid free
