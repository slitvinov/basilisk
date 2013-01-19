#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

#define DX    (L0/_n)
#define XC(i) ((i + 0.5)*DX + X0)
#define YC(j) ((j + 0.5)*DX + X0)
#define XU(i) ((i)*DX + X0)
#define YU(j) YC(j)
#define XV(i) XC(i)
#define YV(j) ((j)*DX + X0)

#define NOT_UNUSED(x) (x = x)

#define VARIABLES \
  double x  = XC(I), y  = YC(J); /* cell center */	\
  double xu = XU(I), yu = YU(J); /* staggered u-coordinates */ \
  double xv = XV(I), yv = YV(J); /* staggered v-coordinates */ \
  /* we need this to avoid compiler warnings */	\
  NOT_UNUSED(x);  NOT_UNUSED(y);		        \
  NOT_UNUSED(xu); NOT_UNUSED(yu);                 \
  NOT_UNUSED(xv); NOT_UNUSED(yv);

typedef struct _Data Data;

enum { right, left, top, bottom };

#define field(data,offset,type) (*((type *)(((char *) &(data)) + (offset))))
typedef size_t var;
#define var(a) offsetof(Data,a)
#define val(a,k,l) field(data(k,l), a, double)
