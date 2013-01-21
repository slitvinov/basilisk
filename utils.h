#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

#define GHOSTS 1        // number of ghost layers

#define NOT_UNUSED(x) (x = x)

#define VARIABLES

typedef struct _Data Data;

enum { right, left, top, bottom };

#define field(data,offset,type) (*((type *)(((char *) &(data)) + (offset))))
typedef size_t var;
#define var(a) offsetof(Data,a)
#define val(a,k,l) field(data(k,l), a, double)

#define pi 3.14159265358979
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
