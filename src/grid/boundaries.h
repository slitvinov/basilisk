// Generic boundaries

typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level)   (const Boundary * b, scalar * list, int l);
  // multigrid only
  void (* restriction) (const Boundary * b, scalar * list, int l);
  // quadtree only
  void (* halo_restriction)      (const Boundary * b, scalar * list, int l);
  void (* halo_restriction_flux) (const Boundary * b, vector * list);
  void (* halo_prolongation)     (const Boundary * b, scalar * list, 
				  int l, int depth);
};

static Boundary ** boundaries = NULL; // list of all boundaries

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = realloc (boundaries, sizeof(Boundary *)*(len + 2));
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      free (b);
  free (boundaries);
  boundaries = NULL;
}

#define boundary_iterate(type,...) {	     \
    Boundary ** _i = boundaries, * _b;	     \
    while ((_b = *_i++))		     \
      if (_b->type)			     \
	_b->type (_b, __VA_ARGS__);	     \
}

// Box boundaries

typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
