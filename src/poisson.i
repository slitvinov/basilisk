%{
  extern void mg_cycle (scalar * a, scalar * res, scalar * da,
			void (* relax) (scalar * da, scalar * res, 
				 int depth, void * data),
			void * data,
			int nrelax, int minlevel);
  
  extern int NITERMAX;
  extern double TOLERANCE;

  typedef struct {
    int i;              // number of iterations
    double resb, resa;  // maximum residual before and after the iterations
    double sum;         // sum of r.h.s.
  } mgstats;

  extern mgstats mg_solve (scalar * a, scalar * b,
			   double (* residual) (scalar * a, scalar * b, 
						scalar ** res,
						void * data),
			   void (* relax) (scalar * da, scalar * res, 
					   int depth, 
					   void * data),
			   void * data);
  
  struct Poisson {
    scalar a, b;
    face vector alpha;
    scalar lambda;
    double tolerance;
  };

  extern mgstats poisson (struct Poisson p);
  extern mgstats project (face vector u, scalar p, 
			  face vector alpha, double dt);
%}
