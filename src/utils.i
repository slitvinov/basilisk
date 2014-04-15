%{
  static void resolution (int n) {
    extern int N;
    N = n;
  }
  extern double interpolate (scalar v, double xp, double yp);

  typedef struct {
    double avg, rms, max, area;
  } norm;
  extern norm normf (scalar f);

  typedef struct {
    double min, max, sum, stddev, area;
  } stats;
  extern stats statsf (scalar f);

  extern void vorticity (const vector u, scalar omega);
%}

extern void resolution (int n);
extern double interpolate (scalar v, double xp, double yp = 0.);

typedef struct {
  double avg, rms, max, area;
} norm;
extern norm normf (scalar f);

typedef struct {
  double min, max, sum, stddev, area;
} stats;
extern stats statsf (scalar f);

extern void vorticity (const vector u, scalar omega);
