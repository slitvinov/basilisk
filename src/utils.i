%{
  static void resolution (int n) {
    extern int N;
    N = n;
  }
  extern double interpolate (scalar v, double xp, double yp);
  static void interpolate1D (scalar v, double * x, double * val, int len) {
    int i;
    for (i = 0; i < len; i++)
      val[i] = interpolate (v, x[i], 0.);
  }
  static void interpolate2D (scalar v, double * x, double * y, double * val, 
			     int len1, int len2) {
    int i;
    for (i = 0; i < len1*len2; i++)
      val[i] = interpolate (v, x[i], y[i]);
  }
%}

%inline %{
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

%apply (double * IN_ARRAY1, int DIM1) {(double * x, int len1)};
%apply (double * ARGOUT_ARRAY1, int DIM1) {(double * val, int len)};
%inline %{
  void _interpolate1D (scalar v, double * x, int len1, double * val, int len) {
    interpolate1D (v, x, val, len1);
  }
%}

%apply (double * IN_ARRAY2, int DIM1, int DIM2) {
  (double * x, int len3, int len4),
  (double * y, int len5, int len6)
}
%apply (double * INPLACE_ARRAY2, int DIM1, int DIM2) {
  (double * val, int len1, int len2)
}
%inline %{
  void _interpolate2D (scalar v, 
                       double * x, int len3, int len4,
                       double * y, int len5, int len6,
                       double * val, int len1, int len2) {
     interpolate2D (v, x, y, val, len3, len4);
  }
%}
