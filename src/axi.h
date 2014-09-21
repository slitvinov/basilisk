#define AXI 1

event init (i = 0) {
  cm = new scalar;
  foreach()
    cm[] = y;
  cm[top] = dirichlet(y);
  cm[bottom] = dirichlet(y);

  fm = new face vector;
  foreach_face()
    fm.x[] = y + 1e-10;

  boundary ({cm, fm});
}
