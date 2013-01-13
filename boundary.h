void boundary_h (Data * m, int n)
{
  /* update stencils */
  foreach (m)
    for (int k = -NS/2; k <= NS/2; k++)
      for (int l = -NS/2; l <= NS/2; l++)
	h(k,l) = hn(k,l);

  int i, j;
  for (i = 1; i <= n; i++) {
    /* symmetry */
    j = 1;
    h(0,-1) = h(0,0);
    j = n;
    h(0,1) = h(0,0);
  }
  for (j = 1; j <= n; j++) {
    /* symmetry */
    i = 1;
    h(-1,0) = h(0,0);
    i = n;
    h(1,0) = h(0,0);
  }
}

void boundary_b (Data * m, int n)
{
  /* update stencils */
  foreach (m)
    for (int k = -NS/2; k <= NS/2; k++)
      for (int l = -NS/2; l <= NS/2; l++)
	b(k,l) = M(i+k,j+l).b[NS/2][NS/2];

  int i, j;
  for (i = 1; i <= n; i++) {
    /* symmetry */
    j = 1;
    b(0,-1) = b(0,0);
    j = n;
    b(0,1) = b(0,0);
  }
  for (j = 1; j <= n; j++) {
    /* symmetry */
    i = 1;
    b(-1,0) = b(0,0);
    i = n;
    b(1,0) = b(0,0);
  }
}

void boundary_u (Data * m, int n)
{
  /* update stencils */
  foreach (m)
    for (int k = -NS/2; k <= NS/2; k++)
      for (int l = -NS/2; l <= NS/2; l++) {
	u(k,l) = un(k,l);
	v(k,l) = vn(k,l);
      }

  int i, j;
  for (i = 1; i <= n; i++) {
    /* solid walls */
    j = 1; v(0,0) = 0.;
    j = n; v(0,1) = 0.;
  }
  for (j = 1; j <= n; j++) {
    /* solid walls */
    i = 1; u(0,0) = 0.;
    i = n; u(1,0) = 0.;
  }
}
