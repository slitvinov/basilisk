/**
# Discharge.h

We define several functions in order to impose a flux at boundaries
*/

/**
## Boundary structure

This structure stores all what we need in order to define the boundary
: a scalar "river" which is equal to zero or 1, it defines the river
location, two doubles “xb" and "yb” to store the x and y coordinates
of the boundary, the double “dx” equal to the maximum spatial
resolution (this can be optimised) and an integer kw which stores the
relative location of the boundary (top, bottom, left and right).
*/
typedef struct 
{
  scalar river;
  double xb;
  double yb;
  double dx;
  int kw;
} defb;


/**
## Real flux

This function return the exact incoming flow on the boundary for an
imposed water height
 */
static double realflux ( defb b, double eta_imp ){
  double q = 0, x = b.xb, y = b.yb, dtmax = 2000;
  scalar lim = b.river;
  for( int j = 0; j < N-1; j++ ){
    if( b.kw == bottom || b.kw == top ) x = b.xb + j * b.dx;
    else y = b.yb + j * b.dx;
    Point point = locate ( x, y );      
    if( lim[] != 0 ){
      double ub;
      if( b.kw == bottom || b.kw == top ) ub = b.kw == bottom ? u.y[] : - u.y[];
      else ub =  b.kw == left ? u.x[] : - u.x[];
      double hn = max ( eta_imp - zb[], 0 ); 
      if( h[] > dry || hn > dry ){ 
	double fh,fu;
	kurganov (hn, h[], ub, ub, b.dx, &fh, &fu, &dtmax);
	q += fh * b.dx;
      }
    }         
  }
  return q;
}

/**
## Approximate flux

The function appflux returns an approximate value of the flux without
calling the kurganov routine. This function was initially built in
order to save computation time, but with lost of precision. In fact,
depending on the topography, this function can be less linear than the
realflux function. So, the false position method is slower when
calling this function than when calling the realflux one.
*/
/*
static double appflux(defb b,double eta_imp){
  double q = 0, x = b.xb, y = b.yb; 
  scalar lim=b.river;
  for(int j=0;j<N-1;j++){
    if(b.kw==top || b.kw==bottom) x=b.xb+j*b.dx;
    else y=b.yb+j*b.dx;
    Point point = locate(x,y);
    if(lim[] != 0){
      double ub;
      if(b.kw<=1) ub= b.kw==0 ? u.y[] : - u.y[];
      else ub=  b.kw==2 ? u.x[] : - u.x[];
      double hn = max(eta_imp-zb[],0); 
      if(hn > dry) q+= hn*ub*b.dx;
    }         
  }
  return q;
}
*/
/**
## Structure definition

 */
static void defbound( int keyw, scalar lim, defb * b){
  b->dx = L0/((double)N);
  b->river = lim;
  b->kw = keyw;
  b->xb = X0 + 0.5 * b->dx;
  b->yb = Y0 + 0.5 * b->dx;
  if( keyw == right ) b->xb += L0 - b->dx;
  else if( keyw == top )  b->yb += L0 - b->dx;
  else if( keyw != left && keyw != bottom ) assert (false);
}

/**
## Altitude

The altitude function find the minimum value of the topography. It
stores this value and also the one of the water height.
 */
static void altitude(defb b,  double * zmin, double * h0){
  double ztemp = *zmin, htemp = 0, x = b.xb, y = b.yb;
  scalar lim = b.river;

#if _OPENMP
#pragma omp parallel for reduction(min:ztemp)
#endif

  for( int j = 0; j < N-1; j++ ){
    if( b.kw == top || b.kw == bottom ) x = b.xb + j * b.dx;
    else y = b.yb + j * b.dx;
    Point point = locate (x,y);
    if( zb[] < ztemp && lim[] != 0 ){
      ztemp = zb[];
      htemp = h[];
    }
  }
  *zmin = ztemp;
  *h0 = htemp;
}

/**
## False position method

Here, we find the water height corresponding to the imposed inflow by
the [false position
method](http://en.wikipedia.org/wiki/False_position_method)
 */
static double metfalsepos(defb bo, double h0, double zmin,double prec, double q_imp){
  double binf, bsup, qinf, newb, qsup;

  // inflow at n-1
  double newq = realflux ( bo, h0+zmin );
  double ratio = (newq-q_imp)/q_imp;
  if( ratio*ratio <= prec*prec ) return h0 + zmin;
  
  // Find sup & inf boundaries
  if( newq > q_imp ){
    bsup = h0 + zmin;
    qsup = newq;
    binf = zmin;
    qinf = 0;
  }
  else{
    // We must take care of the case h0=0
    if( h0 <= dry ) {
      h0 = 1;
      bsup = zmin;
    }
    else  bsup = h0+zmin;
    do{
      // here we find a "good" superior limit
      binf = bsup;
      bsup = bsup + h0;
      qsup = realflux ( bo, bsup );
    } while ( qsup <= q_imp );
    
    qinf = realflux ( bo, binf );
  }

  // The false position method
  do{
    double alpha = ( qsup - qinf ) / ( bsup - binf );
    newb = binf + ( q_imp - qinf )/alpha;
    newq = realflux ( bo, newb );
    if( newq > q_imp ){
      bsup = newb;
      qsup = newq;
    }
    else {
      binf = newb;
      qinf = newq;
    }
    ratio = (newq-q_imp)/q_imp;
  } while( ratio*ratio >= prec*prec );
  
  return newb;
}

/**
## The discharge routine

 */
double discharge (double q_imp, int keyw, scalar limit){  
  
  // Inialisation of the defb structure
  defb bo;
  defbound( keyw, limit, &bo);
  
  // Finding zmin et h
  double h0 = 0, zmin = 200000;
  altitude(bo,&zmin,&h0);
  
  // If imposed inflow <= 0
  if( q_imp <= 0 ) return zmin-0.1;

  // Return eta found by Newton met
  return  metfalsepos(bo,h0,zmin,0.001,q_imp);
}


/**
# Hydrograph

This function return the inflow by linearly interpolating the data.
 */
double hydrograph ( const char * name, double mult ){
  FILE * fp;
  if( ( fp = fopen ( name , "r" ) ) == NULL ) {
    fprintf (stderr,"cannot open hydrograph data file.\n name = %s ",name);
    assert (false);
  }
  double time = 0, timea, q = 0, qa, alpha = 0;
  // We read the data at each call => Can definitively be optimised !
  do{
    qa = q;
    timea = time;
    if ( fscanf(fp,"%lf \t %lf \n",&time,&q) == EOF) break;
    time *= mult;
    if( time - timea != 0) alpha = ( q - qa ) / ( time - timea);
    else alpha = 0;
  } while (time < t);
  fclose(fp);
  /*
    If the solver time is sup to the larger time of the data, return the last value   */
  if ( timea  >= t )  return q;
  // Else, return the interpolation
  else return alpha * ( t - timea ) + qa;
}
