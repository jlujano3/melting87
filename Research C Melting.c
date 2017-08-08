/* This is John Lujano's research with Dr. Tausch for the paper entitled
An Optimization Method for Moving Interface Problems Goverend by the Heat Equation*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define ONETHRD 0.333333333333333333333333333333
#define M_SQPI 1.77245385090552
#define M_SQPINV 0.564189583547756287

/*
 * legendre dynamic algorithm
 * Parameters 
 *    n    degree of legendre polynomial
 *    x    x value
 * Return
 *    y    return values
 */

void leg(int n, double x, double *y){
  int i;
  
  y[0]=1.0;
  if (n==0) return;
  y[1]=x;
  if (n==1) return;
  for (i=2;i<=n;i++){
    y[i]=((2.0*((double)i-1.0)+1.0)*x*y[i-1]-((double)i-1.0)*y[i-2])/(double)i;
  }

}


/*
 * derivative of legendre dynamic algorithm
 * Parameters 
 *    n    degree of legendre polynomial
 *    x    x value
 *    y    legendre values
 * Return
 *    z    return values
 */

void dleg(int n, double x,double *y,double *z){
  int i,j;
  
  z[0]=0.0;
  if (n==0) return;
  z[1]=1.0;
  if (n==1) return;
  for (i=2;i<=n;i++){
    z[i]=0.0;
    for (j=i-1;j>=0;j--){
      if ((i-j)%2==0){
        z[i]=z[i];
      }
      else{
        z[i]+=(2.0*((double)j)+1.0)*y[j];
      }
    }
  }
}

/*
 * Compute the radius at a given time t
 * Parameters 
 *    t    time
 * Return
 *    g    radius
 */

double rKnown(double t){
  
  double g;
  
  g=t;
  //g=2*sqrt(t)*alpha;
  
  return g;
}

/*
 * Compute the derivative of the radius with respect to time
 * Parameters 
 *    t    time
 * Return
 *    g    dr/dt
 */

double dRKnown(double t){
  
  double g;
  
  g=1.0;
  //g=alpha/sqrt(t);
  
  return g;
}

/*
 * Give the value for a test u(x,t)
 * Parameters 
 *    x    the radius at a given time
 *    t    the time
 * Return
 *    g    u
 */

double u(double x,double t){
  
  double g;
  
  g=exp(t-x)-1.0;
  //g=erf(alpha)-exp(x/(2.0*sqrt(t)));
  
  return g;
}

/*
 * Give the value for a test du(x,t)
 * Parameters 
 *    x    the radius at a given time
 *    t    the time
 * Return
 *    g    du
 */

double du(double x,double t){

  double g;

  g = -1.0*(exp(t-x));
  //g = exp(-x*x/(4.0*t))/sqrt(M_PI*t);

  return g;
}

/*
 * Give the value for a test uo(x,t)
 * Parameters 
 *    x    the radius at a given time
 * Return
 *    g    uo
 */

double uo(double x){
  
  double g;
  
  g=0;
  
  return g;
}

/*
 * Give the value for a test duo(x,t)
 * Parameters 
 *    x    the radius at a given time
 * Return
 *    g    duo
 */

double duo(double x){
  
  double g;
  
  g=0;
  
  return g;
}

/*
 * Give the value for a test forcing function force(t)
 * Parameters 
 *    t    the time
 * Return
 *    f    f
 */

double force(double t){

  double g;
  
  g=exp(t);
  //g=1/sqrt(M_PI*t);

  return g;
}

/*
 * Give the value for a test q(t)
 * Parameters 
 *    t    the time
 * Return
 *    q    q
 */

double qKnown(double t){

  double r,r1,q;
  
  r  = rKnown(t);
  r1 = dRKnown(t);
  q = du(r, t) + u(r, t)*r1;

  return q;
}

/*
 * Give the value for kv-> part of the green's function
 * Parameters 
 *    x    radius at the given timestep
 *    xtau the radius spanning up to the given timestep
 *    t    the time at our timestep
 *    tau  the time spanning up to the given timestep
 * Return
 *    g    duo
 */

double kv(double x, double xtau,double t, double tau){
  double g;

  g=0.5*M_SQPINV*(exp(-(x-xtau)*(x-xtau)/(4.0*(t-tau)))+exp(-(x+xtau)*(x+xtau)/(4.0*(t-tau))));
  
  return g;
}

/*
 * Give the value for kv at the current timestep which approaches the value
 * Parameters 
 * Return
 *    g    kvtntn
 */

double kvtntn(double x,double t){
  double g;
  
  //g=0.5*M_SQPINV*(1.0-exp(-(x*x)/(2*t)));
  g=1.0*M_SQPINV;
  
  return g;
}

/*
 * Algorithm for finding the single layer potential in the integral equation
 * Parameters 
 *    n    timestep index
 *    h    dt
 *    r    radius time history
 *    q    history of q values
 * Return
 *         final value
 */

double sLay(int n,double h,double *r,double *q){
  double tn=n*h;
  int i;
  double sum=0.5*kv(r[n],r[0],tn,0.0)*q[0]/sqrt(tn-0.0);
  
  for (i=1; i<n; i++){
    sum += kv(r[n],r[i],tn,i*h)*q[i]/sqrt(tn-i*h);
  }
  
  return sum*h;
}

/*
 * Algorithm for finding the single layer potential in the integral equation evaluated as SG(0,rtau,t,tau)f(t)dt
 * Parameters 
 *    n    timestep index
 *    h    dt
 *    r    radius time history
 *    q    history of q values
 * Return
 *         final value
 */

double zeroLay(int n,double h,double *r,double *q){
  double tn=n*h;
  double t,tau,xtau,x=0.0;
  int i;
  double sum=0.5*kv(0.0,r[0],tn,0.0)*q[0]/sqrt(tn-0.0);
  
  for (i=1; i<n; i++){
    sum += kv(0.0,r[i],tn,i*h)*q[i]/sqrt(tn-i*h);
    xtau=r[i];
    t=tn;
    tau=i*h;
    //printf("kv(0.0,r[i],tn,i*h) %.9f\n",kv(0.0,r[i],tn,i*h));
    //printf("exp(-(x-xtau)*(x-xtau)/(4.0*(t-tau)))->%f\n",exp(-(x-xtau)*(x-xtau)/(4.0*(t-tau))));
    //printf("exp(-(x+xtau)*(x+xtau)/(4.0*(t+tau)))->%f\n",exp(-(x+xtau)*(x+xtau)/(4.0*(t+tau))));
    //printf("x-> %f xtau-> %f \n",x,xtau);
    //printf("t-> %f tau-> %f \n",t,tau);
  }
  return sum*h;
}

/*
 * Algorithm for finding the mu value in the discretization
 * Parameters 
 *    n    timestep index
 *    h    dt
 * Return
 *         final value
 */

double mun(int n,double h){
  double tn=n*h;
  int i;
  double sum=0.5/sqrt(tn-0.0);
  
  for (i=1; i<n; i++){
    sum += 1.0/sqrt(tn-i*h);
  }
  
  return 2.0*sqrt(tn)-sum*h;
}

/*
 * Algorithm for finding the single layer potential in the adjoint integral equation
 * Parameters 
 *    n    timestep index
 *    h    dt
 *    r    radius time history
 *    q    history of q values
 * Return
 *         final value
 */

double sLayAdj(int n, double h, double *r, double *qt){
  double ktntjsum=0.0 , tn=n*h;
  int i;
  
  for (i=1;i<n;i++){
    ktntjsum+=kv(r[n],r[i],tn,i*h)*qt[i]/(sqrt(tn-i*h));
  }
  
  return h*ktntjsum;
}

/*
 * Algorithm for finding the coefficient of q[0] in the adjoint
 * Parameters 
 *    n    timestep index
 *    h    dt
 * Return
 *         final value
 */

double qto(int n, double h){
  double sum=0.0;
  int i;

  for (i=1;i<n;i++){
    sum-=sqrt((double)(n)/(double)(i)-1.0);
  }


  return M_PI*0.5+sum/((double)n);
}

/*
 * Algorithm for finding the mu value in the adjoint integral equation
 * Parameters 
 *    n    timestep index
 *    h    dt
 * Return
 *         final value
 */

double munAdj(int n, double h){
  double tn=n*h , sum=0.0;
  int i;
  
  for (i=1;i<n;i++){
    sum += sqrt(((double)i)/((double)(n-i)));
  }
  
  sum = 0.5*M_PI - sum/(double)n;
  return sum*sqrt(tn);      //jtdeb: was division
}

/*
 * Algorithm for solving the state equation
 * Parameters 
 *    n    number of intervals
 *    tmax maximum times
 *    r    history of radius values
 * Return
 *    qs   the state q values
 */

void solveState(int n, double tmax, double *r,double *qs,double *a,int deg, double *dr){
  int i,j,k,l;
  double *us,*fs;
  double duoVal,dt,totmax,totmaxinv,xk,tk;
  
  us=(double*)malloc( (n+1)*sizeof(double) );
  fs=(double*)malloc( (n+1)*sizeof(double) );
  
  dt=tmax/(double)(n);
  
  for (i=0;i<=n;i++){
    //us[i]=u(r[i],i*dt);
    us[i]=0.0;
    fs[i]=force(i*dt);
  }
  
  duoVal=du(r[0],0);
  //duoVal=duo(r[0]);
  
  qs[0]=duoVal+u(r[0],0)*dr[0];
  
  for (j=1;j<=n;j++){
    qs[j]=-1.0*(sLay(j,dt,r,qs)+zeroLay(j,dt,r,fs)+mun(j,dt)*fs[j])/(mun(j,dt)*kvtntn(r[j],j*dt));
  }
  
  free(us);
}

/*
 * Algorithm for solving the adjoint equation
 * Parameters 
 *    n    timestep index
 *    tmax maximum times
 *    r    history of radius values
 *    qs   state q values
 * Return
 *    qa   the adjoint q values
 */

void solveAdj(int n, double tmax, double *r ,double *qs,double *qarev,double *a, int deg, double *dr){
  int h,i , j,k,l,m;
  double *us, *rrev, *qa,*fsrev;
  double dt,totmax,totmaxinv,xk,tk;
  
  qa=(double*)malloc( (n+1)*sizeof(double) );
  us=(double*)malloc( (n+1)*sizeof(double) );
  rrev=(double*)malloc( (n+1)*sizeof(double) );
  fsrev=(double*)malloc( (n+1)*sizeof(double) );
  
  dt=tmax/(double)(n);
  
  for (i=0;i<=n;i++){
    //us[i]=u(r[i],i*dt);
    us[i]=qs[n-i];
    rrev[i]=r[n-i];
    fsrev[i]=force(((double)n-(double)i)*dt);
  }  

  qa[0]=-1.0*M_SQPINV*qs[n];
  
  for (j=1;j<=n;j++){
    qa[j]=-1.0*(sLay(j,dt,rrev,us)+zeroLay(j,dt,rrev,fsrev)+mun(j,dt)*fsrev[j])/(mun(j,dt)*kvtntn(rrev[j],j*dt));
  }
  
  for (m=0;m<=n;m++){
    qarev[m]=qa[n-m];
    //qarev[m]=qa[m];
  }
  
  free(qa),free(rrev),free(us);
}

/*
 * compute the boundary curve based on the Lagrange coefficients
 * Parameters 
 *    N    degree of Legendre polynomial
 *    M    number of time steps
 *    a    Legendre Polynomial Coefficients
 *    tmax Maximum Time
 * Return
 *    r    The front location as a function of time
 */

void bdry_curve(int N, int M, double *a, double *r, double tmax){
  double *y;
  double dt, xi, ti;
  int i, j;
  
  dt=tmax/(double)M;
  
  y=(double*)calloc(N+1, sizeof(double));
  
  for (i=0;i<M+1;i++){
    ti=i*dt;
    xi=2.0/tmax*(ti-tmax/2.0);
    leg(N,xi,y);
    r[i]=0.0;
    for (j=0;j<N+1;j++){
      r[i]+=a[j]*y[j];
    }
    r[i]*=sqrt(ti);
  }
  free(y);
}


/*
 * compute the derivative of the boundary curve based on the Lagrange coefficients
 * Parameters 
 *    deg  degree of Legendre polynomial
 *    n    number of time steps
 *    a    Legendre Polynomial Coefficients
 *    tmax Maximum Time
 * Return
 *    dr  The front location as a function of time
 */

void d_bdry_curve(int deg, int n, double *a, double *dr, double tmax){
  double *y,*z;
  double dt, xk, tk,totmax,totmaxinv;
  int k, l;
  
  y=(double*)calloc(deg+1, sizeof(double));
  z=(double*)calloc(deg+1, sizeof(double));
  
  totmax=2.0/tmax;
  totmaxinv=1.0/totmax;
  dt=tmax/(double)n;
  
  for (k=0;k<=n;k++){
    tk=k*dt;
    xk=totmax*(tk-totmaxinv);
    leg(deg,xk,y);
    dleg(deg,xk,y,z);
    dr[k]=0.0;
    for (l=0;l<=deg;l++){
      dr[k]+=0.5*a[l]*y[l]/sqrt(tk)+sqrt(tk)*totmax*a[l]*z[l];
    }
  }
  free(y);
  free(z);
}


/*
 * returns the functional and the modified Neumann data on the fixed boundary
 *
 * Parameters
 *   M     number of time steps
 *   tMax  maximal time
 *   r     boundary curve
 *   q     modified Neumann data on the boundary curve (return)
 * Return
 *   fun   Value of the error functional
 */

double compute_functional(int M, double tMax, double *q){
  double fun,dt;
  int i;
  
  dt=tMax/(double)(M);
  
  fun=(0.5*q[0]*q[0]+0.5*q[M]*q[M]);
  
  for (i=1;i<M;i++){
    fun+=q[i]*q[i];
  }
  
  return 0.5*dt*fun;
}

/*
 * compute the shape gradient using the adjoint method.
 *
 * Parameters:
 *   N     degree of Legendre polynomial        
 *   M     number of time steps
 *   a     legendre polynomial coefficients
 *   tMax  maximal time
 *   qS    modified Neumann data computed by compute_functional()
 *   dr    r'
 *   a2    a values defined by t*pi(z(t))
 *   
 * Return:
 *   g     gradient
 */

void compute_gradient(int N, int M, double tMax, double *a, double *qS, double *r,double *g,double *dr){
  
  double *qa,*y,*z;
  double dt,tk,xk,tj,xj,totmax,totmaxinv, hatr, dhatr;
  int i,j,k,l,m,n;
  
  y=(double*)calloc(N+1, sizeof(double));
  z=(double*)calloc(N+1, sizeof(double));
  qa=(double*)calloc(M+1, sizeof(double));
  
  solveAdj(M,tMax,r,qS,qa,a,N,dr);
  
  dt=tMax/(double)M;
  totmax=2.0/tMax;
  totmaxinv=1.0/totmax;
  
  for (n=0;n<N+1;n++){
    g[n]=0.0;
  }

  for (j=0;j<M+1;j++){
    tj=j*dt;
    xj=totmax*(tj-totmaxinv);
    leg(N,xj,y);
    dleg(N,xj,y,z);
    for (i=0;i<N+1;i++){
      hatr = tj*y[i];
      dhatr = totmax*tj*z[i]+y[i];
      if (j==0){
        g[i] += 0.5*(qS[j]*dhatr - (qa[j]+qS[j]*dr[j])*hatr*(qS[j]-dr[j]));
      }
      else if(j==M){
        g[i] += 0.5*(qS[j]*dhatr - qS[j]*dr[j]*hatr*(qS[j]-dr[j]));
      }
      else{
        g[i] += qS[j]*dhatr - (qa[j]+qS[j]*dr[j])*hatr*(qS[j]-dr[j]);
      }
    }
  }
  
  for (m=0;m<N+1;m++){
    hatr=tj*y[m];
    g[m]=g[m]*dt-mun(M,dt)*qa[M]*(qS[M]-dr[M])*hatr;
  }
  free(y),free(z),free(qa);
}

/*
 * compute the shape gradient using the adjoint method.
 *
 * Parameters:
 *   N     degree of Legendre polynomial        
 *   M     number of time steps
 *   tMax  maximal time
 *   
 * Return:
 *   NONE
 */

void test_gradient(int N, int M, double tMax) {
  int l, j,i;
  double fcnl, fcnl1, fcnl2,Gd;
  double *q, *a, *r, *g, *dr;
  double eps = 1e-2;

  q  = (double*)calloc(M+1,sizeof(double));
  a  = (double*)calloc(N+1,sizeof(double));
  g  = (double*)calloc(N+1,sizeof(double));
  r  = (double*)calloc(M+1,sizeof(double));
  dr = (double*)calloc(M+1,sizeof(double));

  /* test curve and exact gradient */
  for (j=0; j<=N; j++ ){
    a[j]=0.1;
  }
  
  bdry_curve(N, M, a, r, tMax);
  d_bdry_curve(N, M, a, dr, tMax);
  solveState(M,tMax,r,q,a,N,dr);
  fcnl = compute_functional(M, tMax, q);
  compute_gradient(N, M, tMax, a, q, r,g,dr);

  for (l=0; l<=N; l++) {
    /* first approximation using eps/2 */
    a[l] += eps;
    bdry_curve(N, M, a, r, tMax);
    d_bdry_curve(N, M, a, dr, tMax);
    solveState(M,tMax,r,q,a,N,dr);
    a[l] -= eps;
    fcnl1 = compute_functional(M, tMax, q);
    
    a[l] -= eps;
    bdry_curve(N, M, a, r, tMax);
    d_bdry_curve(N, M, a, dr, tMax);
    solveState(M,tMax,r,q,a,N,dr);
    a[l] += eps;
    fcnl2 = compute_functional(M, tMax, q);
    
    Gd = (fcnl1 - fcnl2)/(2.0*eps);
    printf("%d grad=%lf fd=%lf ", l, g[l], Gd);

    /* second approximation using eps/4 */
    a[l] += 0.25*eps;
    bdry_curve(N, M, a, r, tMax);
    d_bdry_curve(N, M, a, dr, tMax);
    solveState(M,tMax,r,q,a,N,dr);    
    a[l] -= 0.25*eps;
    fcnl1 = compute_functional(M, tMax, q);
    
    a[l] -= 0.25*eps;
    bdry_curve(N, M, a, r, tMax);
    d_bdry_curve(N, M, a, dr, tMax);
    solveState(M,tMax,r,q,a,N,dr);
    a[l] += 0.25*eps;
    fcnl2 = compute_functional(M, tMax, q);
    
    Gd = 4.0*(fcnl1 - fcnl2)/(2.0*eps);
    printf("fd=%lf ", Gd);
    
    /* third approximation using eps/8 */
    a[l] += 0.125*eps;
    bdry_curve(N, M, a, r, tMax);
    d_bdry_curve(N, M, a, dr, tMax);
    solveState(M,tMax,r,q,a,N,dr);
    a[l] -= 0.125*eps;
    fcnl1 = compute_functional(M, tMax, q);
    
    a[l] -= 0.125*eps;
    bdry_curve(N, M, a, r, tMax);
    d_bdry_curve(N, M, a, dr, tMax);
    solveState(M,tMax,r,q,a,N,dr);
    a[l] += 0.125*eps;
    fcnl2 = compute_functional(M, tMax, q);
    
    Gd = 8.0*(fcnl1 - fcnl2)/(2.0*eps);
    printf("fd=%lf  ", Gd);
    printf("diff = %lf\n", Gd-g[l]);
  }
  exit(1);

}

/*
 * determines the Quasi-Newton-direction
 *
 * Parameters:
 *   a     edge coefficients        
 *   b     gradients
 *   x     direction to be used
 *   m     number of iterations
 *   n     number of unknowns
 *   
 * Return:
 *   y     direction determined
 */

double *quasi_newton(double *a,double *b,double *x,unsigned int m,unsigned int n){
unsigned int	i;
double		c, d, e, *y;

if (m == 0)
{  y = (double*) malloc(n*sizeof(double));
   memcpy(y,x,n*sizeof(double));
   }
else
{  c = d = 0;
   for (i=0; i<n; i++) {  
      e = a[m*n+i]-a[(m-1)*n+i];
      c += e*(b[m*n+i]-b[(m-1)*n+i]);
      d += e*x[i];
      }
   if (c == 0) 	/* Then this was just a gradient step */{
      y = quasi_newton(a,b,x,m-1,n);
      }
   else		/* Then this was a quasi-Newton step */{
      for (i=0; i<n; i++){
        x[i] -= (b[m*n+i]-b[(m-1)*n+i])*d/c;
      }
      e = 0;
      y = quasi_newton(a,b,x,m-1,n);
      for (i=0; i<n; i++){ 
        e += (b[m*n+i]-b[(m-1)*n+i])*y[i];
        
      }
      for (i=0; i<n; i++){
        y[i] += (a[m*n+i]-a[(m-1)*n+i])*(d-e)/c;
      }
      }
   }
return(y);
}


void test_optimization(int numInts, int deg, int tmax, double h_0, int inner_itermax,  int line_search_itermax,double init_guess){
  double  *r;                      /* boundary curve              */
  double  *dr;                     /* derivative of boundary curve*/
  double  *a, *a_new;	             /* current domains             */
  double  *g;		                   /* gradient                    */
  double  *A, *G;		               /* list of domains, gradients  */
  double  df, f, f_new;	           /* values of functionals       */
  double  grad_norm;	             /* norm of gradients           */
  double  h;		                   /* stepsize                    */
  double  *d;		                   /* search direction            */
  double  *q;	                     /* modified Neumann data       */
  double  error;		               /* error of q                  */
  FILE	  *file;	   	             /* protocol file               */
  time_t  t1, t2, t3;              /* time variables              */
  int     i,j,k,l;                 /* indices                     */
  double  eps=1e-12;               /* difference                  */

  time(&t1);
  d = (double*)calloc( deg+1,sizeof(double) );
  a = (double*)calloc( deg+1,sizeof(double) );
  a_new = (double*)calloc( deg+1,sizeof(double) );
  g = (double*)calloc(deg+1,sizeof(double));
  A = (double*)calloc( inner_itermax*(deg+1),sizeof(double) );
  G = (double*)calloc( inner_itermax*(deg+1),sizeof(double) );
  q = (double*)calloc(numInts+1,sizeof(double*)); 
  r = (double*)calloc(numInts+1,sizeof(double*));
  dr = (double*)calloc(numInts+1,sizeof(double*)); 
  
  /* initial guess */
  a[0] = init_guess;
  bdry_curve(deg, numInts, a, r,tmax);
  d_bdry_curve(deg, numInts, a, dr, tmax);
  solveState(numInts,tmax,r,q,a,deg,dr);
  f_new = f = compute_functional(numInts, tmax, q);

  file = fopen("protocol.m","w+");
  if (file == NULL) {  
    printf("ERROR: cannot open file \n");
  }
  else {  
    k = inner_itermax+1;
    fprintf(file,"f = zeros(%d,1);\ng = zeros(%d,1);\na = zeros(%d,%d);\n\n", k,k,deg+1,k);
    fprintf(file,"f(1)=%g;\n",f);
    fprintf(file,"\na(:,1)=[");
    for (i=0; i<=deg; i++) fprintf(file,"\n%g",a[i]);
    fprintf(file,"\n];\n");
    fclose(file);
  }
  
    /* inner loop */
    for (j=0; j<inner_itermax; j++) {  
      /* compute the current gradient */
      time(&t2);
      bdry_curve(deg, numInts, a, r, tmax);
      d_bdry_curve(deg, numInts, a, dr, tmax);
      solveState(numInts,tmax,r,q,a,deg,dr);
      compute_gradient(deg, numInts, tmax, a, q, r,g,dr);
      for ( grad_norm=0, k=0; k<=deg; k++ ) {  
        grad_norm += g[k]*g[k];
      }
      grad_norm = sqrt(grad_norm); 
      for ( error=0, k=0; k<=deg; k++ ) {  
        error += q[k]*q[k];
      }
      error = sqrt(error);
       
      
      printf("Step %d | Functional %g | Norm of Gradient %g | Error q %g\n",j,f,grad_norm,error);
      //printf("%g\n",f);
      //printf("%g\n",grad_norm);
      //printf("%g\n",error);
	 
      /* Line search */
      printf("Line-search: ");
      /* determine the search direction and initial step size */
      memcpy(&A[j*(deg+1)],a,(deg+1)*sizeof(double));
      memcpy(&G[j*(deg+1)],g,(deg+1)*sizeof(double));
      d = quasi_newton(A,G,g,j,deg+1);
      if(j<10){
        h = 0;
        for (k=0; k<=deg; k++) h += d[k]*d[k];
        if (h > 1) {  
          h = sqrt(h);
          for (k=0; k<=deg; k++) d[k] /= h;
        }
      }
      else{
        h=h_0;
      }

      for (k=0; k<line_search_itermax; k++) {   
        
        /* boundary update */
        for (l=0; l<=deg+1; l++) a_new[l] = a[l]-h*d[l];
        bdry_curve(deg, numInts, a_new,r, tmax);
        d_bdry_curve(deg, numInts, a, dr, tmax);
        solveState(numInts,tmax,r,q,a,deg,dr);
        f_new = compute_functional(numInts,tmax,q);
        printf("%d->%g(%g) | ",k,h,f_new);
       
        /* if necessary new step length */
        if ((f_new >= f+eps) && (h != 0))
          {  df = 0;
            for (l=0; l<=deg; l++)  df += h*g[l]*d[l];
            if (df/(df+f_new-f) > 0) h *= 0.5*df/(df+f_new-f); else h = 0;
          }
        else break;
      }
      
      /*j>20 give the procedure time to store information*/
      /*if the optimization stops improving break the loop*/
      if (f_new>f&&j>20){
        break;
      }

      /* print stuff */
      time(&t3);
      printf("%g sec\n",difftime(t3,t2));
      file = fopen("protocol.m","a+");
      if (file == NULL) {  
        printf("ERROR: cannot open file \n");
      }
      else { 
        fprintf(file,"g(%d)=%g;\n\n",j+1,grad_norm);
        fprintf(file,"f(%d)=%g;\n",j+2,f_new);
        fprintf(file,"a(:,%d)=[",j+2);
        for (l=0; l<=deg; l++) fprintf(file,"\n%g",a_new[l]);
        fprintf(file,"\n];\n");
        fclose(file);
      }
         
      /* final boundary update */
      f = f_new;
      memcpy(a,a_new,(deg+1)*sizeof(double));
    } /* inner loop */

  time(&t2);	/* Endzeit */
  printf("OVERALL CPU-TIME: %g sec\n",difftime(t2,t1));
  exit(1);
}



int main(int nargs, char *argv[]){
  double tmax=1.0,dt,eps=1e-4,h_0=1.0,fcnl,fcnl1,fcnl2,Gd,init_guess=0.25;
  int i,j,k,x,l,numInts=10,verbose=0,deg=10,M,inner_itermax=100,outer_itermax=10,line_search_itermax=20,job;
  double *r,*y,*z,*a,*q,*qa,*g,*t,*dr,*fs;
  
  r = (double*)calloc( numInts+1,sizeof(double) );
  fs = (double*)calloc( numInts+1,sizeof(double) );
  dr= (double*)calloc( numInts+1,sizeof(double) );
  t = (double*)calloc( numInts+1,sizeof(double) );
  q = (double*)calloc( numInts+1,sizeof(double) );
  qa = (double*)calloc( numInts+1,sizeof(double) );
  g = (double*)calloc( deg+1,sizeof(double) );
  a = (double*)calloc( deg+1,sizeof(double) );
  y = (double*)calloc( deg+1,sizeof(double) );
  z = (double*)calloc( deg+1,sizeof(double) );

  /* parse the command line */
  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
      case 'v': verbose = atoi( argv[i]+3 );
        break;
      case 'T': tmax = atof( argv[i]+3 );
        break;
      case 'N': numInts = atoi( argv[i]+3 );
        break;
      case 'M': deg = atoi( argv[i]+3 );  
        break;        
      case 'j': job = atoi( argv[i]+3 );
        break;
      }

#if 0
  printf ("Enter a job\n ");
  printf ("0-test gradient\n ");
  printf ("1-test optimization\n ");
  scanf ("%d",&job);
  printf ("You entered: %d\n", job);
  
  printf ("Enter numInts\n ");
  scanf ("%d",&numInts);
  printf ("You entered: %d\n", numInts);
  
  printf ("Enter deg\n ");
  scanf ("%d",&deg);
  printf ("You entered: %d\n", deg);
  
if (job==0){
  test_gradient(deg, numInts, tmax);
}
else if (job==1){
  test_optimization(numInts,deg,tmax,h_0,inner_itermax,line_search_itermax,init_guess);
}
#endif

dt=tmax/(double)numInts;

for(j=0;j<=numInts;j++){
  r[j]=rKnown(j*dt);
  dr[j]=dRKnown(j*dt);
  q[j]=qKnown(j*dt);
  //fs[j]=force(j*dt);
  fs[j]=force(j*dt);
  zeroLay(j,dt,r,fs);
}



#if 0
bdry_curve(deg,numInts,a,r,tmax);
d_bdry_curve(deg, numints, a, dr, tMax);
solveState(numInts,tmax,r,q,a,deg,dr);
solveAdj(numInts,tmax,r,q,qa,a,deg,dr);
compute_gradient(deg,numInts,tmax,a,q,r,g,dr);
#endif


solveState(numInts,tmax,r,q,a,deg,dr);
for (l=0;l<=numInts;l++){
  printf("%f\n",q[l]);
  //printf("%f\n",qa[l]);
  //printf("sLay-> %f mun -> %f zeroLay -> %f force -> %f\n",sLay(l,dt,r,q),mun(l,dt),zeroLay(l,dt,r,fs),force(l*dt));
  //printf("err->%f q->%f qKnown-> %f\n",qKnown(l*dt)-q[l],q[l],qKnown(l*dt));
}


}