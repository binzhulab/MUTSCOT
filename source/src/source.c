/* History: Dec21 2018 Initial coding
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define DEBUG 0

#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}
#define FUZZ 1e-8


/*
void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}
*/


/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

/* Function to allocate a double matrix */
static double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

/* Function to free a matrix */
static void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

/* Function to fill in a matrix from a vector (by column) */
static void fillMat(vec, nr, nc, out)
double *vec, **out;
int nr, nc;
{
  int i, j, col=0, ii;

  ii = 0;
  for (j=0; j<nc; j++) {
    for (i=0; i<nr; i++) {
      out[i][col] = vec[ii];
      ii++;
    }
    col++;
  }

} /* END: fillMat */

static double quadFormSq(vec, mat, n)
int n;
double *vec, **mat;
{
  int i, j;
  double ret, sum;

  ret = 0.0;
  for (i=0; i<n; i++) {
    sum = 0.0;
    for (j=0; j<n; j++) sum += mat[i][j]*vec[j];
    ret += vec[i]*sum;
  }

  return(ret);

} /* END: quadform */

static double getP1(TestStat, nsim, EP, wt, Z, nr, nc, invScore)
int nsim, nr, nc;
double TestStat, *EP, *wt, **invScore, **Z;
{
  int i, j, k, sumps, scrLen;
  double YmEP, Sc1, *Scr, TsN, Y, EPj, ret;

  scrLen = nc + 1;
  Scr    = dVec_alloc(scrLen, 0, 0.0);
  sumps  = 0;

  for (i=0; i<nsim; i++) {
    Sc1 = 0.0;
    for (j=0; j<scrLen; j++) Scr[j] = 0.0;

    for (j=0; j<nr; j++) {
      EPj  = EP[j];
      Y    = rpois(EPj);
      YmEP = Y - EPj; 
      Sc1 += wt[j]*YmEP;
      for (k=0; k<nc; k++) Scr[k+1] += Z[j][k]*YmEP;
    }
    Scr[0] = Sc1;
    TsN = quadFormSq(Scr, invScore, scrLen);
    if (TsN > TestStat + FUZZ) sumps++;
  }
  ret = ((double) sumps)/nsim;
  free(Scr); 

  return(ret);

} /* END: getP1 */

void C_getP(pTestStat, pnsim, wt, EP, Zvec, pnr, pnc, invScoreVec, ret_code, ret_p)
int *pnsim, *pnr, *pnc, *ret_code;
double *pTestStat, *wt, *EP, *invScoreVec, *Zvec, *ret_p;
{
  int nsim, nr, nc, scrLen;
  double TestStat, p, **invScore, **Z;

  *ret_code = 1;
  nsim      = *pnsim;
  nr        = *pnr;
  nc        = *pnc;
  TestStat  = *pTestStat;
  scrLen    = nc + 1;

  invScore  = dMat_alloc(scrLen, scrLen, 0, 0.0);
  Z         = dMat_alloc(nr, nc, 0, 0.0);

  fillMat(invScoreVec, scrLen, scrLen, invScore);
  fillMat(Zvec, nr, nc, Z);

  /* For random number generation */
  GetRNGstate();

  p = getP1(TestStat, nsim, EP, wt, Z, nr, nc, invScore);

  matrix_free((void **) invScore, scrLen);
  matrix_free((void **) Z, nr);

  PutRNGstate();  
  *ret_p    = p;
  *ret_code = 0;
  return;

} /* END: C_getP */

void C_getScoreMat(wt, Zvec, EP, pZnr, pZnc, ret_code, ret_vec)
int *pZnr, *pZnc, *ret_code;
double *wt, *EP, *Zvec, *ret_vec;
{
  /* Compute: t(WZ)%*%diag(EP, length(EP))%*%(WZ)
     Return matrix as a vector (symmetric) !!!
     Zvec must be passed in by column !!!
  */
  int n, nc, ncz, i, j, k;
  double **WZ, sum;

  *ret_code = -1;
  n         = *pZnr;
  ncz       = *pZnc;
  nc        = ncz + 1;
  WZ        = dMat_alloc(n, nc, 0, 0.0);

  /* Get the WZ matrix, first column is wt, other columns are Z = PatCov[, -1] */  
  for (i=0; i<n; i++) WZ[i][0] = wt[i];
  k = 0;
  for (j=1; j<nc; j++) {
    for (i=0; i<n; i++) {
      WZ[i][j] = Zvec[k];
      k++;
    }
  }

  /* Compute the score matrix */
  for (i=0; i<nc; i++) {
    for (j=i; j<nc; j++) {
      sum = 0.0;
      for (k=0; k<n; k++) sum += WZ[k][i]*EP[k]*WZ[k][j]; 
      ret_vec[i*nc + j] = sum;
      if (i != j) ret_vec[j*nc + i] = sum; 
    }
  } 

  matrix_free((void **) WZ, n); 
 
  *ret_code = 0;
  return;

} /* END: C_getScoreMat */

