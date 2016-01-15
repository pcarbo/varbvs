#include <R.h>
#include <Rinternals.h>
#include "misc.h"
#include "varbvs.h"

// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------
// This function is used to implement the R function varbvsnormupdate.
// It is called in R using the .Call interface.
SEXP varbvsnormupdate_Call (SEXP Xp, SEXP sigmap, SEXP sap, SEXP logoddsp, 
			    SEXP xyp, SEXP dp, SEXP alphap, SEXP mup, 
			    SEXP Xrp, SEXP ip) {

  // (1) GET INPUTS AND OUTPUTS
  // --------------------------
  //
  // TO DO: Try to represent input 'i' as an integer.
  //
  double* X       = REAL(Xp);        // Input X.
  double  sigma   = *REAL(sigmap);   // Input sigma.
  double  sa      = *REAL(sap);      // Input sa.
  double* logodds = REAL(logoddsp);  // Input logodds.
  double* xy      = REAL(xy);        // Input xy.
  double* d       = REAL(dp);        // Input d.
  double* alpha   = REAL(alphap);    // Input and output alpha.
  double* mu      = REAL(mup);       // Input and output mu.
  double* Xr      = REAL(Xr);        // Input and output Xr.
  double* i       = REAL(ip);        // Input i.

  // Get the number of samples (n), the number of variables (p), and
  // the number of coordinate ascent updates (m).
  SEXP     d = getAttrib(Xp,R_DimSymbol);
  R_xlen_t n = INTEGER(d)[0];
  R_xlen_t p = INTEGER(d)[1];
  R_xlen_t m = length(ip);

  // (2) CYCLE THROUGH COORDINATE ASCENT UPDATES
  // -------------------------------------------
  for (R_xlen_t j = 0; j < m; j++) {
    R_xlen_t k = (R_xlen_t) i[j];

    // Get the kth column of matrix X.
    const double* xk = getColumn(X,k,n);

    // Perform the update.
    varbvsupdate(xk,xy[k],d[k],sigma,sa,logodds[k],alpha + k,mu + k,Xr,n);
  }
}

// ---------------------------------------------------------------------
// This function is used to implement the R function varbvsbinupdate.
// It is called in R using the .Call interface.
void varbvsbinupdate_Call (SEXP Xp, SEXP sap, SEXP logoddsp, SEXP up, 
			   SEXP xyp, SEXP xup, SEXP dp, SEXP alphap, SEXP mup, 
			   SEXP Xrp, SEXP ip) {

  // (1) GET INPUTS AND OUTPUTS
  // --------------------------
  double* X       = REAL(Xp);        // Input X.
  double  sa      = *REAL(sap);      // Input sa.
  double* logodds = REAL(logoddsp);  // Input logodds.
  double* u       = REAL(up);        // Input u.
  double* xy      = REAL(xyp);       // Input xy.
  double* xu      = REAL(xup);       // Input xu.
  double* d       = REAL(dp);        // Input d.
  double* alpha   = REAL(alphap);    // Input and output alpha.
  double* mu      = REAL(mup);       // Input and output mu.
  double* Xr      = REAL(Xrp);       // Input and output Xr.
  double* i       = REAL(ip);        // Input i.

  // Get the number of samples (n), the number of variables (p), and
  // the number of coordinate ascent updates (m).
  SEXP     d = getAttrib(Xp,R_DimSymbol);
  R_xlen_t n = INTEGER(d)[0];
  R_xlen_t p = INTEGER(d)[1];
  R_xlen_t m = length(ip);

  // (2) CYCLE THROUGH COORDINATE ASCENT UPDATES
  // -------------------------------------------
  for (R_xlen_t j = 0; j < m; j++) {
    R_xlen_t k = (R_xlen_t) i[j];

    // Get the kth column of matrix X.
    const double* xk = getColumn(X,k,n);

    // Perform the update.
    varbvsbinupdate(xk,xy[k],xu[k],d[k],u,sa,logodds[k],alpha + k,mu + k,Xr,n);
  }
}
