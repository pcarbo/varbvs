#ifndef INCLUDE_VARBVSR
#define INCLUDE_VARBVSR

#include <Rinternals.h>

// Declarations for .Call entry points.
SEXP diagsq_Call (SEXP Xp, SEXP ap, SEXP yp);

SEXP diagsqt_Call (SEXP Xp, SEXP ap, SEXP yp);
  
SEXP varbvsnormupdate_Call (SEXP Xp, SEXP sigmap, SEXP sap, SEXP logoddsp,
			    SEXP xyp, SEXP dp, SEXP alphap, SEXP mup,
			    SEXP Xrp, SEXP ip);

SEXP varbvsbinupdate_Call (SEXP Xp, SEXP sap, SEXP logoddsp, SEXP dp,
			   SEXP xdxp, SEXP xyp, SEXP xdp, SEXP alphap,
			   SEXP mup, SEXP Xrp, SEXP ip);

SEXP varbvsbinzupdate_Call (SEXP Xp, SEXP sap, SEXP logoddsp, SEXP dp,
			    SEXP xdxp, SEXP xyp, SEXP dzrp, SEXP alphap,
			    SEXP mup, SEXP Xrp, SEXP ip);

SEXP varbvs_varbvsnormupdate_rcpp (SEXP XSEXP, SEXP sigmaSEXP, SEXP saSEXP,
				   SEXP logoddsSEXP, SEXP xySEXP, SEXP dSEXP,
				   SEXP alphaSEXP, SEXP muSEXP, SEXP XrSEXP,
				   SEXP iSEXP);
 
#endif
