// Part of the varbvs package, https://github.com/pcarbo/varbvs
//
// Copyright (C) 2012-2018, Peter Carbonetto
//
// This program is free software: you can redistribute it under the
// terms of the GNU General Public License; either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANY; without even the implied warranty of
// MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
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

SEXP varbvsmixupdate_Call (SEXP Xp, SEXP sigmap, SEXP sap, SEXP wp, SEXP xyp,
			   SEXP dp, SEXP alphap, SEXP mup, SEXP Xrp, SEXP ip, 
			   SEXP epsp);
 
SEXP varbvs_varbvsnormupdate_rcpp (SEXP XSEXP, SEXP sigmaSEXP, SEXP saSEXP,
				   SEXP logoddsSEXP, SEXP xySEXP, SEXP dSEXP,
				   SEXP alphaSEXP, SEXP muSEXP, SEXP XrSEXP,
				   SEXP iSEXP);

#endif
