
// -----------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  varbvsupdate(X,xy.elems,d.elems,sigma,sa,logodds.elems,
	       alpha,mu,Xr,m,I.elems);
}

// Execute a single iteration of the coordinate ascent updates.
void varbvsupdate (const singlematrix& X, const double* xy, 
		   const double* d, double sigma, double sa, 
		   const double* logodds, doublevector& alpha, 
		   doublevector& mu, doublevector& Xr, mwSize m, 
		   const double* I) {
  // This is storage for columns of X.
  doublevector xk = newdoublevector(n);

  // Repeat for each coordinate ascent update.
  for (mwIndex j = 0; j < m; j++) {


  }

  // Free dynamically allocated memory.
  free(xk);
}
