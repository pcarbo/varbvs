% DIAGPROD(A,B) efficiently computes DIAG(A*B').
function y = diagprod (A, B)
  y = double(sum(A.*B,2));
