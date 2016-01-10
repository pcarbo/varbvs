% NDSUM(X,DIMS) returns NDOP(@SUM,X,DIMS).
function x = ndsum (x, dims)
  x = ndop(@sum,x,dims);
