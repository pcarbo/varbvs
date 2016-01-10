% NDOP(OP,X,DIMS) performs operation OP over all dimensions of X except
% dimensions DIMS. Input OP must be a function handle, and must take
% arguments of the form OP(X,I), where I is a dimension of X. All singleton
% dimensions are removed from the result using SQUEEZE. Note that if the
% result is a row vector, it is turned into a column vector.
function x = ndop (op, x, dims)
  
  % Get all the dimensions other than DIMS.
  dims = setdiff(1:ndims(x),dims);
  n    = length(dims);

  % Repeat for each dimension.
  for i = 1:n
    x = op(x,dims(i));
  end
  
  % Remove singleton dimensions from the final result. If the final
  % result is a row vector, turn it into a column vector.
  x = squeeze(x);
  if size(x,1) == 1
    x = x(:);
  end
