%--------------------------------------------------------------------------
% varbvs.m: One-sentence summary of function goes here.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Overview of function goes here.
%
% USAGE:
%    Summary of usage goes here.
%
% INPUT ARGUMENTS:
% Description of input arguments goes here.
%
% OUTPUT ARGUMENTS:
% Description of output arguments goes here.
%
% DETAILS:
%    Detailed description of function goes here.
%
% LICENSE: GPL v3
%
% DATE: December 28, 2015
%
% AUTHORS:
%    List contributors here.
%
% REFERENCES:
%    List of references goes here.
%
% SEE ALSO:
%    List related functions here.
%
% EXAMPLES:
%    Give some examples here.
%
function y = varbvspredict (fit, X)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % (1) CHECK INPUTS
  % ----------------
