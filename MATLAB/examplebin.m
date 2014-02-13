% In this small example, I simulate a data set from an idealized genetic
% association study in which the genetic markers (single nucleotide
% polymorphisms, or SNPs) are independently distributed, and I assess the
% performance of the variational approximation for mapping associations
% between the genetic markers and the binary trait (e.g. case-control
% status) in this setting.
%
% 
clear

% SCRIPT PARAMETERS.
p = 1e3;  % The number of variables (SNPs).
n = 500;  % The number of samples.

% Set the random number generator seed.
seed = 1;
rng('state',seed);
randn('state',seed);

