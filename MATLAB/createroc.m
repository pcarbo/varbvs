% [FP,TP] = CREATEROC(Y,R) returns points (FP,TP) that make up the ROC
% curve; it calculates the number of false positives FP and the number of
% true positives TP at each threshold. Y is the vector of true binary
% labels, and R is a vector of the same length of predictions. Higher
% values of R correspond to stronger predictions that a label Y equals 1.
function [fp, tp] = createroc (y, r)

  % Sort the ratings from highest to lowest.
  [ans I] = sort(-r);
  y       = y(I);

  % Calculate the number of false positives and the number of true
  % positives for each decision threshold.
  tp = [0; cumsum(y);   sum(y)   ];
  fp = [0; cumsum(1-y); sum(1-y) ];
