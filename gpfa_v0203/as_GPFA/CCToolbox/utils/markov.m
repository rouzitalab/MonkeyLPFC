function [Pi,A] = markov(seq,NumSymbols,norm)
%MARKOV  Fit markov model
%   [Pi,A] = markov(seq,NumSymbols) counts the initial-state transitions
%   in Pi and the state transitions in A. seq is cell array of length
%   N in which each seq{i} is a T_i length state vector. Pi and A are
%   appropriately normalized.
%
%   seq can also be a simple N-by-T matrix which gives N state vectors
%   each of length T.
%
%   [Pi,A] = markov(seq,NumSymbols,'nonorm') is as above but returns
%   raw counts instead of normalized probabilities.

% Scott J Gaffney   11 September 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'markov';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

norm = cexist('norm','norm');
m = NumSymbols;
A = zeros(m,m);
Pi = zeros(m,1);

if (iscell(seq))
  for i=1:length(seq)
    Pi(seq{i}(1)) = Pi(seq{i}(1)) + 1;
    for t=2:length(seq{i})
      A(seq{i}(t-1),seq{i}(t)) = A(seq{i}(t-1),seq{i}(t)) + 1;
    end
  end
else
  [n,T] = size(seq);
  for i=1:n
    Pi(seq(i,1)) = Pi(seq(i,1)) + 1;
    for t=2:T
      A(seq(i,t-1),seq(i,t)) = A(seq(i,t-1),seq(i,t)) + 1;
    end
  end
end

if (strcmp(norm,'norm'))
  Pi = Pi ./ sum(Pi);
  A  = A  ./ (sum(A,2)*ones(1,m));
end