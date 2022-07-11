function [mccv,ops] = make_small(mccvout)
%MAKE_SMALL
%   [MCCV,Ops] = MAKE_SMALL(mccvout)

% Scott J. Gaffney  5 November 2003
% School of Information and Computer Science
% University of California, Irvine

PROGNAME = 'make_small';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

mccv=[];
ops=[];

methods = fieldnames(mccvout);
if (isempty(methods)),  return; end
NumMethods = length(methods);

% grab the options
ops = getfield(mccvout,methods{1});
ops = rmfield(ops,'runs');
if (isfield(ops,'RandPoints'))
  ops.RandPoints = ops.RandPoints(1);
end

% build the test score matrices
for i=1:NumMethods
  scores=[];
  runs = getfield(mccvout,methods{i},'runs');
  if (isfield(runs,'TestSSE'))
    scores = setfield(scores,'sse',[runs.TestSSE]);
  end
  if (isfield(runs,'TestLike'))
    scores = setfield(scores,'like',[runs.TestLike]);
  end
  if (isfield(runs,'bic'))
    scores = setfield(scores,'bic',[runs.bic]);
  end
  mccv = setfield(mccv,methods{i},scores);
end

