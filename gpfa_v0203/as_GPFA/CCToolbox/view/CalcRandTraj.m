function ops = CalcRandTraj(trajs,C,ops)
%CalcRandTraj  Set up options for random trajectory plotting
%   Ops = CalcRandTraj(Trajs,C,Ops)

% Scott J Gaffney   2 May 2005
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'CalcRandTraj';
if (~nargin)
  help(PROGNAME);
  return;
end

ops = cexist('ops',[]);
ops = SetFieldDef(ops,'SelectClass',[]);
ops = SetFieldDef(ops,'PlotNumRand',[]);

% Handle random trajectory plotting
if (~isempty(ops.PlotNumRand))
  if (~isempty(ops.SelectClass))
    cc = [];
    for i=1:length(ops.SelectClass)
      cc = [cc; find(C==ops.SelectClass(i))];
    end
    n = length(cc);
    rperm = randperm(n);
    if (ops.PlotNumRand > n), ops.PlotNumRand=n;  end
    ops.SelectTrajs = cc(sort(rperm(1:ops.PlotNumRand)));
  else
    n = trajslen(trajs);
    rperm = randperm(n);
    if (ops.PlotNumRand > n), ops.PlotNumRand=n;  end
    ops.SelectTrajs = sort(rperm(1:ops.PlotNumRand));
  end
end
