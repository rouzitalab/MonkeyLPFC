function p = plot_splinecoef(B,T,X,ops)
%PLOT_SPLINECOEF  Plot lines representing spline regression coefficients
%   P = PLOT_SPLINECOEF(B,T,X,[OPS]) plots lines representing the
%   spline coeficients.
%
%   B        : spline coefficients (P-by-K-by-D)
%   T        : n-vector of evaluation points where plotting occurs (unless
%              D==2, in which case a 2D plot is generated)
%   X        : basis matrix
%
% Comments:
%  - plots 2D data against itself
%  - higher dimensional data is not separately handled
%  - uses colorplot; passes ops through
%  - defaults colorplot:PlotSymbols to off

% Scott J Gaffney   5 October 2001
% Department of Information and Computer Science
% University of California, Irvine.
%
% Changes
% -----------------


PROGNAME = 'plot_splinecoef';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

ops = cexist('ops',[]);
ops = SetFieldDef(ops,'PlotSymbols',0);

% Handle B and C
if (exist('B')~=1 | isempty(B))
  errorbox('Argument Error: B must be provided.',PROGNAME);
  return;
end
[P,K,D] = size(B);
B = permute(B,[1 3 2]);  % make PxDxK
C   = cexist('C', (1:K)');


%%%%%
%%% Perform the plotting
held = ishold;
p=[];
for k=1:K
  xb = X*B(:,:,k);
  if (size(xb,2)==2)  % Plot 2-Dimensional data against each other
    p{end+1} = colorplot(xb(:,1),xb(:,2),C(k),ops);
  else
    p{end+1} = colorplot(T,xb,C(k),ops);
  end
  hold on;
end

if (~held)
  hold off;
end
