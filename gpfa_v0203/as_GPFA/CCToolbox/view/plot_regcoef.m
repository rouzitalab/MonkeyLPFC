function p = plot_regcoef(B,len,C,ops)
%PLOT_REGCOEF  Plot lines representing regression coefficients
%   P = PLOT_REGCOEF(B,[LEN],[C],[Ops]) plots lines representing the
%   regression coeficients.
%
%   B is assumed to take the form, [P,K,D] = size(B), where
%   P is one plus the order of the regression model, K is the number
%   of regressions run, or the number of classes of separate regressions,
%   and D is the number of output variables.
%
%   Important: If LEN is provided, then both Ops.X and Ops.T are
%   resized according to this. If LEN is not provided, then Ops.len will
%   be used instead. If Ops.len is not set, then LEN will be set
%   using, first, Ops.X if this is provided, second,
%   Ops.T if this is provided, and third, the default value if nothing 
%   else succeeds (see below for defaults).
%
%   Ops.plot_regcoef
%   ---------------------------------
%   .X   - basis matrix to be premultiplied to B in order to give y-values
%   .T   - n-vector of evaluation points where plotting occurs (unless
%          data is 2D, in which case it will use a 2D plot)
%   .len - this behaves exactly as the passed parameter LEN and is 
%          provided here for convience. This parameter is only used
%          if LEN is empty.
%
%   Defaults
%   ---------------------------------
%   B      - none
%   LEN    - 20  (used when .len, .X, and .T are all empty)
%   C      - 1:K
%   OPS
%    .plot_regcoef.X    - regmat(0:LEN-1,size(B,1))
%    .plot_regcoef.T    - 0:LEN-1
%    .plot_regcoef.len  - []
%    .SelectClass - 1:K
%   
%
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
% 14 January 2003 (SJG)
%  - Changed the time axis from 1:len to 0:len-1.
%
% 28 January 2003 (SJG)
%  - Revamped the whole function; it may not work correctly
%    with older sister-code.


PROGNAME = 'plot_regcoef';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

DEFAULT_LEN = 20;

ops = cexist('ops',[]);
ops = SetFieldDef(ops,'PlotSymbols',0);
ops = SetFieldDef(ops,'plot_regcoef',[]);
ops.plot_regcoef = SetFieldDef(ops.plot_regcoef,'X',[]);
ops.plot_regcoef = SetFieldDef(ops.plot_regcoef,'T',[]);
ops.plot_regcoef = SetFieldDef(ops.plot_regcoef,'len',[]);
len = cexist('len',ops.plot_regcoef.len);

% Handle B and C
if (exist('B')~=1 | isempty(B))
  errorbox('Argument Error: B must be provided.',PROGNAME);
  return;
end
[P,K,D] = size(B);
C   = cexist('C', (1:K)');
ops = SetFieldDef(ops,'SelectClass',1:K);

% Handle X, T, and len
X = ops.plot_regcoef.X;
T = ops.plot_regcoef.T(:);
if (~isempty(X))      % if X is provided
  if (~isempty(len))    % And if len is provided, then resize X according to len
    X = X(1:len,:);   
  else                  % else if len isn't provided, then set len according
    len = size(X,1);    % to X
  end
  if (~isempty(T))      % if T is also provided, then make it match size of X
    T = T(1:len);
  else                  % else create it using the size of X
    T = (0:len-1)';
  end
  
elseif (~isempty(T))  % else if X is not provided but T is
  if (~isempty(len))    % And if len is provided, then resize T using len
    T = T(1:len);
  else                  % else if len isn't provided, then set len using T
    len = length(T);
  end
  X = regmat(T,P-1);      % Now we can make X using T
  
else                  % otherwise if neither is provided
  len = DEFAULT_LEN;
  T = (0:len-1)';
  X = regmat(T,P-1);
end



%%%%%
%%% Perform the plotting
held = ishold;
p=[];
B = permute(B,[1 3 2]);  % make PxDxK
for k=1:K
  if (isempty(find(ops.SelectClass==k))), continue; end
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
