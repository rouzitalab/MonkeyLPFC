function [b,lhood,ops] = bic(trajs,ops)
%BIC  Run BIC on a set of curves.
%   [bic_val,lhood,Options] = bic(trajs,[Options])

% Scott J Gaffney   24 September 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'bic';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end



%% Begin Argument Processing
%
State.state  = rand('state');
State.nstate = randn('state');
%
%% End Argument Processing

% default options
ops = cexist('ops',[]);
ops = SetFieldDef(ops,'State',State);
ops = SetFieldDef(ops,'method','lrm');
ops = SetFieldDef(ops,'Krange',[2 3]);
ops = SetFieldDef(ops,'OrderRange',[3]);
ops = SetFieldDef(ops,'MsgHnd',[]);
ops = SetFieldDef(ops,'MsgPrefix','');

% Handle graphical message bar
CreatedMsgBar=0;
if (isempty(ops.MsgHnd))
  ops.MsgHnd = msgbar([],'');
  CreatedMsgBar=1;
end
MsgPrefix = ops.MsgPrefix;

% set the random states
rand('state',ops.State.state);
randn('state',ops.State.nstate);

% Begin BIC analysis
lenK  = length(ops.Krange);
lenP  = length(ops.OrderRange);
lhood = zeros(lenK,lenP);
b     = zeros(lenK,lenP);
for i=1:lenK
  ops.K = ops.Krange(i);
  for j=1:lenP
    ops.order = ops.OrderRange(j);
    ops.MsgPrefix = sprintf('%sOrder %d, K %d, ',MsgPrefix,ops.order,ops.K);
    model = curve_clust(trajs,ops);
    if (isnan(model.TrainLhood)),  b(i,j) = NaN;  continue;  end
    n = length(model.C);
    N = model.NumPoints;
    Nk = model.NumIndParams;
    
    lhood(i,j) = model.TrainLhood;
    b(i,j) = (lhood(i,j) - Nk/2*log(N))/N;
    lhood(i,j) = lhood(i,j)/N;
  end
end
ops.n = n;
ops.N = N;
ops.Nk = Nk;

if (CreatedMsgBar)
  delete(ops.MsgHnd);
  ops.MsgHnd = [];
end



