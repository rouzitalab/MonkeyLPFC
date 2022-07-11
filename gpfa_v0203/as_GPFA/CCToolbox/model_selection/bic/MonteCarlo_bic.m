function bicv = MonteCarlo_bic(trajs,Ops)
%MonteCarlo_bic  Run BIC on random subsets of a large dataset.
%   BICV = MonteCarlo_bic(Trajs,Ops)
%

% Scott J Gaffney   20 November 2003
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'MonteCarlo_bic';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


Ops        = cexist('Ops',[]);
Ops        = SetFieldDef(Ops,'State',[]);  
Ops.State = SetFieldDef(Ops.State,'state',rand('state'));  
Ops.State = SetFieldDef(Ops.State,'nstate',randn('state'));  
Ops.State = SetFieldDef(Ops.State,'sstate',Ops.State.state);  
Ops = SetFieldDef(Ops,'NumRuns',1);
Ops = SetFieldDef(Ops,'MsgHnd',[]);


% Calculate sizes of subsets
n = size(trajs,1);
Ops = SetFieldDef(Ops,'SubsetSize',n);
TrainLen  = Ops.SubsetSize;

% Calculate Random Data Splits
[trash,RandSplitState] = AdvRandState(Ops.State.sstate,n,Ops.NumRuns-1);

% Bring up Message Bar
CreatedMsgBar=0;
if (isempty(Ops.MsgHnd))
  Ops.MsgHnd = msgbar([],'');
  CreatedMsgBar=1;
end

%% Loop over runs
for (r=1:Ops.NumRuns)
  rand('state',RandSplitState(:,r));
  RandIndex = randperm(n);
  
  Ops.MsgPrefix = sprintf('Run %d, ',r);
  TrainTrajs = trajs(RandIndex(1:TrainLen),:);
  [Bic(:,:,r),lhood(:,:,r)] = bic(TrainTrajs,Ops);
end
bicv.bic = Bic;
bicv.lhood = lhood;
bicv.ops = Ops;

if (CreatedMsgBar)
  delete(Ops.MsgHnd);
end




