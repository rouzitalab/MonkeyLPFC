function [sse,Ops,Lhood] = ScoreModels(mods,trajs,Ops)
%ScoreModels  Score a set of trained models (in MCCV or CV format)
%
%   sse = ScoreModels(Models,Trajs,[Options])
%
%   [sse,Ops] = ScoreModels(Models,Trajs,[Options])
%
%   [sse,Ops,Lhood] = ScoreModels(Models,Trajs,[Options])
%
%   Models must be formatted as a struct in the MCCV or CV format.
%

% Scott J Gaffney   10 October 2003
% Department of Information and Computer Science
% University of California, Irvine.
%
% Changes
% ---------------------------------


PROGNAME = 'ScoreModels';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%%
if (isvector(trajs))  % make sure that size(trajs,1) returns correct value
  trajs = trajs(:);
end
if (exist('mods')~=1 | isempty(mods))
  error([PROGNAME,': Models must be provided.']);
end
if (~isfield(mods,'State'))
  error([PROGNAME,': Models must contain the field .State.']);
end
Ops = cexist('Ops',[]);

DoLhood = 0;
if (nargout>2)
  DoLhood=1;
end

mccv=0;
if (isfield(mods,'beta'))
  mccv = 1;
end
%%%
%%% End Argument Processing



R = mods.NumRuns;
[K,Q] = size(mods.runs(1).Models);
sse = zeros(K,Q,R);
if (DoLhood), Lhood = zeros(K,Q,R); end
if (~isfield(Ops,'RandPointsPerRun') | isempty(Ops.RandPointsPerRun))
  Ops.RandPointsPerRun = cell(1,R);
end

if (mccv)
% Calculate sizes of training and testing sets
  TrajsLen = size(trajs,1);
  Ops = SetFieldDef(Ops,'SubsetSize',TrajsLen);
  TestLen = floor(Ops.SubsetSize * mods.beta);
  TrainLen  = Ops.SubsetSize - TestLen;

% Calculate Random Data Splits
  [trash,RandSplitState] = AdvRandState(mods.State,TrajsLen,R-1);
end

%Bring up Message Bar
Options.MsgHnd = msgbar([],'');

%% Loop over different MCCVRuns
for (r=1:R)
  if (mccv)
    rand('state',RandSplitState(:,r));
    RandIndex = randperm(TrajsLen);
    TestTrajs  = trajs(RandIndex(TrainLen+1:Ops.SubsetSize),:);
  else
    TestTrajs = trajs(mods.folds{r});
  end
  
  Ops.RandPoints = Ops.RandPointsPerRun{r};
  for k=1:K
    for q=1:Q
      msgbar(Options.MsgHnd,sprintf('Run %d, K %d, Order %d',r,k,q));
      [sse(k,q,r),Ops] = modelsse(mods.runs(r).Models(k,q),TestTrajs,Ops);
      if (DoLhood)
        Lhood(k,q,r) = model_like(mods.runs(r).Models(k,q),TestTrajs,Ops);
      end
    end
  end
  Ops.RandPointsPerRun{r} = Ops.RandPoints;
end

if (~isempty(Options.MsgHnd))
  delete(Options.MsgHnd);
end

Ops = rmfield(Ops,'RandPoints');


