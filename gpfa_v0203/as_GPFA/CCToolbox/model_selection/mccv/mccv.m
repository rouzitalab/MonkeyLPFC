function MCCV = mccv(trajs,MCCV)
%MCCV  Monte Carlo Cross Validation
%   MCCV = MCCV(TRAJS,MCCV)
%
%   MCCV
%    .State       : same as rand('state')
%    .beta        : test proportion
%    .SubsetSize  : subset of trajs to sample from
%    .NumRuns     : number of runs
%    .Krange      : vector of k-values to use, e.g., [1 3 5 6]
%    .OrderRange  : vector of orders to fit, e.g., [2 4 5]

% Scott J Gaffney   23 October 2001
% Department of Information and Computer Science
% University of California, Irvine.
%
% Changes
% ---------------------------------
% 13 March 2002 (SJG)
% Totally rewritten: placed MCCV runs in the file
%
% 26 March 2002 (SJG)
% Moved the actual data training and testing out to DoDataSplit.
% Changed the argument list to an input structure.
%

PROGNAME = 'mccv';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%%
if (isvector(trajs))  % make sure that size(trajs,1) returns correct value
  trajs = trajs(:);
end
if (exist('MCCV')~=1 | isempty(MCCV))
  error([PROGNAME,': MCCV must be provided.']);
end
MCCV = SetFieldDef(MCCV,'State',rand('state'));
MCCV = SetFieldDef(MCCV,'beta',.3);
MCCV = SetFieldDef(MCCV,'NumRuns',1);
%%%
%%% End Argument Processing



% Create Options Structure for DoDataSplit()
TrajsLen = size(trajs,1);
MCCV = SetFieldDef(MCCV,'SubsetSize',TrajsLen);
Options = MCCV;  % begin with the MCCV struct
Options = rmfield(Options,'NumRuns');
Options = rmfield(Options,'beta');
Options = SetFieldDef(Options,'MsgHnd',[]);

% Calculate sizes of training and testing sets
TestLen = floor(Options.SubsetSize * MCCV.beta);
TrainLen  = Options.SubsetSize - TestLen;

% Calculate Random Data Splits
tstate = rand('state');  % added 4 December 2002
[trash,RandSplitState] = AdvRandState(MCCV.State,TrajsLen,MCCV.NumRuns-1);
rand('state', tstate);  % added 4 December 2002

% Bring up Message Bar
CreatedMsgBar=0;
if (isempty(Options.MsgHnd))
  Options.MsgHnd = msgbar([],'');
  CreatedMsgBar=1;
end

%% Loop over different MCCVRuns
for (r=1:MCCV.NumRuns)
  TempState = rand('state');
  rand('state',RandSplitState(:,r));
  RandIndex = randperm(TrajsLen);
  rand('state',TempState);
  
  Options.MsgPrefix = sprintf('Run %d, ',r);
  TrainTrajs = trajs(RandIndex(1:TrainLen),:);
  TestTrajs  = trajs(RandIndex(TrainLen+1:Options.SubsetSize),:);
  MCCV.runs(r,1) = DoDataSplit(TrainTrajs,TestTrajs,Options);
  clear TrainTrajs TestTrajs;
end

if (CreatedMsgBar)
  delete(Options.MsgHnd);
end





