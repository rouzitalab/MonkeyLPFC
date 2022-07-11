function cvstruct = CV(trajs,cvstruct)
%cvstruct  Run cross validation on a set of curves.
%   cvstruct = CV(TRAJS,cvstruct)
%
%   cvstruct
%     .runs        : contains results from each DoDataSplit() call
%     .folds       : (optional) cell array of indexes for each fold
%     .NumRuns     : number of folds to perform
%     .Krange      : set of K-values to run on each fold
%     .OrderRange  : set of order values to run on each fold
%     .NumEMStarts : number of EM starts for each training session
%     .method      : EM method to run (see Curve_Clust())
%     .zero        : see Curve_Clust()
%     .back        : see Curve_Clust()
%     .State       : call rand('state',cvstruct.State) before we run
%

% Scott J Gaffney   2 February 2003
% Department of Information and Computer Science
% University of California, Irvine.


PROGNAME = 'crossval';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%%
if (isvector(trajs))  % make sure that size(trajs,1) returns correct value
  trajs = trajs(:);
end
if (exist('cvstruct')~=1 | isempty(cvstruct))
  error([PROGNAME,': CV must be provided.']);
end
cvstruct = SetFieldDef(cvstruct,'State',rand('state'));
cvstruct = SetFieldDef(cvstruct,'folds',[]);
%%%
%%% End Argument Processing


% Create Options Structure for DoDataSplit()
Options = cvstruct;  % begin with the CV struct

% Calculate folds
TrajsLen = size(trajs,1);
if (isempty(cvstruct.folds))
  FoldLen = round(TrajsLen / cvstruct.NumRuns);
  rand('state',cvstruct.State);
  RandIndex = randperm(TrajsLen);
  from_i=1;
  for (i=1:cvstruct.NumRuns-1)
    to_i = from_i+FoldLen-1;
    cvstruct.folds{i} = RandIndex(from_i:to_i);
    from_i = to_i+1;
  end
  cvstruct.folds{i+1} = RandIndex(from_i:end);  % get the last one
  clear FoldLen RandIndex from_i to_i i;
end

% Bring up Message Bar
Options.MsgHnd = msgbar([],'');

%% Loop over different CV folds
for (r=1:cvstruct.NumRuns)
  Options.MsgPrefix = sprintf('Run %d, ',r);
  train_index = (1:TrajsLen);
  train_index(cvstruct.folds{r}) = [];  % we train on everything except for which we test
  cvstruct.runs(r,1) = DoDataSplit(trajs(train_index,:),trajs(cvstruct.folds{r},:),Options);
end

if (~isempty(Options.MsgHnd))
  delete(Options.MsgHnd);
end



