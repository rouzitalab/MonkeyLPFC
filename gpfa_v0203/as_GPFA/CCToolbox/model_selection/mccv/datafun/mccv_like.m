function  [Mu,Sigma,trainMu,trainSigma] = mccv_like(mccv)
%MCCV_LIKE   Find mean and standard deviation of likelihood from MCCV runs.
%   [Mu,Sigma,trainMu,trainSigma] = MCCV_LIKE(MCCV_Struct)
%     or
%   [Mu,Sigma,trainMu,trainSigma] = MCCV_LIKE(Array_Struct)

% Scott J Gaffney   13 March 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'mccv_like';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%% Begin Argument Processing
%
Mu=[]; Sigma=[]; trainMu=[]; trainSigma=[];
if (exist('mccv')~=1 | ~isstruct(mccv))
  errorbox('Argument Error: the first argument must be a STRUCT.',PROGNAME);
  return;
end  
if (~isfield(mccv,'TestLike'))
  return;
end
if (isfield(mccv,'runs'))
  mccv = mccv.runs;
end
%
%% End Argument Processing


mccv = mccv(:);
len = length(mccv);
for i=1:len
  like(:,:,i) = mccv(i).TestLike;  % this is already a per point measure
  trainlike(:,:,i) = reshape([mccv(i).Models.TrainLhood],size(mccv(i).TestLike));
  trainlike(:,:,i) = trainlike(:,:,i) / mccv(i).Models(1).NumPoints;
end
clear mccv;

Mu = mean(like,3);
trainMu = mean(trainlike,3);
for j=1:size(like,2)
  Sigma(:,j) = std(like(:,j,:),0,3);
  trainSigma(:,j) = std(squeeze(trainlike(:,j,:)),0,2);
end

