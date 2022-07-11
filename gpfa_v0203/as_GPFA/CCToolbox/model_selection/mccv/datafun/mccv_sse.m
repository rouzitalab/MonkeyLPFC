function  [Mu,Sigma] = mccv_sse(mccv)
%MCCV_SSE   Find mean and standard deviation of SSE from MCCV runs.
%   [Mu,Sigma] = MCCV_SSE(MCCV_Struct)
%     or
%   [Mu,Sigma,trainMu,] = MCCV_SSE(Array_Struct)

% Scott J Gaffney   13 March 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'mccv_sse';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%% Begin Argument Processing
%
Mu=[]; Sigma=[];
if (exist('mccv')~=1 | ~isstruct(mccv))
  errorbox('Argument Error: the first argument must be a STRUCT.',PROGNAME);
  return;
end  
if (~isfield(mccv,'TestSSE'))
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
  like(:,:,i) = mccv(i).TestSSE;  % this is already a per point measure
end

Mu = mean(like,3);
for j=1:size(like,2)
  Sigma(:,j) = std(like(:,j,:),0,3);
end

