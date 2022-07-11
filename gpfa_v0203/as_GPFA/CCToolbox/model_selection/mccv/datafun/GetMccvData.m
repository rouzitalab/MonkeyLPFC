function index = GetMccvData(Mccv,Run,NumTrajs,SetType)
%GetMccvData  Get the training/testing random data split.
%   Index = GetMccvData(Mccv,Run,NumTrajs,SetType) 
%
%   SetType:
%     'train'    : return training set indicies
%     'test'     : return testing set indicies
%     otherwise  : return full random index

% Scott J Gaffney   05 April 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'GetMccvData';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end



%%% Handle Argument Processing
%%%
% Mccv = SetFieldDef(Mccv,'beta',Mccv.Beta);  % deal w/ capitalization error
%%%
%%% End Argument Processing

% Calculate sizes of training and testing sets
TestLen = floor(Mccv.SubsetSize * Mccv.beta);
TrainLen  = Mccv.SubsetSize - TestLen;

% Calculate Random Data Splits
tstate = rand('state');  % added 4 December 2002
NewState = AdvRandState(Mccv.State,NumTrajs,Run-1);

rand('state',NewState);
index = randperm(NumTrajs);
switch (SetType)
case 'test'
  index = index(TrainLen+1:Mccv.SubsetSize);
otherwise
  index = index(1:TrainLen);
end
rand('state', tstate);  % added 4 December 2002
