function [State,RandState] = AdvRandState(State,NumRand,NumCalls)
%AdvRandState  Generate state variables for successive calls to the generator.
%   State = AdvRandState(StartState,NumRand,NumCalls) returns the
%   state resulting from calling RANDPERM(NumRand) exactly NumCalls
%   different times after calling RAND('state',StartState).
%
%     state = rand('state',0);
%     AdvRandState(state,num,N) 
%
%   results in the generator ready to run iteration N+1 using MCCV. 
%   In particular, to begin on iteration 11 you would call 
%
%      state_11 = AdvRandState(state,num,10);
%      mccvout = gomccv(trajs,state_11);
%
%   [State,RandState] = AdvRandState(...) returns all of the NumCalls+1 states
%   in RandState.

% Scott J Gaffney   21 March 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'AdvRandState';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end



%%% Handle Argument Processing
%%%
if (exist('State')~=1)
  error([PROGNAME,': Argument Error: State must be provided.']);
end
if (isempty(State))
  State = rand('state');
end
if (exist('NumRand')~=1 | isempty(NumRand))
  error([PROGNAME,': Argument Error: NumRand must be provided.']);
end
if (exist('NumCalls')~=1 | isempty(NumCalls))
  error([PROGNAME,': Argument Error: NumCalls must be provided.']);
end
%%%
%%% End Argument Processing


% Calculate Random Data Splits
rand('state',State);
RandState(:,1)=State;
for r=2:NumCalls+1
  trash = randperm(NumRand);  % change the state
  clear trash;
  RandState(:,r) = rand('state');  % record the state
end
State = RandState(:,end);