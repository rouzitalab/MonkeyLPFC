function  [mccvout,trajs,Krange,OrderRange] = gomccv(trajs,ops)
%GOMCCV   Simple function to execute MCCV()
%
%   MCCV_Out = GOMCCV(Trajs,[Ops])
%
%   [MCCV_Out, trajs, Krange, OrderRange] = GOMCCV(...)
%
%     state = rand('state',0);
%     AdvRandState(state,num,N) 
%
%   results in the generator ready to run iteration N+1 using MCCV. 
%   In particular, to begin on iteration 11 you would call 
%
%      state_11 = AdvRandState(state,num,10);
%      mccvout = gomccv(trajs,state_11);

% Scott J Gaffney   13 March 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'gomccv';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%% Begin Argument Processing
%
%
%% End Argument Processing

ops = cexist('ops',[]);
ops = SetFieldDef(ops,'method','');
ops = SetFieldDef(ops,'methods',{ops.method});
ops = SetFieldDef(ops,'Krange',[3]);
ops = SetFieldDef(ops,'OrderRange',[2]);
ops = SetFieldDef(ops,'NumRuns',1);
ops = SetFieldDef(ops,'NumEMStarts',5);
ops = SetFieldDef(ops,'SubsetSize',100);
ops = SetFieldDef(ops,'beta',0.3);
ops = SetFieldDef(ops,'zero','nozero');
ops = SetFieldDef(ops,'IterLimit',30);
ops = SetFieldDef(ops,'MinLen',[]);
ops = SetFieldDef(ops,'back','noback');
ops = SetFieldDef(ops,'State',rand('state'));
ops = SetFieldDef(ops,'MaxDec',6);
ops = SetFieldDef(ops,'Sigma',[]);
ops.Sigma = SetFieldDef(ops.Sigma,'Diagonal',1);
ops.Sigma = SetFieldDef(ops.Sigma,'Share',0);
ops = SetFieldDef(ops,'MsgHnd',[]);
ops = SetFieldDef(ops,'MsgPrefix','');
ops = SetFieldDef(ops,'DoSSE',0);


% make ops.method a cell array of strings if it is not already
if (isstr(ops.methods))
  ops.methods = {ops.methods};
end
ops.methods = ops.methods(:);


% Run MCCV over all of the specified methods
NumMethods = length(ops.methods);
for i=1:NumMethods
  ops.method = ops.methods{i};
  result = mccv(trajs,ops);
  setfield(mccvout,ops.method,result);
  ops.State = result.State;
end

