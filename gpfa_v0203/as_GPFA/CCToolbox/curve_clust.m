function MaxModel = curve_clust(varargin)
%Curve_Clust  Cluster trajectories using chosen method
%  See HELP CCToolbox for an explanation of the 'Trajs' argument, the Sequence 
%  format, the Feature Vector format, the Cell format, and all other
%  general concepts referred herein.
%
%  MaxModel = Curve_Clust(Trajs, ModelOptions)
%    - Trajs        : 'Trajs' structure (see HELP CCToolbox)
%    - ModelOptions : options structure used to pass required and/or optional
%                     arguments to specific clustering methods.
% 
%  In the following structure definition, common options are listed which
%  apply to many clustering methods. However, this is not an exhaustive
%  list nor are all required arguments for all methods listed; specific 
%  clustering methods WILL require fields not listed here. see HELP for 
%  the specific method for more information about what arguments are 
%  required and/or recommended. (R) denotes fields which are required
%  for Curve_Clust() itself.
%
%   MODEL_OPTIONS (structure)
%     .method      : (R) select clustering algorithm (see below)
%     .K           : (R) number of clusters
%     .order       : order of model to fit
%     .zero        : see Trajs2Seq()
%     .NumEMStarts : number of EM starts
%     .MsgHnd      : handle to MSGBAR figure (-1 to disable)
%     .MsgPrefix   : string to prepend to MSGBAR call
%
%   METHOD (string)
%     Call this function with the single argument 'methods' to display the
%     current list of acceptable methods. In other words, type 
%     curve_clust('methods') at the matlab prompt.

% Scott J Gaffney   5 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'curve_clust';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end
% Constants

%%% Handle Argument Processing
%%%
args = varargin; clear varargin;
n = length(args);
%
% Check for calling convention
%
% ()
if (n==1 & strcmp(args{1},'methods'))
  ListModels;
  return;
end

% (trajs,Ops)
if (n==2)
  trajs = args{1};
  Ops = args{2};
  
% (X,Y,Seq,Ops)
elseif (n==4)
  trajs.X = args{1};
  trajs.Y = args{2};
  trajs.Seq = args{3};
  Ops = args{4};
else
  error([PROGNAME,': incorrect parameter number/type/order specification.']);
end
%%%
%%% End Argument Processing
Ops = DefaultOptions(Ops);


%% Handle graphical message bar
CreatedMsgBar=0;
if (isempty(Ops.MsgHnd))
  Ops.MsgHnd = msgbar([],'');
  CreatedMsgBar=1;
end
% save the state
save_state = rand('state');  save_nstate = randn('state'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the requested cluster method
%%
% Get handles to the appropriate functions
info = listmodels(Ops.method);
capi_func = str2func(info.CurveAPI);
like_func = str2func(info.likelihood);

% run NumEMStarts random starts
MaxModel.TrainLhood = -inf;
for i=1:Ops.NumEMStarts
  UpdateMessage(Ops,i);
  M = feval(capi_func,trajs,Ops.K,Ops);
  if (isnan(M.TrainLhood)), continue; end  % catch boundary conditions
  if (isempty(M.TrainLhood)), break; end   % catch non-EM algorithms
  if (Ops.ValidateLhood & ~isempty(like_func))
    M.TrainLhood_ppt = feval(like_func,M,trajs,Ops);
    M.TrainLhood = M.TrainLhood_ppt * M.NumPoints;
  end
  if (M.TrainLhood > MaxModel.TrainLhood)
    MaxModel = M;
  end
end
if (MaxModel.TrainLhood==-inf), MaxModel=M; end
MaxModel.state = save_state; MaxModel.nstate = save_nstate;

% Cleanup
if (CreatedMsgBar)
  delete(Ops.MsgHnd);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SetDefaults
%%
function Options = DefaultOptions(Options)
Options = SetFieldDef(Options,'NumEMStarts',1);
Options = SetFieldDef(Options,'order',2);
Options = SetFieldDef(Options,'method','lrm');
Options = SetFieldDef(Options,'K','3');
Options = SetFieldDef(Options,'MsgHnd',[]);
Options = SetFieldDef(Options,'MsgPrefix','');
Options = SetFieldDef(Options,'zero','nozero');
Options = SetFieldDef(Options,'MinLen',[]);
Options = SetFieldDef(Options,'ShowGraphics',0);
Options = SetFieldDef(Options,'ValidateLhood',0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UpdateMessage
%%
function UpdateMessage(Model,i)
MsgPrefix = Model.MsgPrefix;
if (Model.MsgHnd>=0)
  Model.MsgPrefix = sprintf('%sStart %d ',MsgPrefix,i);
elseif (Model.MsgHnd~=-1)
  fprintf('%s: %sStart %d\n',Model.method,Model.MsgPrefix,i);
end







