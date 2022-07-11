function Out = DoDataSplit(TrainTrajs,TestTrajs,Options)
%DoDataSplit  Calculate test statistics over K and Order for a data split.
%   Out = DoDataSplit(TrainTrajs,TestTrajs,Options)
%
%   Required options are listed below, but all of the usual model options
%   are allowed here (e.g., see Curve_Clust()).
%   Options:
%     .Krange       : vector of k-values to use, e.g., [1 3 5 6]
%     .OrderRange   : vector of orders to fit, e.g., [2 4 5]
%     .method       : see ListModels()

% Scott J Gaffney   26 March 2002
% Department of Information and Computer Science
% University of California, Irvine.
%

PROGNAME = 'DoDataSplit';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%%
if (exist('Options')~=1 | isempty(Options))
  error([PROGNAME,': Options must be provided']);
end
Options = SetFieldDef(Options,'NumEMStarts',10);
Options = SetFieldDef(Options,'MsgPrefix','');
Options = SetFieldDef(Options,'MsgHnd',[]);
Options = SetFieldDef(Options,'DoSSE',0);
InMsgPrefix = Options.MsgPrefix;
%%%
%%% End Argument Processing



% Handle graphical message bar
CreatedMsgBar=0;
if (isempty(Options.MsgHnd))
  Options.MsgHnd = msgbar([],'');
  CreatedMsgBar=1;
end

%% Begin the analysis for this data split
more off;
NumOrder=0;
% Loop over all orders
for CurOrder=Options.OrderRange(:)'
  NumOrder = NumOrder+1;
  Options.order = CurOrder;  % set order for fitting a particular model

  %% Loop over all K
  NumK=0;
  for k=Options.Krange(:)'
    NumK = NumK+1;
    
    % Fit model (EM)
    Options.K = k;  % set K for fitting a particular model
    Options.MsgPrefix = sprintf('%sOrder %d, K %d, ',InMsgPrefix,CurOrder,k);
    ModelOut = curve_clust(TrainTrajs,Options);
    
    % Test SSE
    fprintf('Testing....\n');
    if (Options.DoSSE)
      Out.TestSSE(NumK,NumOrder) = modelsse(ModelOut,TestTrajs,Options);
    end

    % Test Like and Train BIC
    Nk = ModelOut.NumIndParams;
    N = ModelOut.NumPoints;
    likefunc = listmodels(Options.method,'like');
    if (~isempty(likefunc))
      Out.TestLike(NumK,NumOrder) = feval(likefunc,ModelOut,TestTrajs,Options);
      Out.bic(NumK,NumOrder) = (ModelOut.TrainLhood- Nk/2*log(N))/N;
    end
    Out.Models(NumK,NumOrder) = ModelOut;  clear ModelOut;
  end
end 

if (CreatedMsgBar)
  delete(Options.MsgHnd);
end
