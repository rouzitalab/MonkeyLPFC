function info = listmodels(method,GetType)
%ListModels  List the current cluster models supported in the CCToolbox.
%
%  ListModels()
%   If no output argument is given, then a list of models and descriptions
%   will be printed to the screen.
%
%  Info = ListModels() 
%   This will output a structure array with the following fields for each
%   supported cluster model in the CCToolbox.
%   .method      : CCToolbox method name for cluster model
%   .filename    : filename (function name) holding the cluster code
%   .description : cell array of strings providing a description of the model
%   .CurveAPI    : Curve API function name
%   .prediction  : posterior (or partial curve) prediction function name
%   .likelihood  : likelihood calculation function name
%   .simfun      : random simulation function 
%
%  Info = ListModels(method) 
%   This will output the same info array as above, but it will only
%   include the info structure that matches the method argument.
%
%  FunctionName = ListModels(method,GetType)
%   Returns the function associated with the clustering algorithm specified
%   in method. If GetType is given, then the actual function returned will
%   depend on the value of GetType (see below). Note that the value of method
%   can be any of the recognized "suffixed" versions as well as the standard
%   name provided in the .method field. Thus, ListModels() can be used to 
%   convert between suffixed and base name versions.
%
%   GetType
%     'like'     : return likelihood calculation function
%     'pred'     : return posterior (or partial curve) prediction function
%     'capi'     : return Curve API function
%     'view'     : return model view function
%     'sim'      : return model simulation function
%     'cluster'  : (otherwise) return cluster function


% Scott Gaffney   1 April 2005

% Default suffix for function type names
LIKE    = '_like';
PRED    = '_pred';
CAPI    = '_capi';

method = cexist('method',[]);
GetType = cexist('GetType',[]);

% Remove suffix
if (~isempty(method))
  j = findstr(method,LIKE);
  if (~isempty(j))
    method = method(1:j-1);
  else
    j = findstr(method,PRED);
    if (~isempty(j))
      method = method(1:j-1);
    else
      j = findstr(method,CAPI);
      if (~isempty(j))
        method = method(1:j-1);
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Start CCToolbox Cluster Models
%%
i=1;
info(i).method = 'lrm';
info(i).filename = 'lrm';
info(i).description = 'Linear Regression Mixtures (LRM) with y = XB + e';
info(i).view = 'lrm_view';
info(i).simfun = 'rm_sim';

i=i+1;
info(i).method = 'lrm_d';
info(i).filename = 'lrm_d';
info(i).description = 'LRM with transformation y=[x]B+d';
info(i).view = 'lrm_view';
info(i).simfun = 'rm_d_sim';

i=i+1;
info(i).method = 'lrm_cd';
info(i).filename = 'lrm_cd';
info(i).description = 'LRM with transformation y=c[x]B+d';
info(i).view = 'lrm_view';
info(i).simfun = 'rm_cd_sim';

i=i+1;
info(i).method = 'lrm_b';
info(i).filename = 'lrm_b';
info(i).description = 'LRM with transformation y=[x+b]B';
info(i).view = 'lrm_view';
info(i).simfun = 'rm_b_sim';

i=i+1;
info(i).method = 'lrm_ab';
info(i).filename = 'lrm_ab';
info(i).description = 'LRM with transformation y=[ax+b]B';
info(i).view = 'lrm_view';
info(i).simfun = 'rm_ab_sim';

i=i+1;
info(i).method = 'lrm_d_b';
info(i).filename = 'lrm_d_b';
info(i).description = 'LRM with transformation y=[x+b]B+d';
info(i).view = 'lrm_view';
info(i).simfun = 'rm_d_b_sim';

i=i+1;
info(i).method = 'lrm_d_ab';
info(i).filename = 'lrm_d_ab';
info(i).description = 'LRM with transformation y=[ax+b]B+d';
info(i).view = 'lrm_view';
info(i).simfun = 'rm_d_ab_sim';

i=i+1;
info(i).method = 'lrm_cd_b';
info(i).filename = 'lrm_cd_b';
info(i).description = 'LRM with transformation y=c[x+b]B+d';
info(i).simfun = 'rm_cd_b_sim';
info(i).view = 'lrm_view';

i=i+1;
info(i).method = 'lrm_cd_ab';
info(i).filename = 'lrm_cd_ab';
info(i).description = 'LRM with transformation y=c[ax+b]B+d';
info(i).view = 'lrm_view';
info(i).simfun = 'rm_cd_ab_sim';

i=i+1;
info(i).method = 'srm';
info(i).filename = 'srm';
info(i).description = 'Spline Regression Mixture (SRM) with y = XB + e';
info(i).view = 'srm_view';
info(i).simfun = 'rm_sim';

i=i+1;
info(i).method = 'srm_d';
info(i).filename = 'srm_d';
info(i).description = 'SRM with transformation y=[x]B+d';
info(i).view = 'srm_view';
info(i).simfun = 'rm_d_sim';

i=i+1;
info(i).method = 'srm_cd';
info(i).filename = 'srm_cd';
info(i).description = 'SRM with transformation y=c[x]B+d';
info(i).view = 'srm_view';
info(i).simfun = 'rm_cd_sim';

i=i+1;
info(i).method = 'srm_b';
info(i).filename = 'srm_b';
info(i).description = 'SRM with transformation y=[x+b]B';
info(i).view = 'srm_view';
info(i).simfun = 'rm_b_sim';

i=i+1;
info(i).method = 'srm_ab';
info(i).filename = 'srm_ab';
info(i).description = 'SRM with transformation y=[ax+b]B';
info(i).view = 'srm_view';
info(i).simfun = 'rm_ab_sim';

i=i+1;
info(i).method = 'srm_d_b';
info(i).filename = 'srm_d_b';
info(i).description = 'SRM with transformation y=[x+b]B+d';
info(i).view = 'srm_view';
info(i).simfun = 'rm_d_b_sim';

i=i+1;
info(i).method = 'gmix';
info(i).filename = 'gmix';
info(i).description = 'Gaussian Mixtures';
info(i).view = 'gmix_view';
info(i).simfun = 'gmix_sim';

% install the default function type names
for j=1:length(info)
  info(j).CurveAPI   = [info(j).filename, CAPI];
  info(j).prediction = [info(j).filename, PRED];
  info(j).likelihood = [info(j).filename, LIKE];
end

i=i+1;
info(i).method = 'kmeans';
info(i).filename = 'kmeans';
info(i).description = 'K-means';
info(i).CurveAPI   = [];
info(i).prediction = 'kmeans_pred';
info(i).likelihood = [];
info(i).view = 'gmix_view';
info(i).simfun = [];
%%
%%%% End CCToolbox Cluster Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% List the models if there are no input args
if (nargin==0)
  if (nargout==0)
    n = length(info);
    for i=1:n
      fprintf('%-11s- %s\n',info(i).method,info(i).description);
    end
    clear info;
  end
  return;
end

for i=1:length(info)
  if (strcmp(method,info(i).method))
    if (isempty(GetType))
      info = info(i);
    else
      switch GetType
        case 'like'
          info = info(i).likelihood;
        case 'pred'
          info = info(i).prediction;
        case 'capi'
          info = info(i).CurveAPI;
        case 'view'
          info = info(i).view;
        case 'simfun'
          info = info(i).simfun;
        otherwise
          info = info(i).filename;
      end
    end
    return;
  end
end
errorbox(sprintf('Unrecognized clustering method: %s',method),...
  'modelsse: Argument Error');
error(sprintf('Unrecognized clustering method: %s',method));
