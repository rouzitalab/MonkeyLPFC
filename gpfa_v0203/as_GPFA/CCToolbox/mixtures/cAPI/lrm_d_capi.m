function M = lrm_d_capi(trajs,K,Ops)
%LRM_D_CAPI  Curve API for LRM_D model.
%
%   Model = LRM_D_CAPI(Trajs,K,Options)
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find

% Scott Gaffney   9 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'lrm_d_capi';
METHOD = PROGNAME;
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

%%% Handle Argument Processing
%%%
n = nargin;
if (n<3)
  error([PROGNAME, ': incorrect number of parameters provided.']);
end
%%
%%% End Argument Processing

M = lrm_d(trajs,K,Ops.order,Ops);

