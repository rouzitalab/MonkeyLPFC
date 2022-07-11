function M = lrm_b_capi(trajs,K,Ops)
%LRM_B_CAPI  Curve API for LRM_B model.
%
%   Model = LRM_B_CAPI(Trajs,K,Options)
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find

% Scott Gaffney   9 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'lrm_b_capi';
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

M = lrm_b(trajs,K,Ops.order,Ops);

