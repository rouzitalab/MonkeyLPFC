function M = srm_d_b_capi(trajs,K,Ops)
%SRM_D_B_CAPI  Curve API for SRM_D_B model.
%
%   Model = SRM_D_B_CAPI(Trajs,K,Options)
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find

% Scott Gaffney   9 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'srm_d_b_capi';
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


% Ops Handling
Ops = SetFieldDef(Ops,'order','4');
Ops = SetFieldDef(Ops,'knots',[]);

M = srm_d_b(trajs,K,Ops.knots,Ops.order,Ops);

