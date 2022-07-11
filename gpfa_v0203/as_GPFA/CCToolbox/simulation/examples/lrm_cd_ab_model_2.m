function M = lrm_cd_ab_model_2()
%LRM_CD_AB_MODEL_2  Example simulation model
%  M = LRM_CD_AB_MODEL_2()

% Scott J Gaffney   2 February 2003
% Department of Information and Computer Science
% University of California, Irvine.


% Alpha: K-1
M.Alpha = [1]';
M.K = length(Alpha);

% Mu: P-K-D
M.Mu(:,1) = [20  0  -8 .8];
M.order = size(M.Mu,1)-1;

% R: K-D
M.R = [.3];

% S: K
M.S = [2];

% U: K-D
M.U = [.3];
     
% T: K-D
M.T = [2];
     
% Sigma: K-D
M.Sigma(1,1) = 1;


if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end