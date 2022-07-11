function M = lrm_cd_ab_model()
%LRM_CD_AB_MODEL  Example simulation model
%  M = LRM_CD_AB_MODEL()

% Scott J Gaffney   8 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

M.Mu(:,1) = [-39.3212   46.7218   -7.4532    0.2686];
M.Mu(:,2) = [26.0459  -20.0342    1.3855    0.0088]; 
M.Mu(:,3) = [-20  14  -2 .03];
M.order = size(M.Mu,1)-1;

M.Alpha = [0.4 0.33 0.27]';
M.K = length(Alpha);

% R: K
M.R = [.04
       .2
       .10];

% S: K
M.S = [1
       2
       4];

% U: K-D
M.U = [.03
       .2
       .1];
     
% T: K-D
M.T = [2
       .1
       1];

M.Sigma(1,1) = [20];
M.Sigma(2,1) = [20];
M.Sigma(3,1) = [20];




if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end