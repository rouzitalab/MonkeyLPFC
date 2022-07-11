function M = lrm_ab_model_2()
%LRM_AB_MODEL_2  Example simulation model
%  M = LRM_AB_MODEL_2()

% Scott J Gaffney   2 February 2003
% Department of Information and Computer Science
% University of California, Irvine.


% Alpha: K-1
M.Alpha = [.6 .4]';
M.K = length(Alpha);

% Mu: P-K-D
M.Mu(:,1,1) = [-8   6.0676   -1.2488    0.0751];
M.Mu(:,1,2) = [-4    8   -1.2488    .04];
M.Mu(:,2,1) = [0 3 0 0];
M.Mu(:,2,2) = [0 1 0 0];
M.order = size(M.Mu,1)-1;

% R: K-D
M.R(1,1) = .1;
M.R(2,1) = .1;

% S: K-D
M.S(1,1) = .1;
M.S(2,1) = .3;

% Sigma: K-D
M.Sigma(1,1) = .5;
M.Sigma(1,2) = .8;
M.Sigma(2,1) = .8;
M.Sigma(2,2) = .8;


if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end