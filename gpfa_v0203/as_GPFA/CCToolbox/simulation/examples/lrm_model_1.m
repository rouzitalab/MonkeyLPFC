function M = lrm_model_1()
%LRM_MODEL_1  Example simulation model
%  M = LRM_MODEL_1()

% Scott J Gaffney   8 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

M.Mu(:,1,1) = [0  8  -0.2];
M.Mu(:,2,1) = [0  3   0.00]; 
M.Mu(:,3,1) = [0 .8   0.00];
M.order = size(M.Mu,1)-1;

M.Alpha = [0.4 0.3 0.3];
M.K = length(Alpha);

M.Sigma(1,1) = [20];
M.Sigma(2,1) = M.Sigma(1);
M.Sigma(3,1) = M.Sigma(1);




if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end