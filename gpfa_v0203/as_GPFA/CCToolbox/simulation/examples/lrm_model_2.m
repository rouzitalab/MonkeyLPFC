function M = lrm_model_2()
%LRM_MODEL_2  Example simulation model
%  M = LRM_MODEL_2()

% Scott J Gaffney   8 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

M.Mu(:,1,1) = [25    2      -0.2];
M.Mu(:,2,1) = [3.6303  2.1858  0.00]; 
M.Mu(:,3,1) = [2.7027  0.4160  0.00];
M.order = size(M.Mu,1)-1;

M.Alpha = [0.4 0.33 0.27];
M.K = length(Alpha);

M.Sigma(1,1) = [100];
M.Sigma(2,1) = M.Sigma(1);
M.Sigma(3,1) = M.Sigma(1);




if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end