function M = gmix_model_2()
%GMIX_MODEL_2  Example simulation model
%  M = GMIX_MODEL_2()

% Scott J Gaffney   8 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

% P-K-D
M.Mu(:,1,1) = [6.68    2      -0.2];
M.Mu(:,2,1) = [3.6303  2.1858  0.00]; 
M.Mu(:,3,1) = [2.7027  0.4160  0.00];
M.order = size(M.Mu,1)-1;  % not used for gmix
M.numdims = 1

M.Alpha = [0.4 0.33 0.27]';
M.K = length(M.Alpha);

% K-D
M.Sigma(1,1) = [1];
M.Sigma(2,1) = [5];
M.Sigma(3,1) = [10];




if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end