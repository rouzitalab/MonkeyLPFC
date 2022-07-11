function M = gmix_model_1()
%GMIX_MODEL_1  Example simulation model
%  M = GMIX_MODEL_1()

% Scott J Gaffney   2 February 2003
% Department of Information and Computer Science
% University of California, Irvine.

% m-by-K  (n_i-by-K)
M.Mu(:,1) = [1 4 3 5 2 1 2 3 7 3]';
M.Mu(:,2) = [4 7 9 0 2 3 4 8 4 5]'; 
M.Mu(:,3) = [8 2 6 1 4 9 5 3 7 2]';
M.order = size(M.Mu,1)-1;  % not used for gmix
M.numdims = 1;

M.Alpha = [0.4 0.33 0.27]';
M.K = length(M.Alpha);

% m-by-m-by-K
M.Sigma(:,:,1) = [4];
M.Sigma(:,:,2) = M.Sigma(:,:,1);
M.Sigma(:,:,3) = M.Sigma(:,:,1);


if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end