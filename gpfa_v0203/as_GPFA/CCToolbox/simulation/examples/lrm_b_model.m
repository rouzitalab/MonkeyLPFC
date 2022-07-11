function M = lrm_b_model()
%LRM_B_MODEL  Example simulation model
%  M = LRM_B_MODEL()

% Scott J Gaffney   8 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

M.Mu(:,1,1) = [5  8  ];
M.Mu(:,2,1) = [0  3  ]; 
M.Mu(:,3,1) = [10 .8  ];
M.order = size(M.Mu,1)-1;

M.Alpha = [0.4 0.3 0.3];
M.K = length(Alpha);

% S: K-D
M.S = [.8
       1.4
       1.6];

 
% Sigma: K-D
M.Sigma(1,1) = [20];
M.Sigma(2,1) = [20];
M.Sigma(3,1) = [20];



if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end