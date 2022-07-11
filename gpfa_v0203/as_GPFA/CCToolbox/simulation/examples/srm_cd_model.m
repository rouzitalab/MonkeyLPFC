function M = srm_cd_model()
%SRM_CD_MODEL  Example simulation model
%  M = SRM_CD_MODEL()

% Scott J Gaffney   2 February 2003
% Department of Information and Computer Science
% University of California, Irvine.

%% Setup spline basis model
% range of x-axis
b_int = 0;    % beginning of data interval
e_int = 20;    % end of data interval
len_int = e_int-b_int+1; % length of interval

degree = 3;  % degree of polynomial pieces, this can be calculated as Kn-L-1
             % where L is size(M.Mu,1)
Kn = max(ceil(len_int/3),1);  % number of knots
knots = linspace(b_int,e_int,Kn);
M.knots = addendpts(knots,degree+1);
Kn = length(M.knots);  % reset M to (possibly) new value
M.order = degree+1;

% Mu: L-K-D or P-K-D
M.Mu(:,1) = [10.0000    8.0000    9.0000    6.6667    6.0000 ...
    4.0000    3.0000    5.0000    2.0000];
M.Mu(:,2) = [4.0000    1.1111    3.3333    6.6667    6.0000 ...
    13.3333   16.6667   18.8889   10.0000];


% Alpha: K-1
M.Alpha = [.6 .4]';
M.K = length(Alpha);

% T: K-D
M.T = [1.4; .6];

% U: K-D
M.U = [.1; .5];

% Sigma: K-D
M.Sigma(1,1) = .05;
M.Sigma(2,1) = .05;


if (nargout==0)
  assignin('caller','SimModel',M);
  clear M;    % keeps default output from printing on-screen
end