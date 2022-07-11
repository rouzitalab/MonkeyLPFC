function M = gmix(X,K,Options)
%GMIX  Run EM for Gaussian mixtures
%   Model = GMIX(X,K,Options) returns the ML estimates 
%   for the mixture coefficients W, the means MU, and the covariance matrices
%   SIGMA. X contains the data points, one per row. K is the number
%   of classes to be found in the data.
%
%   OPTIONS
%     .Sigma
%       .Diagonal       : 1 - use diagonal covariance, 0 - otherwise
%       .Share          : 1 - estimate a single covariance, 0 - otherwise
%
%   w is a K by 1 vector, mu is a D by K matrix, and sigma is
%   a D by D by K matrix, where D is the dimension of the data
%   provided in X. GMIX also returns the log-likelihood as a
%   function of iterations in Lhood.

% Scott Gaffney   10 June 1998
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------
% 31 May 2002 (SJG)
%   Replaced all REPMAT calls and added an initialization
%   routine. Other than that, I didn't really care to
%   change this old code.

PROGNAME = 'gmix';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

Options = cexist('Options',[]);
Options = SetFieldDef(Options,'stopval',1e-5);
Options = SetFieldDef(Options,'minvar',1e-5);
Options = SetFieldDef(Options,'IterLimit',50);
Options = SetFieldDef(Options,'zero','nozero');
Options = SetFieldDef(Options,'Sigma',[]);
Options.Sigma = SetFieldDef(Options.Sigma,'Share',0);
Options.Sigma = SetFieldDef(Options.Sigma,'Diagonal',1);
Options.Sigma = SetFieldDef(Options.Sigma,'Single',0);
Options = SetFieldDef(Options,'NumEMStarts',1);
Options = SetFieldDef(Options,'ShowGraphics',0);

MaxTrys = 5;
for i=1:MaxTrys
  M = em(X,K,Options,PROGNAME);
  if (~isempty(M)), break;  end
  fprintf([PROGNAME, ': starting over due to singularity issues\n\n']);
end

% Calculate number of independent parameters
[D,K] = size(M.Mu);
if (M.Options.Sigma.Share==1), shK = 1; else  shK=K; end
if (M.Options.Sigma.Diagonal==1)
  M.NumIndParams = (K-1) + K*D + shK*D;   % alpha, mu, sigma
else
  M.NumIndParams = (K-1) + K*D + shK*D*(D+1)/2;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM
%
function Model = em(X,K,Options,PROGNAME)

% Preprocessing
[N,D] = size(X);
diagi = (1 + (D*(0:D-1)) + (0:D-1))' * ones(1,K);
diagi = diagi + ones(D,1) * (D^2*(0:K-1));
Pik = zeros(N,K);
num = 0;
Model.method = 'gmix';
Model.version = '1.1';
Model.zero = Options.zero;


% Initialize the algorithm
[w,mu,sigma] = InitEM(X,K,diagi,Options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The EM Procedure
%
%fprintf(['GMIX',': beginning parameter estimation\n']);
lastwarn('');
while (1)
  num = num + 1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Perform the E-Step
  %

  % Enforce a minimum variance value
  badi = diagi(find(sigma(diagi)<Options.minvar));
  if (~isempty(badi))
      fprintf([PROGNAME,': an illegal variance was detected\n']);
      sigma(badi) = Options.minvar;
  end
  
  % Calculate the new posteriors (will be normalized later)
  for c = 1:K
    Pik(:,c) = mvnormpdf(X,mu(:,c),sigma(:,:,c)) * w(c);
  end

  % Calculate the log-likelihood
  s = sum(Pik,2);
  zero = find(s==0);
  if (~isempty(zero))
    fprintf([PROGNAME, ': log(0) detected at iteration %d\n'],num);
    Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*w');
    s(zero) = sum(Pik(zero,:),2);
  end
  loglhood(num) = sum(log(s));
  % fprintf('Iteration %i, Log-likelihood: %f\n',num,loglhood(num));
  
  % Normalize the posteriors before we go anywhere.
  Pik = Pik ./ (sum(Pik,2)*ones(1,K));
  
  
  
  % Should we stop the iterations?
  if (~isempty(lastwarn))
    if (strncmp(lastwarn,'Matrix is sing',14))
      fprintf([PROGNAME,': Sigma is singular, giving up on this attempt.\n']);
      Model=[];
      return;
    else
      lastwarn('');  % reset warning
    end
  end
  if (num>3)
    if (isnan(loglhood(num)))
      fprintf([PROGNAME,': the log-likelihood is equal to NaN; giving up\n\n']);
      loglhood(isnan(loglhood)) = -inf;
      break;
    end
    
    Lag    = loglhood(num)-loglhood(num-1);
    Change = loglhood(num)-loglhood(2);
    if (Change==0)
      %fprintf([PROGNAME,': There was no change from initial position; giving up.\n']);
      break;
    elseif (Lag/Change < Options.stopval | num>=Options.IterLimit)
      break;
    end
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Perform the M-Step
  %
  w = sum(Pik)'./ N;
  mu = (X'*Pik) ./ (ones(D,1)*sum(Pik));
  
  for c = 1:K
    XMu = X - ones(N,1)*mu(:,c)';
    PikXMu = Pik(:,c)*ones(1,D) .*XMu;
    sigma(:,:,c) = PikXMu'*XMu;
    if (~Options.Sigma.Share)
      sigma(:,:,c) = sigma(:,:,c)./sum(Pik(:,c));
    end
    if (Options.Sigma.Diagonal)
      sigma(:,:,c) = diag(diag(sigma(:,:,c)));
    end
    if (Options.Sigma.Single)
      sigma(:,:,c) = diag(ones(D,1)*sum(diag(sigma(:,:,c)))./D);
    end
  end 
  if (Options.Sigma.Share)
    sigma = sum(sigma,3)./N;
    sigma = sigma(:,:,ones(1,K));
  end
    
  
  if (Options.ShowGraphics & mod(num,2))
    [trash,C] = max(Pik,[],2);
    M.Mu = mu; M.Options = Options; M.C = C; M.method = 'gmix';
    M.Options.MinLen = [];
    Options = viewmodel(M,X,Options);
  end
  
  
end % {while}



% Prepare Outputs
Model.K = K;
Model.Mu = mu;
Model.Sigma = sigma;
Model.Alpha = w;
Model.Lhood = loglhood;
Model.NumPoints = prod(size(X));
Model.TrainLhood = Model.Lhood(end);
Model.TrainLhood_ppt = Model.Lhood(end)./Model.NumPoints;
Model.Pik = Pik;
[trash,Model.C] = max(Pik,[],2);
Model.Options = Options;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InitEM
%
function [kw,kmu,ksigma] = InitEM(X,K,diagi,Options)
[N,D] = size(X);

KMEANS = 1;
if (KMEANS == 1)
    [kmu, klabels] = kmeans(X,K);
else
    [kmu, klabels] = randmeans(X,K);
end

for c = 1:K
    kw(c,1) = length(find(klabels==c)) ./ N;
end


ksigma = zeros(D,D,K);
for c = 1:K
  [trash, ksigma(:,:,c)] = normmle(X(find(klabels==c),:));
end
badi = diagi(find(ksigma(diagi)<Options.minvar));
if (~isempty(badi))
  ksigma(badi) = Options.minvar;
end


if (Options.Sigma.Share)
  ksigma = repmat(sum(ksigma,3)/N,[1 1 K]);
end
if (Options.Sigma.Diagonal | Options.Sigma.Single)
  for (k=1:K)
    if (Options.Sigma.Single)
      ksigma(:,:,k) = diag(ones(D,1)*sum(diag(ksigma(:,:,k)))./D);
    else
      ksigma(:,:,k) = diag(diag(ksigma(:,:,k)));
    end
  end
end


