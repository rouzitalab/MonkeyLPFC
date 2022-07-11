function [X,Y] = alignmat(M,x,Seq,Y)
%ALIGNMAT  Generate aligned curves according to model.
%
%   [X,Y] = AlignMat(M,x,Seq,Y) align data according to model M.

% Scott J Gaffney   10 September 2003
% School of Computer Science
% University of California, Irvine.

PROGNAME = 'AlignMat';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

if (nargin<4), Y=[]; end
n = length(Seq)-1;
lens = diff(Seq);
N = length(x);


% Note: the model is  y = e(ax-b) + f

switch (M.method)
case {'lrm','srm'}
  X = x;

case {'lrm_cd','lrm_d','srm_cd','srm_d'}
  X = x;
  D = size(M.Mu,3);
  if (isfield(M,'Ef'))  % do we have a trained model?
    [n,K,shareD] = size(M.Ef);
    for d=1:D
      shd = min(d,shareD);
      indx = sub2ind([n,K,shareD],(1:n)',M.C,shd(ones(n,1)));
      Y(:,d) = Y(:,d)-copy(M.Ef(indx),lens);
    end
  else   % we must have a true model
    [n,shareD] = size(M.f);
    for d=1:D
      shd = min(d,shareD);
      Y(:,d) = Y(:,d)-copy(M.f(:,shd),lens);
    end
  end  
  
case {'lrm_ab','srm_ab'}
  if (isfield(M,'Ea'))  % do we have a trained model?
    [n,K,D] = size(M.Ea);
    X = zeros(N,D);
    for d=1:D
      indx = sub2ind([n,K,D],(1:n)',M.C,d(ones(n,1)));
      X(:,d) = copy(M.Ea(indx),lens).*x -copy(M.Eb(indx),lens);
    end
  else  % we must have a true model
    [n,D] = size(M.a);
    X = zeros(N,D);
    for d=1:D
      X(:,d) = copy(M.a(:,d),lens).*x -copy(M.b(:,d),lens);
    end
  end
  
case {'lrm_b','srm_b'}
  if (isfield(M,'Eb'))  % do we have a trained model?
    [n,K,D] = size(M.Eb);
    X = zeros(N,D);
    for d=1:D
      indx = sub2ind([n,K,D],(1:n)',M.C,d(ones(n,1)));
      X(:,d) = x - copy(M.Eb(indx),lens);
    end
  else   % we must have a true model
    [n,D] = size(M.b);
    for d=1:D
      X(:,d) = x - copy(M.b(:,d),lens);
    end
  end
  
case {'lrm_d_b','srm_d_b','lrm_cd_b','srm_cd_b'}
  if (isfield(M,'Eb'))  % do we have a trained model?
    [P,K,D] = size(M.Mu);
    [n,K,shareB] = size(M.Eb);
    [n,K,shareF] = size(M.Ef);
    X = zeros(N,shareB);
    for d=1:D
      shb = min(d,shareB);
      shf = min(d,shareF);
      indxb = sub2ind([n,K,shareB],(1:n)',M.C,shb(ones(n,1)));
      indxf = sub2ind([n,K,shareF],(1:n)',M.C,shf(ones(n,1)));
      X(:,shb) = x - copy(M.Eb(indxb),lens);  % this gets recalc'ed if sharing
      Y(:,d) = Y(:,d)-copy(M.Ef(indxf),lens);
    end
  else   % we must have a true model
    [P,K,D] = size(M.Mu);
    [n,shareB] = size(M.b);
    [n,shareF] = size(M.f);
    X = zeros(N,shareB);
    for d=1:D
      shb = min(d,shareB);
      shf = min(d,shareF);
      X(:,shb) = x - copy(M.b(:,shb),lens);  % this gets recalc'ed if sharing
      Y(:,d) = Y(:,d)-copy(M.f(:,shf),lens);
    end
  end  

case {'lrm_d_ab','srm_d_ab','lrm_cd_ab','srm_cd_ab'}
  if (isfield(M,'Eb'))  % do we have a trained model?
    [P,K,D] = size(M.Mu);
    [n,K,shareB] = size(M.Eb);
    [n,K,shareF] = size(M.Ef);
    X = zeros(N,shareB);
    for d=1:D
      shb = min(d,shareB);
      shf = min(d,shareF);
      indxb = sub2ind([n,K,shareB],(1:n)',M.C,shb(ones(n,1)));
      indxf = sub2ind([n,K,shareF],(1:n)',M.C,shf(ones(n,1)));
      X(:,shb) = copy(M.Ea(indxb),lens).*x - ...
                 copy(M.Eb(indxb),lens);  % this gets recalc'ed if sharing
      Y(:,d) = Y(:,d)-copy(M.Ef(indxf),lens);
    end
  else   % we must have a true model
    [P,K,D] = size(M.Mu);
    [n,shareB] = size(M.b);
    [n,shareF] = size(M.f);
    X = zeros(N,shareB);
    for d=1:D
      shb = min(d,shareB);
      shf = min(d,shareF);
      X(:,shb) = copy(M.a(:,shb),lens).*x - ...
                 copy(M.b(:,shb),lens);  % this gets recalc'ed if sharing
      Y(:,d) = Y(:,d)-copy(M.f(:,shf),lens);
    end
  end  

end