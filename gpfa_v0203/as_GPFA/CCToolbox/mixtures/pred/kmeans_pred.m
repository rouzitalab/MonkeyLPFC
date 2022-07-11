function [Yhat,Pmod] = kmeans_pred(M,pt,x,Y,Seq,maxx)
%KMEANS_PRED  Make a partial curve prediction with Kmeans.
%
% IMPORTANT:
%  This function does not work for n > 1, where n=length(Seq)-1.
%  For more information about this function see LRM_PRED.
%
%   [Yhat,PostModel] = KMEANS_PRED(M,pt,x,Y,Seq,['max']) when pt~=[]
%
%   Note that Kmeans does not use the 'max' argument. Kmeans
%   always picks the 'max-class'.

% Scott Gaffney   10 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'kmeans_pred';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

maxx = 1;  % kmeans does not have a 'sum' option


[P,K,D] = size(M.Mu);
M.Mu = permute(M.Mu,[1 3 2]);
n = length(Seq)-1;
N = Seq(end)-1;
if (n>1)
  error('This function does not work for more than one curve at a time.\n');
end
mlen = max(diff(Seq));


if (isempty(x))
  Pmod.Pik = M.Alpha';
else
  indx = x+1;
  y = Y(:,:,ones(1,K));
  mu = M.Mu(indx,:,:);
  Pk = permute(sum(sum((y-mu).^2,1),2),[1 3 2]);
  Pmod.Pik = Pk ./ sum(Pk); 
  Pmod.Pik = 1./Pmod.Pik;              % convert to prob. scale using...
  Pmod.Pik = Pmod.Pik./sum(Pmod.Pik);  % ...inverse-proportion to error
  [trash, Pmod.C] = max(Pmod.Pik,[],2);
end


% Simply return if no prediction is requested
Yhat = [];
if (isempty(pt))
  return;
end

% Generate prediction at pt
Yhat = M.Mu(pt+1,:,Pmod.C);







