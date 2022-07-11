function [Yhat,Pmod] = gmix_pred(M,pt,x,Y,Seq,maxx)
%GMIX_PRED  Make a partial curve prediction with a Gaussian Mixture model.
%
% IMPORTANT:
%  This function does not work for n > 1, where n=length(Seq)-1.
%  For more information about this function see LRM_PRED.
%
%   [Yhat,PostModel] = GMIX_PRED(M,pt,x,Y,Seq,['max']) when pt~=[]
%
%   If you pass the string 'max' as the last argument, then Yhat is
%   calculated from the class w/ maximum membership probability instead
%   of summing across Pik as in the default case.

% Scott Gaffney   10 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'gmix_pred';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

maxx = cexist('maxx',0);
if (isstr(maxx) & strcmp(maxx,'max'))
  maxx = 1;
else
  maxx = 0;
end


[P,K,D] = size(M.Mu);
n = length(Seq)-1;
N = Seq(end)-1;
if (n>1)
  error('This function does not work for more than one curve at a time.\n');
end
mlen = max(diff(Seq));


if (isempty(x))
  Pmod.Pik = M.Alpha';
else
  Pk = ones(1,K);
  indx = x+1;
  for k = 1:K
    for d=1:D
      Pk(k) = Pk(k)*mvnormpdf(Y(:,d)',M.Mu(indx,k,d),M.Sigma(indx,indx,k,d));
    end  
  end
  Pk = Pk .*M.Alpha';
  Pmod.Lhood_ppt = log(sum(Pk))./prod(size(Y));
  Pmod.Pik = Pk ./ sum(Pk); 
  [trash, Pmod.C] = max(Pmod.Pik,[],2);
end


% Simply return if no prediction is requested
Yhat = [];
if (isempty(pt))
  return;
end

% Generate prediction at pt
if (maxx)
  [trash, k] = max(Pmod.Pik);
  Yhat = permute(M.Mu(pt+1,k,:),[1 3 2]);
else
  Yhat = Pmod.Pik*permute(M.Mu(pt+1,:,:),[2 3 1]);
end






