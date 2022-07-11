function s = ssum(x,seq,m)
%SSUM  Sequence sum
%   S = SSUM(X,Seq,[MaxLenSeq]) returns the sequence sum 
%   in S for vectors X and Seq. Seq is in the same form as that returned 
%   in Trajs2Seq(). You can provide the maximum sequence length in 
%   MaxLenSeq if you have it handy.


lens = diff(seq);
n = length(lens);
if (nargin<3)
  m = max(lens);
end

j = ones(1,seq(end)-1+n-1);
j(seq(2:end-1)+(0:n-2)) = -lens(1:end-1);
j(:) = cumsum(j);
j(seq(2:end-1)+(0:n-2)) = [];

i = zeros(1,seq(end)-1);
i(seq(2:end-1)) = m;
i(:) = cumsum(i) + j;

%p = full(sum(sparse(i,j,x,n,m,length(x)),2));
p = zeros(m,n);
% p(sub2ind([n,m],i,j)) = x;
p(i) = x;
s = sum(p',2);