function [s,scl] = ssprod(x,seq,m)
%SSPROD  Scaled sequence product (see SPROD for non-scaled version)
%   [P,SCL] = SSPROD(X,Seq,[MaxLenSeq]) returns the scaled sequence product 
%   in P and the scale factors in SCL for vector X and vector Seq. Seq is
%   in the same form as that returned in Trajs2Seq(). You can provide the 
%   maximum sequence length in MaxLenSeq if you have it handy. The
%   sequences are currently scaled by their MEAN.
%
%   EXAMPLES:
%   To get the log product you perform this set of commands.
%     [s,scl] = ssprod(x,seq);
%     logProd = log(s) + diff(seq(:)).*log(scl);

%   The actual product is available as follows.
%     [s,scl] = ssprod(x,seq);
%     actualProd = (scl.^diff(seq(:))).*s;
%
%   Of course you can also get the log product as...
%     sumLog = ssum(log(x),Seq);   % ...thus the above example is simply
%   for illustration purposes.

% Scott J Gaffney  2 September 2003
% School of Information and Computer Science
% University of California, Irvine


% perform basic SPROD initialization
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

p = ones(m,n);
p(i) = x;


% now we scale the columns before we take the product
scl = mean(p);
p = p./ (ones(m,1)*scl);
s = prod(p',2);
scl = scl.';