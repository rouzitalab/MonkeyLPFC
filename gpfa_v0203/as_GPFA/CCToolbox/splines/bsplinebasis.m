function S = bsplinebasis(knots,K,X,Dn,output)
%BSPLINEBASIS  Calculate the B-spline basis functions at various X values.
%   S = BSPLINEBASIS(KNOTS,K,X,DERIV_ORDER,OUTPUT)
%     - KNOTS  : knots for spline
%     - K      : order of B-splines
%     - X      : values at which to evaluate the basis functions
%     - DERIV  : order of derivative (0=none,1=first_deriv,2=second_deriv,etc.)
%
%   If OUTPUT equals:
%   -------------------------
%   'block'    : S will be in block diagonal form.
%   'sparse'   : S will be a sparse matrix.
%   Otherwise  : S will be a standard matrix.

% Scott J Gaffney   29 January 2002
% Department of Information and Computer Science
% University of California, Irvine.
%
% This work is psuedo-derivative of several different other works
% including companion code from de Boor, and Ramsay and Silverman.
%
% de Boor, C. (1978). "A practical guide to splines."
%    New York, NY: Springer-Verlag.
%
% Ramsay, J., & Silverman, B. W. (1997). "Functional data analysis." 
%    New York, NY: Springer-Verlag.
%
%
% Changes
% -----------------
% 12 September 2003 (SJG)
%  - added derivative capability


PROGNAME = 'bsplinebasis';
if (~nargin)
  help(PROGNAME);
  return;
end

% Constants
STANDARD    = 1;
BLOCK       = 2;
SPARSE      = 3;

BLOCK_STR    = 'block';
SPARSE_STR   = 'sparse';


% Handle Preprocessing
S=[];
if (~isempty(find(diff(knots)<0)))
   errorbox('Argument Error: KNOTS must be nondecreasing.',PROGNAME);
   return;
end

M = length(knots);
L = M-K;
if (L < 1)
  errorbox('There are no B-splines for the given input.',PROGNAME)
  return;
end

Dn = cexist('Dn',0)+1; % 1 plus the order of derivative

if (exist('output')~=1 | isempty(output))
  output=STANDARD;
else
  switch (output)
  case BLOCK_STR
    output = BLOCK;
  case SPARSE_STR
    output = SPARSE;
  otherwise
    output = STANDARD;
  end
end



% Add end knots to satisfy end conditions
N = length(X);
X = X(:);
[addknots,NewLeft] = addendpts(knots,K); 
addknots = addknots(:);
addlen = length(addknots)-K;

% Find the knot to the left of each point
knt = max(K,findfloor(addknots(1:addlen),X));



% Perform the recurrence for all X.
B = zeros(Dn*N,K); % B comes from Bik = w_{i,k}*B_{i,k-1} + w'_{i,k}*B_{i+1,k-1}
                   % where w'_{i,k} = 1 - w_{i+1,k}.
                   % Note there is at most K supporting B-splines for any 
                   % of the N x's. Thus we need not concern ourselves with 
                   % all of the L B-splines for each x. Instead we just 
                   % calculate the K nonzero B-splines for each x. 
                   % Hence B is of size N-by-K. If we have derivatives then we
                   % need to increase N appropriately using Dn.
dindx = Dn*[1:N]';  % this gives us access to the requested (deriv) values

B(:,1) = ones(Dn*N,1);  % set Bik to 1 for k=1 since B_{i,1}==1 for nonzero B
% Calculate the w_ik and B_ik starting from order 2 on up to K
for (k=1:K-Dn)        % note that k==j means calculate B_{i,j+1}
  last = zeros(N,1);
  for (i=1:k)  % calculate supporting w_ik from far left to right
    inRight  = addknots(knt+i)   - X;          % w'_{knt-k+i-1,k+1}
    inLeft   = X - addknots(knt+i-k);          % w_{knt-k+i,k+1}
    prevB    = B(dindx,i)./(inRight+inLeft);  % B_{knt-k+i,k}/D, D gives norm
    B(dindx,i) = last + (inRight.*prevB);     % B_{knt-k+i-1,k+1} = ...
    last(:)  = inLeft.*prevB;                 % w*Bk
  end
  B(dindx,k+1) = last; % initialization for the k+2nd order B-splines or
                       % finish the last column for the k+1st order B-splines
end


% Slide values down (left) in the matrix as we sweep through again.
% This is only necessary for derivative values.
for (dk=1:Dn-1)
  k = K-Dn+dk;
  indx = dindx-1;
  last = zeros(N,1);
  for (i=1:k)  % calculate supporting w_ik from far left to right
    inRight  = addknots(knt+i)   - X;          % w'_{knt-k+i-1,k+1}
    inLeft   = X - addknots(knt+i-k);          % w_{knt-k+i,k+1}
    prevB    = B(dindx,i)./(inRight+inLeft);  % B_{knt-k+i,k}/D, D gives norm
    B(indx,i) = last + (inRight.*prevB);      % save in previous block
    last(:)  = inLeft.*prevB;                 % w*Bk
  end
  B(indx,k+1) = last;  % initialization for the k+2nd order B-splines or ...
  dindx = indx;   % keep going down (left)
end


% Now we can calculate the derivatives (via differencing)
for dk=Dn-1:-1:1
   k = K-dk;
   dindx = (dk:Dn-1)'*ones(1,N) + ones(Dn-dk,1)*indx';  % repmat
   for i=k:-1:1
      oper = ones(Dn-dk,1)*(addknots(knt+i)'-addknots(knt+i-k)')./k; % repmat
      B(dindx,i)   = -B(dindx,i)./oper(:);
      B(dindx,i+1) =  B(dindx,i+1) - B(dindx,i);
   end
end




% Catch points outside support region
nosupport = find(X<knots(1) | X>knots(M));
if (~isempty(nosupport))
  % build index that gets rid of successive entries created by Dn
 indx = (1-Dn:0)'*ones(1,length(nosupport))+ Dn*ones(Dn,1)*nosupport(:)';
 B(indx(:),:) = zeros(Dn*length(nosupport),K);
end

% Adjust index to deal with a negative offset.
knt = knt + max(0,-NewLeft);


% Handle output issues
switch (output)
  
case BLOCK
  [BB,nb,rows,last] = makeblock(B(Dn*(1:N),:),knt,NewLeft,addlen,N,K,L);
  S = [41 nb rows K last L-sum(last) BB(:)'];
  
case SPARSE
  [BB,nb,rows,last] = makeblock(B(Dn*(1:N),:),knt,NewLeft,addlen,N,K,L);
  S = makesparse(BB,nb,rows,last,N,K,L);
  
otherwise
  S = makestandard(B(Dn*(1:N),:),knt,NewLeft,addlen,N,K,L);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% makeblock
%%
function [R,nb,rows,last] = makeblock(R,leftknot,NewLeft,addlen,N,k,L)
if (NewLeft>0)
  width = 2*k-1;
  cc = zeros(N*width,1);
  knotsIndx = min(k,leftknot-NewLeft); 
  cc(repmat(([1-N:0]+N*knotsIndx).',1,k) + repmat(N*[0:k-1],N,1)) = R;
  R(:) = cc(repmat([1-N:0].',1,k) + repmat(N*(k+[0:k-1]),N,1));
  leftknot = leftknot + k - knotsIndx;
end
dl = diff(leftknot);
index = [0 find(dl>0) N];
rows=diff(index);
nb=length(index)-1;
last=dl(index(2:nb));
if (NewLeft<0)
  nb=nb+1; 
  rows = [0 rows]; 
  last = [-NewLeft last];
end
addr = addlen-L-NewLeft;
if (addr<0)
  nb=nb+1;
  rows=[rows 0];
  last=[last -addr];
end

  

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% makesparse
%%
function S = makesparse(R,nb,rows,last,N,k,L)
nr = (1:N).'; 
nc = 1:k; 
nrnc = N*k;
ncc = zeros(1,N);
ncc(1+cumsum(rows(1:(nb-1)))) = last;
ncc = reshape(cumsum(ncc),N,1);
ijs = [reshape(nr(:,ones(1,k)),nrnc,1), ...
    reshape(ncc(:,ones(1,k))+nc(ones(N,1),:),nrnc,1), reshape(R,nrnc,1)];
index = find(ijs(:,2)>L);
if (~isempty(index))
  ijs(index,:) = [];
end
S = sparse(ijs(:,1),ijs(:,2),ijs(:,3),N,L);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% makestandard
%%
function S = makestandard(B,leftknot,NewLeft,addlen,N,k,L)
width = max([L,addlen])+ 2*k-2;

D = repmat([1-N:0].',1,k)+ repmat(N*leftknot.',1,k) + repmat(N*[-k+1:0],N,1);
c = zeros(N*width,1);
c(D) = B;  % c(D)(i,j) == c(D(i,j))

D = (1-N:0)'*ones(1,L) + N*ones(N,1)*(max(0,NewLeft)+[1:L]);
S = reshape(c(D),N,L);

  
  
  
  