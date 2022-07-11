function [Lhood,other] = srm_like(varargin)
%SRM_LIKE  Calculate log-likelihood with SRM model.
%
%   [Lhood,Other] = SRM_LIKE(M,Trajs,[Options])
%    - M       : trained model
%    - Trajs   : 'Trajs' structure; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   [Lhood,Other] = SRM_LIKE(M,X,Y,Seq,[Options])
%    - M       : trained model
%    - X,Y,Seq : curves in Sequence format; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   Other
%     .C   : classification labels

% Scott Gaffney   15 February 1999
% DataLab@UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.

PROGNAME = 'srm_like';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%%
args = varargin; clear varargin;
n = length(args);
trajs=[]; X=[]; Y=[]; Seq=[]; Ops=[];
%
% Check for calling convention
%
% SRM_LIKE(M,Trajs,[Options])
if (n<4)
  M = args{1};
  trajs = args{2};
  if (n>2)
    Ops = args{3};
  end
  
% SRM_LIKE(M,x,Y,Seq,[Options])
else
  M = args{1};
  X = args{2};
  Y = args{3};
  Seq = args{4};
  if (n>4)
    Ops = args{5};
  end
end
%%
%%% End Argument Processing

if (isempty(Y))
  [Y,X,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
end
if (size(X,2)==1)
  X = bsplinebasis(M.knots,M.order,X);
end

func = listmodels('lrm','like');
[Lhood,other] = feval(func,M,X,Y,Seq,Ops);