function [trajs,M] = rm_cd_b_sim(M,NumSamps,len,ops)
%RM_CD_B_SIM  Simulate from a *RM_CD_B model.
%   [TRAJS,M] = RM_CD_B_SIM(M,NUM_SAMPS,[LEN],[Ops])

% Scott J Gaffney   4 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'rm_cd_b_sim';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

len = cexist('len',[]);
ops = cexist('ops',[]);
ops = SetFieldDef(ops,'SampleAt',[]);
ops = SetFieldDef(ops,'SampMinLen',1);
ops = SetFieldDef(ops,'SampMaxLen',10);
ops = SetFieldDef(ops,'DoPerturb',0);
ops = SetFieldDef(ops,'Perturb_Std',.1);


% Handle specified sampling vector
if (~isempty(ops.SampleAt))
  lensamp = length(ops.SampleAt);
  ops.SampMaxLen = lensamp;
end

% Handle length specifier
if (isempty(len))
  len = floor(rand(NumSamps,1)*(ops.SampMaxLen+1-ops.SampMinLen)) ...
    + ops.SampMinLen;
elseif (length(len)==1)
  len = len(ones(NumSamps,1));
end

% Stratify the sample according to the priors
C = pmfrnd(M.Alpha, NumSamps);
C = sort(C);   C=C(:);  % sort them just because

% Set the sample range
if (isempty(ops.SampleAt))
  samp = (0:max(len)-1)'; 
else
  if (any(len>lensamp))
    fprintf('SIMRM: WARNING; invalid lengths were detected, reducing max.\n');
    len(find(len>lensamp)) = lensamp;
  end
  samp = ops.SampleAt(:);
end

% messy way to select the method
DoSRM = 0;
if (isfield(M,'method'))
  if (M.method(1)=='s')
    DoSRM = 1;
  end
elseif (isfield(M,'knots'))
  DoSRM = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform the sampling
[P,K,D] = size(M.Mu);
M.b = zeros(NumSamps,1);
M.e = zeros(NumSamps,D);
M.f = zeros(NumSamps,D);
if (ops.DoPerturb | ~isempty(ops.SampleAt))
  xx = zeros(sum(len),1);
  start=1;
end

for i=1:NumSamps
  ni = len(i);
  x = samp(1:ni);
  
  % handle random perturbation of the sampling value
if (ops.DoPerturb | ~isempty(ops.SampleAt))
    indx = start:start+ni-1;  start = start+ni;
    if (ops.DoUniformPerturb)
      x = rand(ni,1)*samp(ni);
    elseif (ops.DoPerturb)
      x = x + randn(ni,1)*ops.Perturb_Std;
    end
    x = sort(x);
    xx(indx) = x;
  end
  
  M.b(i,1) = randn(1).*sqrt(M.S(C(i)));      % b ~ N(0,s^2)
  M.e(i,:) = randn(1,D).*sqrt(M.U(C(i),:)) + 1;
  M.f(i,:) = randn(1,D).*sqrt(M.T(C(i),:));
  zero = find(M.e(i,:)<0);
  if (~isempty(zero)), M.e(i,zero)=-M.e(i,zero);  end
  for d=1:D
    Eps = randn(ni,1)*sqrt(M.Sigma(C(i),d));
    if (DoSRM)
      X = bsplinebasis(M.knots,M.order,x-M.b(i));
    else
      X = regmat(x-M.b(i),P-1);
    end
    trajs{i}(:,d) = M.e(i,d).*X*M.Mu(:,C(i),d) + M.f(i,d) + Eps;
  end
end
M.TrueC = C;
if (ops.DoPerturb | ~isempty(ops.SampleAt))
  M.x = xx;
end
