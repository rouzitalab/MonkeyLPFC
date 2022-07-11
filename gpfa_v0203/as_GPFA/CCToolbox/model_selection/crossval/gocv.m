function  [cvout,trajs,Krange,OrderRange] = gocv(trajs,ops)
%GOCV   Simple function to execute cross validation with curve data.
%
%   CV_Out = GOCV(Trajs,[Ops])
%
%   [CV_Out, trajs, Krange, OrderRange] = GOCV(...)

% Scott J Gaffney   5 February 2003
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'gocv';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


ops.Krange      = [1 2 3];
ops.OrderRange  = [2];
ops.NumRuns     = 10;
ops.NumEMStarts = 5;
ops.IterLimit   = 30;
ops.method      = 'lrm';
ops.zero        = 'zero';
ops.back        = 'noback';
ops             = SetFieldDef(ops,'State',[]);
ops.R.Diagonal  = 1;

if (~( strcmp(ops.method,'lrm') | strcmp(ops.method,'srm') ))
  ops.R.Diagonal  = 0;
  ops.R.Share     = 0;
end


% run lrm_ta
ops.method = 'lrm_ta';
ltacv = crossval(trajs,ops);
save lttcv_results lttcv;
clear lttcv;
end



if (0)
% run lrm
ops.method = 'lrm';
lcv = crossval(trajs,ops);
save lcv_results lcv;
folds = lcv.folds;
clear lcv;

% run lrm_at
ops.method = 'lrm_at';
latcv = crossval(trajs,ops);
save latcv_results latcv;
clear latcv;

% run lrm_aa
ops.method = 'lrm_aa';
laacv = crossval(trajs,ops);
save laacv_results laacv;
clear laacv;

% run lrm_tt
ops.method = 'lrm_tt';
lttcv = crossval(trajs,ops);
save lttcv_results lttcv;
clear lttcv;
end

cvout=[];

% Krange = ops.Krange;
% OrderRange = ops.OrderRange;