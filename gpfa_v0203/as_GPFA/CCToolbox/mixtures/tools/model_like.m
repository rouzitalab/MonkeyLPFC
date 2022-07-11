function [Lhood,other] = model_like(M,trajs,Ops)
%MODEL_LIKE  Calculate log-likelihood for any CCToolbox model.
%
%   [Lhood,other] = MODEL_LIKE(M,trajs,[Options])
%
%   Options
%     .NumSamps : number of samples to take for initial integration
%     .Scale    : (default 1.5) scale factor multiplied by NumSamps to 
%                 get the number of samples taken each time extra integration 
%                 attempts are taken. Note that the original samples are saved 
%                 from each failed attempt so that if .Scale is 1, 
%                 there will be exactly NumSamps samples taken on 
%                 the next attempt which gives a total of 2*NumSamps 
%                 samples (due to the saved samples on the first iteration).
%     .VarScale : (default 1) scale factor multiplied by the variance
%                 terms in the sampled distributions on extra integrations.
%                 This is a trick which allows for a wider range of sampling
%                 when an outlier is causing integration problems. The default
%                 is to not scale the variance terms at all. If you scale
%                 the variance, then you introduce a small bias.
%     .MaxAtts  : (default 3) this gives the maximum attempts that are
%                 tried during failed integrations.
%     .MsgBox
%       'nogui' - will not display progress box
%       else    - (default) progress box will be displayed

% Scott Gaffney   10 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'model_like';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

Lhood=[]; other=[];

% Handle invalid model
if (isnan(M.TrainLhood)),  return; end

if (nargin<3)
  Ops = [];
end
Ops = SetFieldDef(Ops,'NumSamps',[]);
Ops = SetFieldDef(Ops,'Scale',[]);
Ops = SetFieldDef(Ops,'VarScale',[]);
Ops = SetFieldDef(Ops,'MaxAtts',[]);
Ops = SetFieldDef(Ops,'MsgBox',[]);


% calculate the lhood
CalcLike = listmodels(M.method,'like');
if (~isempty(CalcLike))
  [Lhood,other] = feval(CalcLike,M,trajs,Ops);
end
