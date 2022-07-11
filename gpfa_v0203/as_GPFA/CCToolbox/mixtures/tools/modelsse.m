function [sse,Ops] = modelsse(M,trajs,Ops)
%MODELSSE  Calculate Sum of Squared Error (SSE) over test set of curves
%
%   [SSE,Ops] = MODELSSE(M,trajs,[Options]) when .Method='sse'
%
%   PredModel = MODELSSE(M,trajs,[Options]) when .Method='pred'
%
%   Options
%     .Method
%       'sse'  - (default) return Sum of Squared Error (per point)
%       'pred' - return hidden variable predictions given all curves;
%                this overrides all other options!
%
%     .PredType
%       'fwbk' - (default) predict each point given all others (forward and 
%                backward)
%       'forw' - predict points using only past (forward)
%
%     .Points  (note: this works in conjunction with .PredType)
%       'all'  - (default) predict all points
%       'half' - predict last half given first half
%       'last' - predict last .NumPredPoints points, given rest of curve
%       'rand' - predict .NumPredPoints randomly chosen points
%
%     .NumPredPoints: (default 1) max number of points predicted for EACH curve
%     .RandPoints: n-by-1 cell matrix, each cell gives a vector of indicies
%                  referencing the points in the i-th curve that are to be
%                  predicted. If this is passed-in then these points
%                  are used when .Points='rand'. This matrix is provided
%                  on output in Ops so that the same random points can be
%                  used again.
%
%     .Membership
%       'sum-out' - (default) sum prediction across Pik (membership values)
%       'max'     - predict according to class with largest Pik (classify)
%
%     .MsgBox
%       'nogui' - will not display progress box
%       else    - (default) progress box will be displayed
%
%     .ShowPrediction
%       1  - graphically display prediciton results
%       0  - (default) no graphics output
%
%     .predfig - you may provide a figure handle for .ShowPrediction if
%                so desired (not required).

% Scott Gaffney   10 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'modelsse';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

% Constants
SSE    = 1;
PRED   = 2;
FORW   = 1;
FWBK   = 2;
ALL    = 1;
HALF   = 2;
LAST   = 3;
RAND   = 4;
SUM    = 1;
MAX    = 2;
GUI    = 1;
NOGUI  = 2;


%%% Handle Argument Processing
%%%
if (nargin<3)
  Ops = [];
end
Ops = SetFieldDef(Ops,'Method','sse');
Ops = SetFieldDef(Ops,'PredType','fwbk');
Ops = SetFieldDef(Ops,'Points','all');
Ops = SetFieldDef(Ops,'NumPredPoints',1);
Ops = SetFieldDef(Ops,'RandPoints',[]);
Ops = SetFieldDef(Ops,'Membership','sum-out');
Ops = SetFieldDef(Ops,'MsgBox',1);
Ops = SetFieldDef(Ops,'ShowPrediction',0);
Ops = SetFieldDef(Ops,'predfig',[]);

% Convert to integer options
switch (Ops.Method)
  case 'pred'
    Method = PRED;
  otherwise
    Ops.Method = 'sse';
    Method = SSE;
end
switch (Ops.PredType)
  case 'forw'
    PredType = FORW;
  otherwise
    Ops.PredType = 'fwbk';
    PredType = FWBK;
end
switch (Ops.Points)
  case 'half'
    Points = HALF;
  case 'last'
    Points = LAST;
  case 'rand'
    Points = RAND;
  otherwise
    Points = ALL;
    Ops.Points = 'all';
end
switch (Ops.Membership)
  case 'max'
    Membership = MAX;
  otherwise
    Membership = SUM;
    Ops.Membership = 'sum';
end
switch (Ops.MsgBox)
  case 'nogui'
    MsgBox = NOGUI;
  otherwise
    MsgBox = GUI;
    Ops.MsgBox = 'gui';
end
%%%
%%% End Argument Processing

% Handle invalid model issues
if (isfield(M,'TrainLhood') & isnan(M.TrainLhood))
  sse = NaN;
  if (Method==PRED),  sse = M; end
  return;
end

% Preprocessing
Predict = str2func(listmodels(M.method,'pred'));   % get function handle

% build data
[Y,x,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);

% handle Posterior prediction right away and get out of here
if (Method==PRED)
  [trash, sse] = feval(Predict,M,[],x,Y,Seq);
  return;
end

% Graphics
if (MsgBox==GUI)
  MsgHnd = msgbar([],'');
end
isCreated = 0;
if (Ops.ShowPrediction)
  fpos = [0.0029    0.6038    0.3329    0.3267];
  if (isempty(Ops.predfig))
    isCreated = 1;
    Ops.predfig = figure('units','normalized','position',fpos);
  end
end
  
% Loop over each curve
n = length(Seq)-1;
if (isempty(Ops.RandPoints)),  Ops.RandPoints = cell(1,n); end
sse=0;
numpts=0;
for i=1:n
  indx = Seq(i):Seq(i+1)-1;
  ni = length(indx);
  if (MsgBox==GUI),  msgbar(MsgHnd,sprintf('Curve %d\n',i)); end
  if (Ops.ShowPrediction)
    figure(Ops.predfig);  cla;
    plot(x(indx),Y(indx,1),'b:');  hold on;
  end

  % select prediction points
  if (Points==ALL)
    ppts = 1:ni;
  elseif (Points==HALF)
    ppts = ceil(ni/2):ni;
  elseif (Points==LAST)
    MinNum = min(ni,Ops.NumPredPoints);
    ppts = ni-MinNum+1:ni;
  else   % RAND
    if (isempty(Ops.RandPoints{i}))
      ppts = randperm(ni);
      ppts = ppts(1:min(ni,Ops.NumPredPoints));
      Ops.RandPoints{i} = ppts;
    else
      ppts = Ops.RandPoints{i};
    end
    ppts = ppts(:)';
  end
  
  % cycle through the prediction points in this curve
  numpts = numpts + length(ppts);
  for j=ppts
    if (PredType==FORW)
      ndx = indx(1:j-1);
    else  % FWBK
      ndx = indx;
      ndx(j) = [];
    end
    
    % predict hidden data
    Yhat = feval(Predict,M,x(indx(j)),x(ndx),Y(ndx,:),[1 length(ndx)+1]);
    
    if (Ops.ShowPrediction)
      plot(x(indx(j)),Yhat(1),'rx');
    end
    
    % and record error in prediction
    sse = sse + sum((Y(indx(j),:)-Yhat).^2);
  end
end

sse = sse ./numpts;
Ops.TotalNumPts = numpts;

if (MsgBox==GUI)
  delete(MsgHnd);
end
if (isCreated)
  delete(Ops.predfig);
end

