function [y,x] = trajs2cell(trajs)
%Trajs2Cell  Extract curves in 'Trajs' format into Cell format
%
%   [y,x] = Trajs2Cell(trajs)

% Scott J Gaffney   20 October 2003
% Department of Information and Computer Science
% University of California, Irvine.
%
% Changes
% -----------------------------------


PROGNAME = 'Trajs2Cell';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

y=[];  x=[];

% trajs is already in Cell format
if (iscell(trajs))
  y = trajs;
  return;
end

% Feature Vector format
if (~isstruct(trajs) & isnumeric(trajs))
  y = fvec2cell(trajs);
  return;
end

% Handle 'Trajs' structure
if (isfield(trajs,'Y'))
  trajs = SetFieldDef(trajs,'X',[]);
  if (iscell(trajs.Y))
    y = trajs.Y; 
    x = trajs.X;
  end
  if (~isfield(trajs,'Seq') | isempty(trajs.Seq))
    y = fvec2cell(trajs.Y);
  else
    [y,x] = seq2cell(trajs.Y,trajs.Seq,trajs.X);
  end
  return;
end

% Cycs format
if (isfield(trajs,'trajs'))
  y = trajs.trajs;
end