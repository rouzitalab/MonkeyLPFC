function n = trajslen(trajs)
%TrajsLen  Calculate number of curves/trajectories in 'Trajs' structure
%
%   n = TrajsLen(trajs)

% Scott J Gaffney   20 October 2003
% Department of Information and Computer Science
% University of California, Irvine.
%
% Changes
% -----------------------------------


PROGNAME = 'trajslen';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

n=0;

% trajs is already in Cell format
if (iscell(trajs))
  n = length(trajs);
  return;
end

% Feature Vector format
if (~isstruct(trajs) & isnumeric(trajs))
  n = size(trajs,1);
  return;
end

% Handle 'Trajs' structure
if (isfield(trajs,'Y'))
  if (~isfield(trajs,'Seq') | isempty(trajs.Seq))
    n = trajslen(trajs.Y);
  else
    n = length(Seq)-1;
  end
  return;
end

% Cycs format
if (isfield(trajs,'trajs'))
  n = length(trajs.trajs);
  return;
end