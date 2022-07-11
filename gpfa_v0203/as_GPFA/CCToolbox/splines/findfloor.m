function flr = findfloor(grid, pts)
%FINDFLOOR  Finds the index of the floor of PTS(i) within the GRID.
%   FINDFLOOR(GRID,PTS)

% Scott J Gaffney   29 January 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'findfloor';
if (~nargin)
  help(PROGNAME);
  return;
end

if (any(diff(grid)<0))
  errorbox('Argument Error: GRID must be nondecreasing.',PROGNAME);
  return;
end

G = length(grid);
P = length(pts);
grid = grid(:)';
pts = pts(:)';

[trash,sin] = sort([grid pts]);
flr = find(sin > G);
flr = flr - (1:P);

[trash,pin] = sort(pts);
flr = flr(unsort(pin));