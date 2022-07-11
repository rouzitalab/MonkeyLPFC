function ah = gaa(fh)
%GAA  Get all of the axes handles in a figure.
%   ah = GAA([fh])
%
%   Returns all of the axes handles that are in the children list
%   for the figure fh. Uses the current figure if fh is not provided.
%

% Scott Gaffney   29 January 2003
% Department of Information and Computer Science
% University of California, Irvine

fh = cexist('fh',gcf);

ah = findobj(get(fh,'children'),'type','axes');