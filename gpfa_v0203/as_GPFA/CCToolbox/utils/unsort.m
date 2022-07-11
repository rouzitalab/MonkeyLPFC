function index = unsort(sin)
%UNSORT  Unshuffles a sorted index so that you can unsort a sorted array.
%   INDEX = UNSORT(SIN) returns the sorted index SIN in a form
%   useful to unsort the associated sorted array.
%
%   Example:
%
%   X = rand(1,10);
%   [SortedX,SortedIndex] = sort(X);
%   UnsortedIndex = unsort(SortedIndex);
%   Y = SortedX(UnsortedIndex);
%   [X;Y]

% Scott J Gaffney   29 January 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'unsort';
if (~nargin)
  help(PROGNAME);
  return;
end

[trash,index] = sort(sin);