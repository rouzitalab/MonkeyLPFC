function test = usage(name,str,test)
%USAGE  Outputs usage message in standard format
%   U = USAGE(NAME,STR,[TEST]) returns 1 if a usage message
%   was output, otherwise returns 0. NAME is the name of the
%   function for whom we output the usage message in STR. STR
%   may contain any input accepted by FPRINTF.
%
%   If TEST is provided, then we only output a usage message
%   if TEST is not zero.

% Scott J Gaffney   11 October 2001
% Department of Information and Computer Science
% University of California, Irvine.


if (exist('test') ~= 1 | isempty(test))
  test=1;
end

if (test)
  fprintf(['\nUSAGE: %s', str, '\n\n'],name);
end