function s = vec2str(x)
%VEC2STR  Convert entries of a vector into cell array of strings
%   S = VEC2STR(x)

% Scott Gaffney   10 February 2004
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'vec2str';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


for i=1:length(x)
  s{i} = num2str(x(i));
end