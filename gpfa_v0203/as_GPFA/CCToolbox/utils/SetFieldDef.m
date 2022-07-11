function [S,isSet] = SetFieldDef(S,FieldName,Value)
%SETFIELDDEF  Set default values in a structure.
%   [S,isSet] = SETFIELDDEF(S,FieldName,Value)

% Scott Gaffney   11 January 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'SetFieldDef';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

if (~isfield(S,FieldName) | isempty(getfield(S,FieldName)))
  S = setfield(S,FieldName,Value);
  isSet = 1;
end
