function [y,x] = seq2cell(Y,Seq,X)
%Seq2Cell  Convert curves in Sequence form to Cell form
%
%   [y,x] = Seq2Cell(Y,Seq,[X])

% Scott J Gaffney   20 October 2003
% Department of Information and Computer Science
% University of California, Irvine.
%
% Changes
% -----------------------------------


PROGNAME = 'Seq2Cell';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

x=[];
X = cexist('X',[]);
n = length(Seq)-1;
y = cell(1,n);
if (nargout>1)
  x = cell(1,n);
end

if (isempty(x))
  for i=1:n
    y{i} = Y(Seq(i):Seq(i+1)-1,:);
  end
else
  for i=1:n
    y{i} = Y(Seq(i):Seq(i+1)-1,:);
    if (isempty(X))
      ni = Seq(i+1)-Seq(i);
      x{i} = (0:ni-1)';
    else
      x{i} = X(Seq(i):Seq(i+1)-1);
    end
  end
end

