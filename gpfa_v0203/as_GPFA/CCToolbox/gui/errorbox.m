function errorbox(msg,title)
%ERRORBOX  Display a modal dialog box that suspends execution
%   ERRORBOX(MSG,TITLE)

% Scott J Gaffney   5 October 2001
% Department of Information and Computer Science
% University of California, Irvine.


if (exist('msg')~=1)
  msg = '';
end
if (exist('title')~=1)
  title = '';
end

uiwait(msgbox(msg,title,'error','modal'));