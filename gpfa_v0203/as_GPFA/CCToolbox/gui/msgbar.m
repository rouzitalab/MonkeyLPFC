function [hfig] = msgbar(hfig,msg,append)
%MSGBAR  Create/Update a message bar useful as a stripped-down WAITBAR.
%   HFig = MSGBAR(HFig,Msg,'app') 
%
%  INPUTS:
%    HFig   : Handle to figure window
%    Msg    : String message
%    Append : appends the string Msg to the message bar if Append='app'
%
%  OUTPUTS
%    HFig   : Handle to figure window

% Scott Gaffney   18 March 2002
% DataLab@UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.

PROGNAME = 'msgbar';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%% Begin Argument Processing
%
if (exist('hfig')~=1)
  hfig = [];
end
if (exist('msg')~=1 | isempty(msg))
  msg='';
end
%
%% End Argument Processing


if (isempty(hfig))
  %set(0,'Units','normalized');
  %ScreenSize = get(0,'ScreenSize');
  width = 0.2571;
  height = 0.0714;
  
  % Calculate a centered position
  %pos = [ScreenSize(3)/2-width/2 ScreenSize(4)/2-height/2 width height];
  pos = [0.2750  0.4857  width  height];
  
  hfig = figure(...
    'Units','normalized', ...
    'Position',pos, ...
    'Resize','off', ...
    'CreateFcn','', ...
    'NumberTitle','off', ...
    'IntegerHandle','off', ...
    'MenuBar','none', ...
    'Tag','SJG_MsgBar', ...
    'Interruptible','off', ...
    'Visible','on');
  
  haxes = axes('Visible','off','Units','normalized','Position',[0 0 1 1], ...
    'Parent',hfig);
  htext = text(.1,.6,'','Parent',haxes,'Tag','SJG_MsgText');
end

htext = findobj(hfig,'Tag','SJG_MsgText');
if (exist('Append')==1 & strcmp(Append,'app'))
  set(htext,'String', [get(htext,'String'),msg] , ...
    'FontName','Times','FontWeight','normal','FontSize',12);
else
  set(htext,'String',msg,'FontName','Times','FontWeight','normal','FontSize',12);
end

drawnow;


