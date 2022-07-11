function [p,cinx] = colorplot(x,y,C,ops)
%COLORPLOT  General replacement for plot().
%
%   P = COLORPLOT(X,Y,[C],[Options])
%
%   P = COLORPLOT(X,Y,[C],linespec)
%
%   ColorString = COLORPLOT('colors')
%
%   Color = COLORPLOT(ClassIndex,[Options])
%
%   COLORPLOT('darken',h,spec) will darken recognized colors using
%   its internal infrastructure. h gives the handles to the objects
%   to darken and spec is one of the letters from the call colorplot('colors')
%   that determines the darkened color.
%
%   ops
%     .PlotLines
%     .PlotSymbols
%     .PlotColors
%     .linespec
%     .symbolspec
%     .colorspec
%     .MarkerSize
%     .Colors
%     .Symbols
%     .Styles
%

   
% Scott J Gaffney   5 October 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'colorplot';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


% COLORS = 'rgbkmcy';
COLORS = 'rgbkmc';
STYLES = {'-','-.','--',':'};
SYMBOLS = 'xosdv^<>*+ph';


%%% Begin Argument Processing
%%
if (nargin==2 & isstruct(y))
  ops = cexist('y',[]);
else
  ops = cexist('ops',[]);
end
if (~isempty(ops) & isstr(ops))
  ops.linespec = ops;
end
ops = SetFieldDef(ops,'linespec',[]);
ops = SetFieldDef(ops,'colorspec',[]);
ops = SetFieldDef(ops,'symbolspec',[]);
ops = SetFieldDef(ops,'PlotSymbols',1);
ops = SetFieldDef(ops,'PlotLines',1);
ops = SetFieldDef(ops,'PlotColors',1);
ops = SetFieldDef(ops,'PlotDayMarkers',0);
ops = SetFieldDef(ops,'MarkerSize',[]);
ops = SetFieldDef(ops,'PlotMarkerFaceColor',0);
ops = SetFieldDef(ops,'Colors',COLORS);
ops = SetFieldDef(ops,'Symbols',SYMBOLS);
ops = SetFieldDef(ops,'Styles',STYLES);
COLORS = ops.Colors;
SYMBOLS = ops.Symbols;
STYLES = ops.Styles;

if ((nargin==1 & ~isstr(x)) | (nargin==2 & isstruct(y)))
  p = MakeSpec(ops,x,COLORS,STYLES,SYMBOLS);
  return;
end
if (nargin==1 & isstr(x))  % 'colors'
  p = COLORS;  
  return;
end
if (nargin==3 & isstr(x))  % 'darken'
  p = DarkenColors(y,C);
  return;
end
if (exist('C')~=1 | isempty(C))
  C = 1;
end
% if (any(C > length(COLORS)))
%   warning(['Color label is larger than the number of available colors; ', ...
%             're-using colors. Colors will not be unique.']);
% end
%%
%%% End Argument Processing



p=[];
gcf;
hold on;
lenC = length(C);

% If we aren't plotting lines and we have one class, then
% duplicate it for all of the points.
if (~ops.PlotLines & lenC==1)
  lenC = length(x);
  C = ones(lenC,1)*C(1);
end

% Plot the points with one call
if (length(x)>1 & lenC==1)
  spec = MakeSpec(ops,C,COLORS,STYLES,SYMBOLS);
  p = plot(x,y,spec);
  
  if (ops.PlotDayMarkers)
    x_days = x(1:4:end);
    y_days = y(1:4:end);
    cops = ops;  cops.PlotSymbols=1;  cops.symbolspec='o';  cops.PlotLines=0;
    spec = MakeSpec(cops,C,COLORS,STYLES,SYMBOLS);
    p(end+1:end+length(x_days)) = plot(x_days,y_days,spec);
  end
  
  c = DarkenColors(p,spec);
  if (ops.PlotMarkerFaceColor)
    set(p,'MarkerFaceColor',c);
  end

% Or plot the points separately
else
  for i=1:length(x)
    spec = MakeSpec(ops,C(i),COLORS,STYLES,SYMBOLS);
    p(end+1) = plot(x(i),y(i),spec);
    c = DarkenColors(p(end),spec);  
    if (ops.PlotMarkerFaceColor)
      set(p(end),'MarkerFaceColor',c);
    end
  end
end

if (~isempty(ops.MarkerSize))
  set(p,'MarkerSize',ops.MarkerSize);
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MakeSpec
%
function [spec,ops] = MakeSpec(ops,c,COLORS,STYLES,SYMBOLS)
if (ops.PlotColors)
  if (isempty(ops.colorspec))
    ops.colorspec  = COLORS(mod(c-1,length(COLORS))+1);
  end
else
  ops.colorspec = [];
end
if (ops.PlotSymbols)
  if (isempty(ops.symbolspec))
    ops.symbolspec = SYMBOLS(mod(c-1,length(SYMBOLS))+1);
  end
else
  ops.symbolspec = [];
end
if (ops.PlotLines)
  if (isempty(ops.linespec))
    ops.linespec = STYLES{mod(c-1,length(STYLES))+1};
  end
else
  ops.linespec = [];
end

spec = [ops.colorspec, ops.linespec, ops.symbolspec];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DarkenColors
%
function c = DarkenColors(h,spec)
DarkGreen = [.0 .8 .4];
DarkCyan  = [.0 .0 .4];
if (~isempty(findstr(spec,'g')))
  c = DarkGreen;
  set(h,'color',c);
elseif (~isempty(findstr(spec,'c')))
  c = DarkCyan;
  set(h,'color',c);
else
  c = get(h(1),'color');
end
