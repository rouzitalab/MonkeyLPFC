function h = plot_seqdata_aa(x,X,Y,Seq,Mu,e,C,ops)
%PLOT_SEQDATA_AA  Plot sequence data using the Seq vector.
%
%   H = PLOT_SEQDATA_AA([x],X,Y,Seq,Mu,e,[C],[Options])
%     x,Y,Seq - same as returned by Trajs2Seq()
%     X       - regression/spline matrix
%     Mu      - model regression coefficients
%     e       - model parameters for scaling in space (see, e.g., lrm_cd)
%     C       - class labels (scalar, vector, or matrix)
%     Options - standard options structure
%
%   If x is empty, then an adaptive x will be used.
%   
%   Options
%   -------------------------
%   SelectClass - {int}, only plots sequences in class SelectClass 
%                        (default: all)
%   PlotTerminus  - {0,1}, plot a circle at the end point? (no)
%   TerminusMarkerSize  - {int}, size of marker (current default)
%
% Comments:
%  - uses colorplot; passes ops through.

% Scott Gaffney   11 January 2002
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -------------------------------


PROGNAME = 'plot_seqdata_aa';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

%% Set defaults
ops = cexist('ops',[]);
ops = SetFieldDef(ops,'SelectClass',[]);
ops = SetFieldDef(ops,'SelectTrajs',[]);
ops = SetFieldDef(ops,'PlotTerminus',0);
ops = SetFieldDef(ops,'TerminusMarkerSize',[]);
C = cexist('C',1);
n = length(Seq)-1;

% if C is a scalar, then use for each sequence
% if its a matrix then assume it is a membership matrix
if (prod(size(C))==1)
  C = ones(n,1)*C;
elseif (prod(size(C))~=length(C))
  [trash,C] = max(C,[],2);
end


%% Perform the plotting
sops = ops;  sops.PlotSymbols=1;  sops.symbolspec='o'; sops.PlotLines=0;
if (~isempty(ops.TerminusMarkerSize))
  sops.MarkerSize = ops.TerminusMarkerSize;
end
[N,P,shareX] = size(X);
[ne,shareK,shareE] = size(e);
D = size(Y,2);
held = ishold;
h=[];

for i=1:n
  if (isempty(ops.SelectClass) | any(ops.SelectClass==C(i))) % plot this class?
    if (isempty(ops.SelectTrajs) | any(ops.SelectTrajs==i))    % plot this traj?
      indx = Seq(i):Seq(i+1)-1;
      shk = min(shareK,C(i));
      ytil = Y(indx,1)+ (1-e(i,shk,1)).*(X(indx,:,1)*Mu(:,C(i),1));
      if (D==2)
        shx = min(shareX,2);
        she = min(shareE,2);
        xtil = ytil;
        ytil = Y(indx,2)+ (1-e(i,shk,she)).*(X(indx,:,shx)*Mu(:,C(i),2));
      else
        if (isempty(x))
          xtil = (0:length(indx)-1)';
        else
          xtil = x(indx);
        end
      end
      h{end+1} = colorplot(xtil,ytil,C(i),ops);
      set(h{end},'Tag',num2str(i));
      if (ops.PlotTerminus),  
        colorplot(xtil(end),ytil(end),C(i),sops); 
      end
      hold on;
    end
  end
end

if (~held)
  hold off;
end