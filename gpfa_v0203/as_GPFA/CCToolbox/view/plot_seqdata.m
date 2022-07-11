function h = plot_seqdata(x,Y,Seq,C,ops)
%PLOT_SEQDATA  Plot sequence data using the SEQ vector.
%
%   H = PLOT_SEQDATA_AA([x],Y,Seq,Mu,e,[C],[Options])
%     x,Y,Seq - same as returned by Trajs2Seq()
%     Mu      - model regression coefficients
%     e       - model parameters for scaling in space (see, e.g., lrm_cd)
%     C       - class labels (scalar, vector, or matrix)
%     Options - standard options structure
%
%   If x is empty, then an adaptive x will be used.
%
%   Options
%   -------------------------
%   SelectClass - {int}, only plots sequences in class SelectClass (default: all)
%   PlotTerminus  - {0,1}, plot a circle at the end point? (no)
%
% Comments:
%  - uses colorplot; passes ops through.

% Scott Gaffney   11 January 2002
% Department of Information and Computer Science
% University of California, Irvine

% Changes
% -------------------------------
% 28 January 2003 (Scott Gaffney)
%  - Changed the time axis to 0:len-1 from 1:len.
%  - Revamped the whole file, including parameters. It might not be compatible 
%    with old sister-code.


PROGNAME = 'plot_seqdata';
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

%% Handle C parameter
% if C is a scalar, then use for each sequence
if (prod(size(C))==1)
  C = ones(n,1)*C;
% else if C is a matrix, then assume it is a membership matrix
elseif (prod(size(C))~=length(C))
  [trash,C] = max(C,[],2);
end


%% Perform the plotting
sops = ops;  sops.PlotSymbols=1;  sops.symbolspec='o'; sops.PlotLines=0;
if (~isempty(ops.TerminusMarkerSize))
  sops.MarkerSize = ops.TerminusMarkerSize;
end
held = ishold;
h = [];
for i=1:n
  if (isempty(ops.SelectClass) | any(ops.SelectClass==C(i))) % plot this class?
    if (isempty(ops.SelectTrajs) | any(ops.SelectTrajs==i))    % plot this traj?
      indx = Seq(i):Seq(i+1)-1;
      if (isempty(x)) % use adaptive x?
        t = (0:length(Y(indx,1))-1)';
      else
        t = x(indx);
      end
      h{end+1} = colorplot(t,Y(indx,:),C(i),ops);
      if (ops.PlotTerminus),  colorplot(t(end),Y(indx(end),:),C(i),sops); end
      set(h{end},'Tag',num2str(i));
      hold on;
    end
  end
end

if (~held)
  hold off;
end