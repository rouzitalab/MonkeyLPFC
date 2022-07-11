function showmccv(mccv,Legend)
%SHOWMCCV  View MCCV clustering results
%   SHOWMCCV(mccv,[Legend])

% Scott J Gaffney   15 May 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'showmccv';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


ops.PlotLines = 0;
ops.PlotSymbols =1;
ops.symbolspec = '*';


if (exist('Legend')~=1 | isempty(Legend))
  no_legend = 1;
else
  no_legend = 0;
end

% calculate means and deviations
[mu,sigma] = mccv_like(mccv);
[muSSE,sigmaSSE] = mccv_SSE(mccv);

if (~isempty(mu))
  figure;
  orderHnd = [];
  K=mccv.Krange(:);
  Q = length(mccv.OrderRange);
  for q=1:Q
    spec = colorplot(q);
    spec = [spec(1), '-'];
    h = errorbar(K,mu(:,q),sigma(:,q),spec);
    set(h,'LineWidth',2);
    %h = colorplot(K,Pred_Error(:,i),i,ops);
    orderHnd(end+1) = h(1);
    hold on;
  end
  
  set(gca,'xtick',K);
  %set(gca,'xlim', [K(1)-1 K(end)+1]);
  title('MCCV Results');
  ylabel('Test log-likelihood');
  xlabel('Number of components');
  if (~no_legend)
    legend(orderHnd,Legend);
  end
end


if (~isempty(muSSE))
  figure;
  orderHnd = [];
  K=mccv.Krange(:);
  Q = length(mccv.OrderRange);
  for q=1:Q
    spec = colorplot(q);
    spec = [spec(1), '-'];
    h = errorbar(K,muSSE(:,q),sigmaSSE(:,q),spec);
    set(h,'LineWidth',2);
    %h = colorplot(K,Pred_Error(:,i),i,ops);
    orderHnd(end+1) = h(1);
    hold on;
  end
  
  set(gca,'xtick',K);
  %set(gca,'xlim', [K(1)-1 K(end)+1]);
  title('MCCV Results');
  ylabel('Test SSE');
  xlabel('Number of components');
  if (~no_legend)
    legend(orderHnd,Legend);
  end
end


