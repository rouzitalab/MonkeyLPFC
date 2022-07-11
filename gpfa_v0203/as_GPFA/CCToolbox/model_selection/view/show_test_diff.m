function table = show_test_diff(mccv,ops)
%SHOW_TEST_DIFF
%   TableOfValues = SHOW_TEST_DIFF(mccv,[ops])

% Scott J. Gaffney  5 November 2003
% School of Information and Computer Science
% University of California, Irvine

PROGNAME = 'show_test_diff';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

% get the methods
methods = fieldnames(mccv);
if (isempty(methods)),  return; end
NumMethods = length(methods);

% grab the options
indK = [];  K = [];
ops = cexist('ops',[]);
if (~isstruct(ops) & isvector(ops))
  indK = ops;  ops=[];
end
if (isempty(ops))
  obj = getfield(mccv,methods{1});
  if (isfield(obj,'runs'))
    ops = rmfield(obj,'runs');
  elseif (isfield(obj,'ops'))
    ops = obj.ops;
  end
end

ops = SetFieldDef(ops,'PlotSSE',1);
ops = SetFieldDef(ops,'PlotLike',1);
ops = SetFieldDef(ops,'PlotBIC',1);
ops = SetFieldDef(ops,'MarkerSize',10);
ops = SetFieldDef(ops,'FontSize',14);
ops = SetFieldDef(ops,'LineWidth',.5);

% Set the K plot points
if (isfield(ops,'Krange'))
  K = ops.Krange';
else
  % K wasn't specified, so we must calculate it from the first data field
  obj = getfield(mccv,methods{1});
  if (isfield(obj,'runs'))
    obj = getfield(mccv,methods{1},'runs');
  end
  if (isfield(obj,'sse'))
    data = [obj.sse];
  elseif (isfield(obj,'TestSSE'))
    data = [obj.TestSSE];
  elseif (isfield(obj,'like'))
    data = [obj.like];
  elseif (isfield(obj,'TestLike'))
    data = [obj.TestLike];
  elseif (isfield(obj,'bic'))
    data = [obj.bic];
  end
  K = 1:size(data,1);   % set K to the size of the first data field
end
if (isempty(indK))
  indK = 1:length(K);  % set indK to the max if it wasn't specified
end
K = K(indK);  % choose the K that we want to plot

% test for presence of labels
if (isfield(ops,'labels'))
  labels = ops.labels;
else
  labels = methods;
end

if (isfield(ops,'OrderRange') & prod(size(ops.OrderRange))~=1)
  error('show_full: only works with a scalar OrderRange\n');
end


% build the test score matrices
sse_mn=[];  sse_std=[];  like_mn=[];  like_std=[];  bic_mn=[];  bic_std=[];
for i=1:NumMethods
  obj = getfield(mccv,methods{i});
  if (isfield(obj,'runs'))
    obj = getfield(mccv,methods{i},'runs');
  end

  % SSE
  if (isfield(obj,'sse') | isfield(obj,'TestSSE'))
    if (isfield(obj,'sse')), data = [obj.sse];  else data = [obj.TestSSE]; end
    if (size(data,1)~=length(indK)),  data = data(indK,:); end
    sse_mn = [sse_mn mean(data,2)];
    sse_std = [sse_std std(data,0,2)];
  end

  % TestLike
  if (isfield(obj,'like') | isfield(obj,'TestLike'))
    if (isfield(obj,'like')), data = [obj.like]; else data = [obj.TestLike]; end
    if (size(data,1)~=length(indK)),  data = data(indK,:); end
    like_mn = [like_mn mean(data,2)];
    like_std = [like_std std(data,0,2)];
  end

  % BIC
  if (isfield(obj,'bic'))
    data = [obj.bic];
    if (size(data,1)~=length(indK)),  data = data(indK,:,:); end
    if (ndims(data)==3)
      bic_mn = [bic_mn mean(data,3)];
      bic_std = [bic_std std(data,0,3)];
    else
      bic_mn = [bic_mn mean(data,2)];
      bic_std = [bic_std std(data,0,2)];
    end
  end
end

fpos = [0.0821    0.5210    0.4000    0.4000];
apos = [0.1768    0.1500    0.7125    0.7690];
lpos = [0.6788    0.7727    0.1929    0.1126];


figHnds=[];
% plot SSE
if (~isempty(sse_mn) & ops.PlotSSE)
  figHnds(end+1) = figure('units','norm','posi',fpos);
  methodHnd = [];
  Q = size(sse_mn,2);
  if (isempty(ops)), K = (1:size(sse_mn,1))'; end
  for q=1:Q
    spec = colorplot(q,ops);
    h = errorbar(K,sse_mn(:,q),sse_std(:,q),spec);
    set(h,'Tag',labels{q});
    set(h,'LineWidth',ops.LineWidth);
    set(h,'MarkerSize',ops.MarkerSize);
    methodHnd(end+1) = h(2);  % get the lines, not the error bars
    set(h(1),'Tag','ErrorBar');
    hold on;
  end
  set(gca,'posi',apos);
  set(gca,'xtick',K);
  set(gca,'XLim',[K(1)-1 K(end)+1]);
  set(gca,'FontSize',ops.FontSize);
  ylabel('Test SSE');
  xlabel('Number of components');
  h = legend(methodHnd,labels);
  set(h,'posi',lpos);
  h = get(h,'Children');
  set(h(end),'Interpreter','none');
  legend boxoff;
  set(gcf, 'PaperPositionMode', 'auto');
  set(gcf, 'Renderer', 'painters');
end


% plot log-likelihood
if (~isempty(like_mn) & ops.PlotLike)
  figHnds(end+1) = figure('units','norm','posi',fpos);
  methodHnd = [];
  Q = size(like_mn,2);
  if (isempty(ops)), K = (1:size(like_mn,1))'; end
  labels = labels(1:Q);
  for q=1:Q
    spec = colorplot(q,ops);
    h = errorbar(K,like_mn(:,q),like_std(:,q),spec);
    set(h,'Tag',labels{q});
    set(h,'LineWidth',ops.LineWidth);
    set(h,'MarkerSize',ops.MarkerSize);
    methodHnd(end+1) = h(2);  % get the lines, not the error bars
    set(h(1),'Tag','ErrorBar');
    hold on;
  end
  set(gca,'posi',apos);
  set(gca,'xtick',K);
  set(gca,'XLim',[K(1)-1 K(end)+1]);
  set(gca,'FontSize',ops.FontSize);
  ylabel('Test Log-likelihood');
  xlabel('Number of components');
  h = legend(methodHnd,labels);
  set(h,'posi',lpos);
  h = get(h,'Children');
  set(h(end),'Interpreter','none');
  legend boxoff;
  set(gcf, 'PaperPositionMode', 'auto');
  set(gcf, 'Renderer', 'painters');
end


% plot BIC
if (~isempty(bic_mn) & ops.PlotBIC)
  figHnds(end+1) = figure('units','norm','posi',fpos);
  methodHnd = [];
  Q = size(bic_mn,2);
  if (isempty(ops)), K = (1:size(bic_mn,1))'; end
  labels = labels(NumMethods-Q+1:end);
  for q=1:Q
    spec = colorplot(q,ops);
    h = errorbar(K,bic_mn(:,q),bic_std(:,q),spec);
    set(h,'Tag',labels{q});
    set(h,'LineWidth',ops.LineWidth);
    set(h,'MarkerSize',ops.MarkerSize);
    methodHnd(end+1) = h(2);  % get the lines, not the error bars
    set(h(1),'Tag','ErrorBar');
    hold on;
  end
  set(gca,'posi',apos);
  set(gca,'xtick',K);
  set(gca,'XLim',[K(1)-1 K(end)+1]);
  set(gca,'FontSize',ops.FontSize);
  ylabel('BIC');
  xlabel('Number of components');
  h = legend(methodHnd,labels);
  set(h,'posi',lpos);
  h = get(h,'Children');
  set(h(end),'Interpreter','none');
  legend boxoff;
  set(gcf, 'PaperPositionMode', 'auto');
  set(gcf, 'Renderer', 'painters');
end


if (nargout>0)
  table.sse_mn  = sse_mn;
  table.sse_std = sse_std;
  
  table.like_mn  = like_mn;
  table.like_std = like_std;
  
  table.bic_mn  = bic_mn;
  table.bic_std = bic_std;
end