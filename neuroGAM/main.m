clear;
clc;
path = 'G:\Projects\neuroGAM\data';
data_list = dir(fullfile(path, '*.mat'));
data_path = append(path, '\',data_list(1).name);
data1 = load(data_path);
pos = data1.gaze;
v = diff(pos,1,2);
pos = pos(:,1:end-1,:);
speed = sqrt(v(:,:,1).^2 + v(:,:,2).^2);
target = data1.target;
saccade = data1.saccade;
spikes = data1.spikes;
spikes(isnan(spikes))=0;
tmp = data1.outcome;
flag = find(tmp~=-1);
tmp = data1.color;
color = double.empty();
color(length(tmp)) = 0;
for i=1:length(tmp)
    if isequal(tmp{i},'r')
        color(i)=1;
    elseif isequal(tmp{i},'g')
        color(i)=2;
    elseif isequal(tmp{i},'b')
        color(i)=3;
    end
end
tp = size(spikes,2);
ch = size(spikes,3);
pos = pos(flag,:,:);
speed = speed(flag,:);
color = color(flag);
target = target(flag);
saccade = saccade(flag);
spikes = spikes(flag,:,:);
tmp = permute(spikes,[2,1,3]);
tmp = downsample(tmp,10);
spikes = permute(tmp, [2,1,3]);
spikes(:,end+1,:)=0;
spikes = abs(spikes);
spikes = uint16(spikes);
for i=1:size(pos,1)
    for j=1:size(pos,2)
        if isnan(pos(i,j,1))
            pos(i,j,1) = pos(i,j-1,1);
        end
        if isnan(pos(i,j,2))
            pos(i,j,2) = pos(i,j-1,2);
        end
        if isnan(speed(i,j))
            speed(i,j) = speed(i,j-1);
        end
    end
end
%%%
tmp = sum(spikes,3);
yt = reshape(tmp, [size(tmp,1)*size(tmp,2), 1]);
tmp = reshape(pos, [size(pos,1)*size(pos,2), size(pos,3)]);
tmp2 = reshape(speed, [size(speed,1)*size(speed,2), 1]);
tmp3 = [];
tmp4 = [];
for i=1:size(spikes,1)
tmp3((i-1)*size(spikes,2)+1:i*size(spikes,2))=saccade(i);
tmp4((i-1)*size(spikes,2)+1:i*size(spikes,2))=color(i);
end
tmp3 = tmp3.';
tmp4 = tmp4.';
xt = {tmp; tmp2};
% xt = xt.';
prs.varname = {'Position', 'Speed'};
prs.vartype = {'2D', '1D'};
prs.basistype = {'boxcar', 'boxcar'};
prs.nbins = {[10,10],10};
prs.binrange = [];
prs.nfolds = 10;
prs.dt = 0.01;
prs.filtwidth = 10;
prs.linkfunc = 'log';
prs.lambda = {10, 10};
prs.alpha = 0.05;
prs.varchoose = [0, 0];
prs.method = 'Forward';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[models, x, nvars] = BuildGAM(xt,yt,prs);
% PlotGAM(models,prs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nvars = length(xt);
% 
% %% load analysis parameters
% prs = struct2cell(prs);
% [xname,xtype,basistype,nbins,binrange,nfolds,dt,filtwidth,linkfunc,lambda,alpha,varchoose,method] = deal(prs{:});
% method = lower(method);
% 
% %% define undefined analysis parameters
% if isempty(alpha), alpha = 0.05; end
% if isempty(lambda), lambda = cell(1,nvars); lambda(:) = {5e1}; end
% if isempty(linkfunc), linkfunc = 'log'; end
% if isempty(filtwidth), filtwidth = 3; end
% if isempty(nfolds), nfolds = 10; end
% if isempty(nbins)
%     nbins = cell(1,nvars); nbins(:) = {10}; % default: 10 bins
%     nbins(strcmp(xtype,'2D')) = {[10,10]};
% end
% if isempty(binrange), binrange = []; end
% if isempty(dt), dt = 1; end
% % define bin range
% if isempty(binrange)
%     binrange = mat2cell([min(cell2mat(xt.'));max(cell2mat(xt.'))],2,strcmp(xtype,'2D')+1);
%     binrange(strcmp(xtype,'event')) = {[-0.36;0.36]}; % default: -360ms to 360ms temporal kernel
% end
% % express bin range in units of dt for temporal kernels
% indx = find(strcmp(xtype,'event'));
% for i=indx, binrange{i} = round(binrange{i}/dt); end
% 
% %% compute inverse-link function
% if strcmp(linkfunc,'log')
%     invlinkfunc = @(x) exp(x);
% elseif strcmp(linkfunc,'identity')
%     invlinkfunc = @(x) x;
% elseif strcmp(linkfunc,'logit')
%     invlinkfunc = @(x) exp(x)./(1 + exp(x));
% end
% fprintf(['...... Fitting ' linkfunc '-link model\n']);
% 
% %% encode variables in 1-hot format
% x = cell(1,nvars); % 1-hot representation of xt
% basis = cell(1,nvars); % bin centres
% nprs = cell(1,nvars); % number of parameters (weights)
% for i=1:nvars
% %     [x{i},xc{i},nprs{i}] = Encode1hot(xt{i}, xtype{i}, binrange{i}, nbins{i});
%     [x{i},basis{i},nprs{i}] = RecodeInput(xt{i}, xtype{i}, binrange{i}, nbins{i}, basistype{i}, dt);
%     Px{i} = sum(x{i}~=0)/size(x{i},1); 
% end
% 
% %% define filter to smooth the firing rate
% t = linspace(-2*filtwidth,2*filtwidth,4*filtwidth + 1);
% h = exp(-t.^2/(2*filtwidth^2));
% h = h/sum(h);
% 
% switch method
%     case {'forward','backward'}
%         %% define model combinations to fit
%         Model = DefineModels(nvars,1:nvars,varchoose);
% %         Model = Model(end-nvars:end); % fit only nvars+1 models (full model & models missing one variable)
%         %% fit all models
% %         models = FitModels(Model,x,xtype,nprs,yt,dt,h,nfolds,lambda,linkfunc,invlinkfunc);        
% %         %% select best model
% %         fprintf('...... Performing model selection\n');
% %         testFit = cell2mat(models.testFit); nrows = size(testFit,1);
% %         LLvals = reshape(testFit(:,4),nfolds,nrows/nfolds); % 4th column contains likelihood values
% %         if strcmp(method,'forward'), models.bestmodel = ForwardSelect(Model,LLvals,alpha);
% %         elseif strcmp(method,'backward'), models.bestmodel = BackwardEliminate(Model,LLvals,alpha); end
%     case {'fastforward'}
%         Model = DefineModels(nvars,1:nvars,varchoose);
%         models.class = Model; 
%         for n = 1:length(Model), models.testFit{n,1} = nan(nfolds,7); models.trainFit{n,1} = nan(nfolds,7); models.wts{n,1} = nan(1,sum(cell2mat(nprs).*models.class{n})); end
%         models = FastForwardSelect(Model,models,x,xtype,nprs,yt,dt,h,nfolds,lambda,linkfunc,invlinkfunc,alpha);
%     case {'fastbackward'}
%         Model = DefineModels(nvars,1:nvars,varchoose);
%         models.class = Model; 
%         for n = 1:length(Model), models.testFit{n,1} = nan(nfolds,7); models.trainFit{n,1} = nan(nfolds,7); models.wts{n,1} = nan(1,sum(cell2mat(nprs).*models.class{n})); end
%         models = FastBackwardEliminate(Model,models,x,xtype,nprs,yt,dt,h,nfolds,lambda,linkfunc,invlinkfunc,alpha);
% end
% 
% nModels = length(Model);
% models.class = Model; models.testFit = cell(nModels,1); models.trainFit = cell(nModels,1); models.wts = cell(nModels,1);
% n = 1;
% X = cell2mat(x(Model{n})); % X{i} stores inputs for the i^th model
% Xtype = xtype(Model{n});
% Nprs = cell2mat(nprs(Model{n}));
% Lambda = lambda(Model{n});
% fprintf('\t- Fitting model %d of %d\n', bin2dec(num2str(Model{n})), bin2dec(num2str(ones(1,length(Model{n})))));
% nchunks = 5; % divide data into these many chunks (5 is large enough for experiments that last <3 hrs)
% [~,nprs] = size(X); % number of parameters to fit
% nsections = nfolds*nchunks;
% 
% % divide the data up into nfolds*nchunks pieces
% y = yt;
% edges = round(linspace(1,numel(y)+1,nsections+1));
% 
% % initialize outputs
% testFit = nan(nfolds,7); % to hold 7 values: var ex, var ex pseudo, correlation, llh increase, mse, # of spikes, length of test data
% trainFit = nan(nfolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of train data
% paramMat = nan(nfolds,nprs);
% response.true_test = [];
% response.pred_test = [];
% k=1;
% test_ind = cell2mat(arrayfun(@(j) edges((j-1)*nfolds + k):edges((j-1)*nfolds + k + 1)-1, 1:nchunks,'UniformOutput',false));
%     
% test_spikes = y(test_ind); %test spiking
% smooth_spikes_test = convn(test_spikes,h,'same'); %returns vector same size as original
% smooth_fr_test = smooth_spikes_test./dt;
% test_X = X(test_ind,:);
% 
% % get training data
% train_ind = setdiff(1:numel(y),test_ind);
% train_spikes = y(train_ind);
% smooth_spikes_train = convn(train_spikes,h,'same'); %returns vector same size as original
% smooth_fr_train = smooth_spikes_train./dt;
% train_X = X(train_ind,:);    
% data{1} = train_X; data{2} = train_spikes;
% opts = optimset('Gradobj','on','Hessian','on','Display','off','Algorithm','trust-region');
% init_param = 1e-3*randn(nprs, 1);
% param = fminunc(@(param) log_link(param,data,Xtype,Nprs,Lambda),init_param,opts);
% fr_hat_test = invlinkfunc(test_X * param)/dt;
% smooth_fr_hat_test = convn(fr_hat_test,h,'same'); %returns vector same size as original
% 
% % variance explained
% sse = sum((smooth_fr_hat_test-smooth_fr_test).^2);
% sst = sum((smooth_fr_test-mean(smooth_fr_test)).^2);
% varExplain_test = 1-(sse/sst);
% 
% % linear correlation
% correlation_test = corr(smooth_fr_test,smooth_fr_hat_test,'type','Pearson');
% 
% % log-likelihood increase from "mean firing rate model" - NO SMOOTHING
% r_test = invlinkfunc(test_X * param); n = test_spikes; meanFR_test = nanmean(test_spikes);
% log_llh_test_model = nansum(r_test-n.*log(r_test)+log(factorial(n)))/sum(n);
% log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
% log_llh_test_best = nansum(n-n.*log(n)+log(gamma(n+1)))/sum(n);
% varExplainPseudo_test = (log_llh_test_model - log_llh_test_mean)/(log_llh_test_best - log_llh_test_mean);
% log_llh_test = (-log_llh_test_model + log_llh_test_mean); % nats/spike
% log_llh_test = log_llh_test/log(2); % convert to bits/spike
% 
% % mean-squared-error
% mse_test = nanmean((smooth_fr_hat_test-smooth_fr_test).^2);
% testFit(k,:) = [varExplain_test varExplainPseudo_test correlation_test log_llh_test mse_test sum(n) numel(test_ind)];







function [J,G,H] = roughness_penalty(param,vartype,lambda)
if strcmp(vartype,'2D')
    numParam = numel(param);
    D1 = spdiags(ones(sqrt(numParam),1)*[-1 1],0:1,sqrt(numParam)-1,sqrt(numParam));
    DD1 = D1'*D1;
    M1 = kron(eye(sqrt(numParam)),DD1); M2 = kron(DD1,eye(sqrt(numParam)));
    M = (M1 + M2);
    % compute J, G, and H
    J = lambda*0.5*param'*M*param;
    G = lambda*M*param;
    H = lambda*M;
elseif strcmp(vartype,'1Dcirc')
    numParam = numel(param);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    % to correct the smoothing across first and last bin
    DD1(1,:) = circshift(DD1(2,:),[0 -1]);
    DD1(end,:) = circshift(DD1(end-1,:),[0 1]);
    % compute J, G, and H
    J = lambda*0.5*param'*DD1*param;
    G = lambda*DD1*param;
    H = lambda*DD1;
elseif any(strcmp(vartype,{'1D','event'}))
    numParam = numel(param);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    % compute J, G, and H
    J = lambda*0.5*param'*DD1*param;
    G = lambda*DD1*param;
    H = lambda*DD1;
elseif strcmp(vartype,'0D')
    J = lambda*sum(abs(param));
    G = lambda*sign(param);
    H = diag(zeros(1,length(param)));
end
end