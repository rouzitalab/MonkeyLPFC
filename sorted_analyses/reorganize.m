s = wf.trialInfo.trialStartTimes;
e = wf.trialInfo.trialEndTimes;
nTrials=length(s);
units = sua;
labels=conditions;
clc
sps=1000;
% idx = [];
tr = [];
% count=0;
for i = 1:nTrials
    if e(i)<(s(i)+2.2*sps)
%         count = count+1;
        tr = [tr, i];
%         idx=[idx, ceil(s(i)-s(1)):ceil(e(i)-s(1))];
%         s(i+1)=s(i);
%         e(i+1)=e(i+1)-e(i)+s(i)-1;
    end
end
% units(:,idx)=[];
% labels(tr)=[];
% nTrials = nTrials - count;
% e(tr)=[];
% s(tr)=[];

X = zeros(nTrials, ceil(1.5*sps),nUnits);
for i = 1:nTrials
    fprintf("Trial #%i\n",i);
    id_s = ceil(s(i)+0.7*sps-s(1));
    id_e = ceil(s(i)+2.2*sps-s(1));
    X(i,:,:) = transpose(units(:,id_s:id_e-1));
end

X(tr,:,:)=[];
labels(tr)=[];

%% Shuffle Data
idx = randperm(length(labels));
X = X(idx,:,:);
labels = labels(idx);

X = reshape(X, size(X,1), size(X,2)*size(X,3));
% xtr=reshape(X(500:end,:,:),size(X(500:end,:,:),1), size(X(500:end,:,:),2)*size(X(500:end,:,:),3));
% xts=reshape(X(1:500,:,:),size(X(1:500,:,:),1), size(X(1:500,:,:),2)*size(X(1:500,:,:),3));
% ytr=labels(500:end);
% yts=labels(1:500);


% mdl = fitcecoc(X, labels);
% cvmdl = crossval(mdl);
% generror = kfoldLoss(cvmdl);



