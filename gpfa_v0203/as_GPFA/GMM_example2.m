% GMM_example2

% First, create some data from a mixture of two bivariate Gaussian distributions using the mvnrnd function:
nSamp=10;
MU1 = [-2 2];
SIGMA1 = [1 0; 0 1];
MU2 = [3 0];
SIGMA2 = [1 0; 0 1];
MU3 = [-3 -3];
SIGMA3 = [1 0; 0 1];
x1=mvnrnd(MU1,SIGMA1,nSamp);
x2=mvnrnd(MU2,SIGMA2,nSamp);
x3=mvnrnd(MU3,SIGMA3,nSamp);
X = [x1;x2;x3];

origidx = [ones(nSamp,1); 2*ones(nSamp,1); 3*ones(nSamp,1)];
origidx = shuffle(origidx);
origGroup{1}=ones(nSamp,1)';
origGroup{2}=2*ones(nSamp,1)';
origGroup{3}=3*ones(nSamp,1)';

scatter(X(:,1),X(:,2),10,'.')
%%

%%%%

AIC = zeros(1,4);
obj = cell(1,4);
for k = 1:4
    obj{k} = gmdistribution.fit(X,k);
    AIC(k)= obj{k}.AIC;
end

[minAIC,numComponents] = min(AIC);
numComponents
%%%%%
%%

% Next, fit a two-component Gaussian mixture model:

options = statset('Display','final');
obj = gmdistribution.fit(X,3,'Options',options);% was 2
hold on
h = ezcontour(@(x,y)pdf(obj,[x y]),[-8 6],[-8 6]);
hold off

%%
%%

idx = cluster(obj,X);
cluster1 = X(idx == 1,:);
cluster2 = X(idx == 2,:);
cluster3 = X(idx == 3,:);
unorderedGroupEst{1} = origidx(idx == 1)';
unorderedGroupEst{2} = origidx(idx == 2)';
unorderedGroupEst{3} = origidx(idx == 3)';
figure,
hold on,
h1 = scatter(cluster1(:,1),cluster1(:,2),10,'r.');
h2 = scatter(cluster2(:,1),cluster2(:,2),10,'g.');
h3 = scatter(cluster3(:,1),cluster3(:,2),10,'b.');
AXIS([-8 6 -8 6]);
hold off
% legend([h1 h2],'Cluster 1','Cluster 2', 'Location','NW');
% legend([h1 h2],'Cluster 1','Cluster 2');

%%
% x = meas;
[orderedGroupEst,map, cm] = as_classify(origGroup,unorderedGroupEst);

H = as_infoClust1(cm);
Pc = sum(diag(cm)) / sum(sum(cm));
