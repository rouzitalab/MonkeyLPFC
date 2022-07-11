% TrajectoryClustering_example2.m
% First, create some data from a mixture of two bivariate Gaussian distributions using the mvnrnd function:

% clear all;close all;

sigma = 10; %noise variance
x=[0:.1:3.5];
y1=x.^2;
y2=3-x.^2;
y3=x.^3;
y4=2*x;

ops.K=4;
ops.zero='none';
ops.order=3;
ops.method = 'gmix'; %lrm, kmeans , gmix - ListModels();
ops.IterLimit=200;
ops.NumEMStarts=4;
ops.ShowGraphics	=0;
  
origidx=zeros(40,1);
origGroup{1}=1*ones(10,1)';
origGroup{2}=2*ones(10,1)';
origGroup{3}=3*ones(10,1)';
origGroup{4}=4*ones(10,1)';

Trajs.X=x;
for i = 1:40
    origGroupNum = ceil(i/10);
    origidx(i,1)=origGroupNum;
    switch origGroupNum
        case 1
            yi=y1+sigma*randn(1,length(x));
        case 2
            yi=y2+sigma*randn(1,length(x));
        case 3
            yi=y3+sigma*randn(1,length(x));
        case 4
            yi=y4+sigma*randn(1,length(x));
    end
    % Trajs.Y{i}=[yi;x]';
    Trajs.Y{i}=[x;yi]';
end




model = curve_clust(Trajs,ops);
showmodel(model,Trajs);
AIC=2*model.K - 2*model.TrainLhood_ppt;
idx=model.C;


%%
% x = meas;
unorderedGroupEst{1} = origidx(idx == 1)';
unorderedGroupEst{2} = origidx(idx == 2)';
unorderedGroupEst{3} = origidx(idx == 3)';
unorderedGroupEst{4} = origidx(idx == 4)';

[orderedGroupEst,map, cm] = as_classify(origGroup,unorderedGroupEst);

H = as_infoClust1(cm);
Pc = sum(diag(cm)) / sum(sum(cm));