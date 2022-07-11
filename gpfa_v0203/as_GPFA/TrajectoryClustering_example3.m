% TrajectoryClustering_example3.m
% First, create some data from a mixture of two bivariate Gaussian distributions using the mvnrnd function:

% clear all;close all;
% example3 uses data from simNeuron fed to GPFA

load sim1

ops.K=8;
ops.zero='none';
ops.order=3;
ops.method = 'gmix';
ops.IterLimit=200;
ops.NumEMStarts=8;
ops.ShowGraphics	=1;
  
origidx=[];
origGroup{1}=1*ones(10,1)';
origGroup{2}=2*ones(10,1)';
origGroup{3}=3*ones(10,1)';
origGroup{4}=4*ones(10,1)';
origGroup{5}=5*ones(10,1)';
origGroup{6}=6*ones(10,1)';
origGroup{7}=7*ones(10,1)';
origGroup{8}=8*ones(10,1)';

% Trajs.X=x;
for i = 1:8
    origidx=[origidx;origGroup{i}];
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
unorderedGroupEst{5} = origidx(idx == 5)';
unorderedGroupEst{6} = origidx(idx == 6)';
unorderedGroupEst{7} = origidx(idx == 7)';
unorderedGroupEst{8} = origidx(idx == 8)';

[orderedGroupEst,map, cm] = as_classify(origGroup,unorderedGroupEst);

H = as_infoClust1(cm);
Pc = sum(diag(cm)) / sum(sum(cm));