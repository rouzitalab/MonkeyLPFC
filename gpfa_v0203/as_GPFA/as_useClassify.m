% as_useClassify.m
% July 2012
% use on synthetic data to classify based on number of spikes in time
% interval
% revised Oct 2012 to allow folding with bootstrap test data adn train 
% data 
% use upper section for dat1 (spike data)

if 1
N=size(dat1,2);
numIter = 2*N;
NumTest = round(N/10);
PC = [];

for i = 1:numIter
    z=randperm(N);
    ind_test=z(1:NumTest);
    ind_train = z;
    ind_train(ind_test)=[];
    
    
    
    % ind_train = 1:60;
    % ind_test  = 61:80;
    SPEC=[];
    MEAS=[];
    for i=1:size(dat1,2)
        dat1(i).totSpikes=sum(dat1(i).spikes')'; % make spike counts from spike trains
        SPEC=[SPEC;dat1(i).dir];
        MEAS=[MEAS;dat1(i).totSpikes'];
    end
    N = size(MEAS,1);
    % below may classify with 'linear', 'quadratic',
    % 'diagLinear', 'diagQuadratic', or 'mahalanobis'., and may add 'prior
    % argument
    ldaClass = classify(MEAS(ind_test,1:size(MEAS,2)),MEAS(ind_train,1:size(MEAS,2)),SPEC(ind_train),'linear');% 'mahalanobis'
    bad = ~strcmp(ldaClass,SPEC(ind_train));
    ldaResubErr = sum(bad) / N; % compute the resubstitution error, which is the misclassification error
    % (the proportion of misclassified observations) on the training set
    [ldaResubCM,grpOrder] = confusionmat(SPEC(ind_test),ldaClass); %ldaResubCM is the confusion matrix
    Pc=sum(diag(ldaResubCM))/sum(sum(ldaResubCM)); % Pc is probability correct
    PC= [PC;Pc];
end
end
%% use this for FA data
shuffleIdx=0;
dat1DIR=[];
for i=1:size(dat1,2)
    dat1DIR(i,:)=[dat1(i).trialId dat1(i).dir];
end

SID=[];
SX=[];
for i=1:size(seqTrain,2)
    SID(i)=dat1DIR(seqTrain(i).trialId,2);
    % SX(i,:)=mean(seqTrain(i).xorth')';
      SX(i,:,:)=seqTrain(i).xorth;
end
if shuffleIdx
    SID=SID(randperm(size(SID,2)));
end

N=size(SID,2);
numIter = 2*N;
NumTest = 1; %round(N/20);
PC = [];

for i = 1:numIter
    z=randperm(N);
    ind_test=z(1:NumTest);
    ind_train = z;
    ind_train(ind_test)=[];
    
    
    
    % ind_train = 1:60;
    % ind_test  = 61:80;
    
    
    
    SPEC=SID;
    MEAS=SX;
    j=1;
    while any(hist(SPEC(ind_train),8)==1)
       
        z=randperm(N);
        ind_test=z(1:NumTest);
        ind_train = z;
        ind_train(ind_test)=[];
        j=j+1;
        if j==100
            error('Need to have at last two points per category in training');
        end
    end
    
    
    
    % ind_train = 1:60;
    % ind_test  = 61:80;
    
    
    
    SPEC=SID;
    MEAS=SX;
    % MEAS(:,[5:8])=[];    
        
        
        
        
    N = size(MEAS,1);
    % below may classify with 'linear', 'quadratic',
    % 'diagLinear', 'diagQuadratic', or 'mahalanobis'., and may add 'prior
    % argument
    ldaClass = classify(MEAS(ind_test,1:size(MEAS,2)),MEAS(ind_train,1:size(MEAS,2)),SPEC(ind_train),'diagQuadratic');% 'mahalanobis'
    bad = ~strcmp(ldaClass,SPEC(ind_train));
    ldaResubErr = sum(bad) / N; % compute the resubstitution error, which is the misclassification error
    % (the proportion of misclassified observations) on the training set
    [ldaResubCM,grpOrder] = confusionmat(SPEC(ind_test),ldaClass); %ldaResubCM is the confusion matrix
    Pc=sum(diag(ldaResubCM))/sum(sum(ldaResubCM)); % Pc is probability correct
    PC= [PC;Pc];
end
mean(PC)
