function dat1=as_convertres2dat2(filename,analysisBlock,excludeNoiseUnits,analysisGroup)
% analysisBlock may be a scalar or vector
% gpfa requires matrix dimension agreement.  Can take max length or min
% length
% corrects a problem in as_convertres2dat1 in which spike pattern
% consisting of zeros was used (instead of spike time).  This created empty
% spike trains when the spike pattern toggle in analysis parameters is
% turned off
if nargin<4
    analysisGroup=[];
end

unitDefined =0 ;  % for now all units
analysisBase=[1:8];
analysisPeriods = [];
jPeriod=1;
load(filename, '-mat');
while isempty(res.aPeriod(jPeriod).unitExcluded)
    jPeriod=jPeriod+1;
end
unitsExcluded=res.aPeriod(jPeriod).unitExcluded; % Since the units are the same for every analysis period
                                            % Therefore arbitrarily chose 1
                                            % or 2 or 3

nTrodes =32;
if excludeNoiseUnits
    unitsExcluded(:,1)=-2;
end
unitIndices = find(unitsExcluded==0);
unitSubs    = ind2sub(size(unitsExcluded),unitIndices);
totNumUnits = length(find(unitsExcluded==0));
[i_UnitIndex j_UnitIndex]    = ind2sub(size(unitsExcluded),unitIndices);
for i_block = 1:length(analysisBlock)
    analysisPeriods = [analysisPeriods (analysisBlock(i_block)-1 )*8+analysisBase];
end

NewRes = makeResMuchSmaller( res,analysisPeriods, unitsExcluded );


% memory bullshit
if 1
clear res
res=NewRes;
clear NewRes
save Temp
clear all
load Temp
end

trialID = [];
durations = [];  % NOte HACK= ideally all durations are same length
for i_Period = 1: length(analysisPeriods)
    trialsInaPeriod = res.aPeriod(analysisPeriods(i_Period)).trial ;
    aPer            = analysisPeriods(i_Period)* ones(size(trialsInaPeriod));
    aBlock          = ceil(aPer./8);
    trialIndex      = 1:length(aPer);
    trialID         = [trialID; [trialsInaPeriod aPer res.aPeriod(analysisPeriods(i_Period)).duration' trialIndex' aBlock]];
    durations = [durations res.aPeriod(analysisPeriods(i_Period)).duration];
end
duration = max(durations);  %%% CONSIDER USING MIN
       
dat1Length = size(trialID,1);

for itrial = 1: dat1Length
    
    dat1(itrial).trialId = trialID(itrial, 1);
    aPeriodforthistrial = trialID(itrial,2);
    % this is when you use actual trial length, but we need idententical trial length -->
    % spikes =  zeros(totNumUnits,trialID(itrial,3));
    spikes =  zeros(totNumUnits,duration);
    for u = 1:totNumUnits
        % spikes(iUnit,res.aPeriod(4).trode(30).unit(1).spikeTimes(4))=1;
        temp = res.aPeriod(aPeriodforthistrial).trode(i_UnitIndex(u)).unit(j_UnitIndex(u)).spikeTimes{trialID(itrial,4)}';
        % spikes(u,:) = temp(1,1:trialID(itrial,3));  % florian used max size of all trials to allow spike density funtion, but we dont have this restriction here
        spikes(u,1:trialID(itrial,3)) = logical(zeros(1,trialID(itrial,3)));
        spikes(u,temp)                = 1;
        clear temp
    end
    dat1(itrial).spikes = spikes;
    clear spikes
    dat1(itrial).analysisPeriod = trialID(itrial,2);
    if ~ isempty(analysisGroup)
        
        if ismember(dat1(itrial).analysisPeriod,analysisGroup{1}.aPeriods)
            dat1(itrial).analysisPeriod=1;
        elseif ismember(dat1(itrial).analysisPeriod,analysisGroup{2}.aPeriods)
            dat1(itrial).analysisPeriod=2;
        else 
            dat1(itrial).analysisPeriod=-1; % HACK this is a flag;
            
        end
    end
    
    dat1(itrial).analysisBlock  = trialID(itrial,5);
    dat1(itrial).dir            = trialID(itrial,2) - 8* (trialID(itrial,5) - 1);
end

j=1;
for i =1:dat1Length
    if dat1(i).analysisPeriod==-1
    else
        dat2(j) = dat1(i);
        j=j+1;
    end
end
try
dat1=dat2;
catch
    error('These analysis periods do not exist in your dataset');
end
dummy=3;

% process allowing only fixed intervals
k=1;

% for now just does trial intervals for 2 analysis groups
for i = 1:size(dat1,2)
    if dat1(i).analysisPeriod == 1 
        if isempty(analysisGroup{1}.trialIntervals)
            dat3(k)=dat1(i);
            k=k+1;
        else
        for j = 1:size(analysisGroup{1}.trialIntervals,1)
            if any(dat1(i).trialId== [analysisGroup{1}.trialIntervals(j,1):analysisGroup{1}.trialIntervals(j,2)])
                dat3(k)=dat1(i);
                k=k+1;
            end
        end
        end
    end
    
    if dat1(i).analysisPeriod == 2
        if isempty(analysisGroup{2}.trialIntervals)
            dat3(k)=dat1(i);
            k=k+1;
        else
        for j = 1:size(analysisGroup{2}.trialIntervals,1)
            if any(dat1(i).trialId== [analysisGroup{2}.trialIntervals(j,1):analysisGroup{2}.trialIntervals(j,2)])
                dat3(k)=dat1(i);
                k=k+1;
            end
        end
        end
    end
    
end
dat1=dat3;