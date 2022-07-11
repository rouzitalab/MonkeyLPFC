function dat1=as_convertres2dat1(res,analysisBlock,excludeNoiseUnits)
% analysisBlock may be a scalar or vector
% gpfa requires matrix dimension agreement.  Can take max length or min
% length

unitDefined =0 ;  % for now all units
analysisBase=[1:8];
analysisPeriods = [];
jPeriod=1;
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
        temp = res.aPeriod(aPeriodforthistrial).trode(i_UnitIndex(u)).unit(j_UnitIndex(u)).pattern(trialID(itrial,4),:);
        % spikes(u,:) = temp(1,1:trialID(itrial,3));  % florian used max size of all trials to allow spike density funtion, but we dont have this restriction here
        spikes(u,1:trialID(itrial,3)) = temp(1,1:trialID(itrial,3));
        ss = sum(temp(1,1:trialID(itrial,3)));
        if ss~=0
            [itrial u ss]
        end
    end
    dat1(itrial).spikes = spikes;
    clear spikes
    dat1(itrial).analysisPeriod = trialID(itrial,2);
    dat1(itrial).analysisBlock  = trialID(itrial,5);
    dat1(itrial).dir            = trialID(itrial,2) - 8* (trialID(itrial,5) - 1);
end

dummy=3;

    
