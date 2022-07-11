function dat1=as_convertres2dat1a(res,analysisBlock,errorBlock, excludeNoiseUnits)
% a - makes analysisBlock a matrix, and incudes errorblock 0 for not error,
% 1 for error
% analysisBlock may be a scalar or vector
% gpfa requires matrix dimension agreement.  Can take max length or min
% length

unitDefined =0 ;  % for now all units
analysisBase=[1:8];
analysisPeriods = [];
unitsExcluded=res.aPeriod(1).unitExcluded; % Since the units are the same for every analysis period
                                            % Therefore arbitrarily chose 1
nAB = length(analysisBlock);
if nAB ~= length(errorBlock)
    error('analysisBlock must have same number of entries as errorBlock');
end


nTrodes =32;
if excludeNoiseUnits
    unitsExcluded(:,1)=-2;
end
unitIndices = find(unitsExcluded==0);
unitSubs    = ind2sub(size(unitsExcluded),unitIndices);
totNumUnits = length(find(unitsExcluded==0));
[i_UnitIndex j_UnitIndex]    = ind2sub(size(unitsExcluded),unitIndices);

trialIndex = 1;
for iAB = 1:nAB
    for i_block = 1:length(analysisBlock(iAB))
        analysisPeriods = [analysisPeriods (analysisBlock(iAB)(i_block)-1 )*8+analysisBase];
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
    
    for itrial = 1: dat1Length % use itrial to determine the number of loops, but not as the index!
        dat1(trialIndex).trialId = trialID(trialIndex, 1);
        aPeriodforthistrial = trialIndex(itrial,2);
        % this is when you use actual trial length, but we need idententical trial length -->
        % spikes =  zeros(totNumUnits,trialID(itrial,3));
        spikes =  zeros(totNumUnits,duration);
        for u = 1:totNumUnits
            % spikes(iUnit,res.aPeriod(4).trode(30).unit(1).spikeTimes(4))=1;
            temp = res.aPeriod(aPeriodforthistrial).trode(i_UnitIndex(u)).unit(j_UnitIndex(u)).pattern(trialID(trialIndex,4),:);
            % spikes(u,:) = temp(1,1:trialID(itrial,3));  % florian used max size of all trials to allow spike density funtion, but we dont have this restriction here
            spikes(u,1:trialID(trialIndex,3)) = temp(1,1:trialID(trialIndex,3));
            
        end
        dat1(trialIndex).spikes = spikes;
        clear spikes
        dat1(trialIndex).analysisPeriod = trialID(trialIndex,2);
        dat1(trialIndex).analysisBlock(iAB)  = trialID(trialIndex,5);
        dat1(trialIndex).dir            = trialID(trialIndex,2) - 8* (trialID(trialIndex,5) - 1);
        dat1(trialIndex).iserrortrial = errorBlock(iAB);
        trialIndex=trialIndex+1;
    end
end
dummy=3;

    
