function [ NewRes ] = makeResMuchSmaller( res,analysisPeriods, unitsExcluded )
% written ajs Aug 2011

n_aPar =  size(res.aPeriod,2);
for i_aPar = 1:n_aPar
    if any(analysisPeriods==i_aPar)  %% Note, this still copies over empty analysis periods
        NewRes.aPeriod(i_aPar).trial = res.aPeriod(i_aPar).trial;
        NewRes.aPeriod(i_aPar).trialTimeWindow = res.aPeriod(i_aPar).trialTimeWindow;
        NewRes.aPeriod(i_aPar).fileTimeWindow = res.aPeriod(i_aPar).fileTimeWindow;
        NewRes.aPeriod(i_aPar).duration = res.aPeriod(i_aPar).duration;
        NewRes.aPeriod(i_aPar).unitExcluded = res.aPeriod(i_aPar).unitExcluded;
        for t = 1:size(unitsExcluded,1)
            for u = 1:size(unitsExcluded,2)
                if unitsExcluded(t,u) == 0 && ~isempty (res.aPeriod(i_aPar).trial)
                    NewRes.aPeriod(i_aPar).trode(t).unit(u).spikeTimes = res.aPeriod(i_aPar).trode(t).unit(u).spikeTimes;
                    NewRes.aPeriod(i_aPar).trode(t).unit(u).numSpikes  = res.aPeriod(i_aPar).trode(t).unit(u).numSpikes;
                    NewRes.aPeriod(i_aPar).trode(t).unit(u).meanFR     = res.aPeriod(i_aPar).trode(t).unit(u).meanFR;
                elseif unitsExcluded(t,u) ~= 0 && ~isempty (res.aPeriod(i_aPar).trial)
                    NewRes.aPeriod(i_aPar).trode(t).unit(u).spikeTimes = [];
                    NewRes.aPeriod(i_aPar).trode(t).unit(u).numSpikes  = [];
                    NewRes.aPeriod(i_aPar).trode(t).unit(u).meanFR     = [];
                    
                end
            end
        end
        
    else
        %%% NewRes.aPeriod(i_aPar)=[];
    end
end



