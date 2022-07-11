function saccade = deriveSaccades(gazeData, varargin)
%saccade = deriveSaccades(degData)
% extracts saccades based on algorithm presenented in
% Armstrong et al Neuron 50, 791?798, June 1, 2006
% Velocity check, amplitude check, stat check, refractory check.
% degData is positional data in degrees.
% fs (1000 is default)
% states is a int vector of states matching DEF.EYE_State

%From the paper:
%"Microsaccades coinciding with microstimulation were detected off-line
%using an iterative algorithm based on the intersection of a velocity
%threshold, an amplitude threshold, and statistically significant 
%deflections in the x or y position. The velocity threshold flagged time 
%points at which the instantaneous velocity was above 10 deg/s for a 
%minimum of 10 ms (Bair and O?Keefe, 1998). Two moving windows of 50 ms 
%separated by 25 ms were iterated in *5 ms steps over the x and y 
%components of eye position. At each step, a two sample Kolmogorov-Smirnov
%test (p < 0.01) compared the x and y components of the two periods. 
%If either the x or y component differed significantly, the time point
%at the end of the first window was marked. The amplitude of each eye
%movement was approximated by the displacement of the median x and y 
%components between the two windows. Points with amplitudes >0.1º were 
%flagged (Bair and O?Keefe, 1998). The first of consecutive time points 
%which passed all three criteria was considered a saccade start time. 
%Successive saccades were constrained to start a minimum of 50 ms after 
%any previous saccade."
%
%*They sampled eye position at 200 Hz, so 5 ms steps is once per sample,
%hence the ability to require consecutive samples to meet the criteria.
%
%Adam & Florian made some significant changes.
%1) Their stat check used different windows.
%2) Adam canceled the stat check.
%3) They added an amplitude check on first vs last sample in putative saccade.
%4) They grouped saccades if there were <refrac samples between end_n-1 and
%start_n whereas the paper grouped saccades if their starts were within
%refrac of each other.
%5) They kept the 1 saccade in the group with the largest sum(amp) whereas
%the paper kept the first saccade of the group.
%6) The peak velocity was identified and a window before the peak and after
%the peak were used to calculate the starting point, the ending point, the
%direction, and the amplitude. Only saccades with amplitude within the
%range [7.3 19] were saved.
%
%This last point is a little surprising, because it isn't guaranteed the
%windows before and after the peak velocity were the start and end points.
%
%The current implementation is the same as in the paper except:
%-the saccade start and end point are calculated from the velocity
%crossings + 3 samples.
%-(3) from above, but using a range.
%-(5) from above.

% Constants
global DEF
velThresh = 30;  % deg/sec % 50   10;
minVelDur = 4;  % consecutive msec required for velocity>thresh % 6 4 ; 
winSize = 50;  % ms
winSep = winSize + 25;  % ms gap between windows.
winStep = 1;  % ms
pThresh = 0.01;
winAmpThresh = 0.2;  % Win1 to Win2 in deg. Threshold seems low.
refrac = 75;  % minimum distance between 2 saccades in msec, changed from 150 ajs aug 20111
sacAmpRange = [3 27];  % Empirically. 8 < 6. 3 > 24.
% totAmpThresh = 3;

% Check optional arguments.
%states = DEF.EYE_State_undefined + zeros( [size(gazeData.sample),1], DEF.errorcodes);
if ~isempty(varargin) && ~isempty(varargin{1})
    gaze_bool = varargin{1};
else
    gaze_bool = true(size(gazeData.sample));
end

%Do not use sampleRate, but use gazeData.sample.
gazeMsec = gazeData.sample(gaze_bool);
velNorm = 1000*gazeData.velNorm(gaze_bool);

% Find putative saccades based on thresholding instantaneous velocity
velBool = velNorm > velThresh; %Samples with putative saccades
dVelBool = [nan;diff(velBool)]; %+1 is rising, -1 is falling.
sacStarts = find(dVelBool == 1); %Where 
sacStops = find(dVelBool == -1) - 1;

if ~isempty(velBool) && velBool(1)
%First sample was over threshold. i.e., record began mid-saccade.
    sacStops(1) = []; %Start sample unknown, so delete corresponding stop.
end
if ~isempty(velBool) && velBool(end) %Last sample was over threshold. i.e., record ended mid-saccade.
    sacStarts(end) = []; %Saccade end unknown, so get rid of corresponding start.
end

% Ignore saccades that are faster than minVelDur
shortSacBool = (gazeMsec(sacStops) - gazeMsec(sacStarts)) < minVelDur;
sacStarts = sacStarts(~shortSacBool);
sacStops = sacStops(~shortSacBool);

% Ignore saccades with bad states.
badstates = [DEF.EYE_State_blink DEF.EYE_State_invalid DEF.EYE_State_offscreen];
badStateBool = ismember(gazeData.state(gaze_bool), badstates);
badSacBool = true(size(sacStarts));
for sac_ix = 1:length(sacStarts)
    badSacBool(sac_ix) = any(badStateBool(sacStarts(sac_ix):sacStops(sac_ix)));
end
sacStarts = sacStarts(~badSacBool);
sacStops = sacStops(~badSacBool);

% Adam - Keep saccades with first and last positions distance in range.
degData = gazeData.degree(gaze_bool, :);
%[nDeg, nDim] = size(degData);
sacAmp = sqrt(sum((degData(sacStops,:) - degData(sacStarts,:)).^2, 2));

% Chad - Calculate maximum distance between all points and check amplitude.
% sacAmp = NaN(size(sacStarts));
% for sac_ix = 1:length(sacAmp)
%     sacdeg = double(degData(sacStarts(sac_ix):sacStops(sac_ix),:));
%     
%     %Reduce the points to the vertices of the convex outer hull.
%     if exist('convhull', 'builtin')==5
%         try
%             k = convhull(sacdeg, 'simplify', true);
%             sacdeg = sacdeg(k(1:end-1),:);
%         catch err
%             fprintf([err.identifier ' Saccade %i has no outer hull.\n']);
%         end
%     end
% 
%     %find the distance between all pairs.
%     K=sacdeg*sacdeg';
%     d=diag(K);
%     one=ones(length(d),1);
%     D=sqrt(d*one'+one*d'-2*K);
% 
%     sacAmp(sac_ix) = max(max(D));    
% end

sacAmpBool = sacAmp >= sacAmpRange(1) & sacAmp <= sacAmpRange(2);
sacStarts = sacStarts(sacAmpBool);
sacStops = sacStops(sacAmpBool);

% Step paired windows through each saccade.
% Move saccade start to the end of the first window of the pair that:
% (1) window_median vector amplitude > medAmpThresh
% (2) X || Y positions are statistically different
% (3) the next window-pair also meets these criteria.
newSacStarts = nan(size(sacStarts));
for sac_ix = 1:length(sacStarts) %Saccades that meet velocity criteria.
    
    % This function is called once per trial, but it need not be.
    if mod(sac_ix, 2500)==0
        fprintf('.\n')
    end
    if mod(sac_ix, 100)==0
        fprintf('.');
    end
    
    % Define the windowing.
    % Use gazeMsec
    
    % First win1 stops on first putative saccade sample
    % Last win1 stops on last putative saccade sample.
    % Interval is winStep (in msec). This might not line up perfectly.
    win1stops_ms = gazeMsec(sacStarts(sac_ix)):winStep:gazeMsec(sacStops(sac_ix));
    [~, win1stops] = ismember(win1stops_ms, gazeMsec);
    for ww = 1:length(win1stops)
        if win1stops(ww) == 0
            win1stops(ww) = find(abs(gazeMsec - win1stops_ms(ww)) == min(abs(gazeMsec - win1stops_ms(ww))), 1, 'first');
        end
    end
    win1stops = unique(win1stops);
    win1stops_ms = gazeMsec(win1stops);
    win1starts_ms = win1stops_ms - winSize;
    win2stops_ms = win1stops_ms + winSep;
    win2starts_ms = win2stops_ms - winSize;
    
    win1Bool = bsxfun(@ge, gazeMsec, win1starts_ms');
    win1Bool = win1Bool & bsxfun(@le, gazeMsec, win1stops_ms');
    win2Bool = bsxfun(@ge, gazeMsec, win2starts_ms');
    win2Bool = win2Bool & bsxfun(@le, gazeMsec, win2stops_ms');
    
    % Eliminate windows that include the very edge of gazeMsec (maybe out
    % of range) or that are not at all in gazeMsec
    badBool = win1Bool(1, :) | win1Bool(end, :) ...
        | win2Bool(1, :) | win2Bool(end, :) ...
        | ~any(win1Bool) | ~any(win2Bool);
    win1Bool = win1Bool(:, ~badBool);
    win2Bool = win2Bool(:, ~badBool);
    
    % Because the windows are not guaranteed to have the same number of
    % samples, I cannot vectorize the code. :(
    nWindows = size(win1Bool, 2);
    
    % Calculate the amplitude for each window pair.
    
    winok = false(nWindows, 1);
    for ww = 1:nWindows
        % median degree of each window
        degree1 = degData(win1Bool(:, ww), :);
        degree2 = degData(win2Bool(:, ww), :);
        xyxy = [median(degree1)' median(degree2)'];
        amp = sqrt(sum(diff(xyxy').^2));
        winok(ww) = amp > winAmpThresh;
        if winok(ww)  % If the amplitude is large enough
            winok(ww) = kstest2(degree1(:, 2), degree2(:, 2), pThresh);  % And the amps are stat. diff in y
            if ~winok(ww)
                winok(ww) = kstest2(degree1(:, 1), degree2(:, 1), pThresh);  % Or the amps are stat. diff in x.
            end
        end
    end
    doubleok = winok & [winok(2:end);false];  % Must have 2 consecutive samples.
    
    if any(doubleok)
        newSacStarts(sac_ix) = win1stops(find(doubleok, 1, 'first'));
    end
end
newSacStops = sacStops(~isnan(newSacStarts));
newSacStarts = newSacStarts(~isnan(newSacStarts));
if size(newSacStarts, 1) ~= size(newSacStops, 1)
    newSacStops = newSacStops';
end

% Group saccades that are within each other's refrac period.
groupedWithPrev = [false; diff(newSacStarts)<refrac];  % Original - If time between starts < refrac
% groupedWithPrev = [false;(newSacStarts(2:end)-newSacStops(1:end-1))<refrac];  % Adam - If time between n stop and n+1 start < refrac
groupStarts = find(diff(groupedWithPrev)==1);
groupStops = find(diff(groupedWithPrev)==-1);
if groupedWithPrev(end)
    groupStops = [groupStops;length(groupedWithPrev)];
end

% Only keep one sac from each group.
keepSac = true(size(newSacStarts));
for group_ix = 1:length(groupStarts)
    keepSac(groupStarts(group_ix):groupStops(group_ix)) = false;
    
    %Original - keep the first of the group. Fails on Trial 3.
%     keepSac(groupStarts(group_ix)) = true;
    
    %Adam - Only keep the saccade with the maximum summed-amp from each group.
    gSacStarts = newSacStarts(groupStarts(group_ix):groupStops(group_ix));
    gSacStops = newSacStops(groupStarts(group_ix):groupStops(group_ix));
    sacAmp = NaN(size(gSacStarts));
    for sac_ix = 1:length(sacAmp)
        sacAmp(sac_ix) = sum(velNorm(gSacStarts(sac_ix):gSacStops(sac_ix)));
    end
    [~, maxIdx] = nanmax(sacAmp);
    keepSac(groupStarts(group_ix) + maxIdx - 1) = true;
end
newSacStarts = newSacStarts(keepSac);
newSacStops = newSacStops(keepSac);

%Now we need to save the saccades.

%saccade strucutre is intended to relate eData to gazeData and to store 
%saccade-specific information.
saccade = struct(...
    'expTrial', [], ... %The experimental trial id #
    'gazeSampleStart', num2cell(newSacStarts), ... The sample number in gazeData
    'gazeSampleEnd', num2cell(newSacStops), ...
    'startPt', [], ... [X Y] screen coordinates in deg.
    'endPt', [], ... [X Y] screen coordinates in deg.
    'vectDir', [], ... rad2deg(atan2((endY-startY),(endX-startX)))
    'vectAmp', []); %norm of endPt - startPt
% 'trialIdx', [], ...The trial this saccade belongs to.
%     'trialSampleStart', [], ... The sample number, w.r.t. trial onset, that this saccade starts

for sac_ix = 1:length(saccade)
    start_ix = saccade(sac_ix).gazeSampleStart;
%     start_deg = degData(start_ix - floor(winSize/2) : start_ix + ceil(winSize/2),:);
    start_deg = degData(start_ix - min(3, start_ix-1) : start_ix, :);
    start_pt = nanmean(start_deg);
    stop_ix = min([saccade(sac_ix).gazeSampleEnd size(degData,1)-3]);
%     stop_deg = degData(stop_ix - floor(winSize/2) : stop_ix + ceil(winSize/2),:);
    stop_deg = degData(stop_ix:stop_ix + min(3, length(degData)-stop_ix), :);
    stop_pt = nanmean(stop_deg);
    %Adam used a different window for calculating start and stop points,
    %centred on the peak velocity.
    
    saccade(sac_ix).startPt = start_pt;
    saccade(sac_ix).endPt = stop_pt;
    saccade(sac_ix).vectDir = rad2deg(atan2(( stop_pt(2)-start_pt(2) ),( stop_pt(1)-start_pt(1) )));
    saccade(sac_ix).vectAmp = sqrt(sum((stop_pt - start_pt).^2));
end