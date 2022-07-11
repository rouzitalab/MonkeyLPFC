function eData = consolidateBehavAndGaze(eData, gazeData)
%eData = consolidateBehavGaze(eData, gazeData)
%Replaces eData.trial(n).eyePosition with the better data from gazeData.
%NOTE: It now pulls the whole trial data, not just data 100 msec before
%flipScreen(2)
% preBuffer = 100;  % msec before flipScreen(2) to include in saccade detection.
% endBuffer = 200;  % msec extra data to detect trial-ending saccades.

% e2g = 1; %Factor to convert eData times to gazeData times.
% if strcmpi(eData.trialTimeUnits, 's') && strcmpi(gazeData.units,'samples')
%     e2g = gazeData.sampleRate;
% end

% We need timing information from each eData.trial, but they have very
% different structures.
% See PTB/checkPTBCharacteristics.m

% We want to match the trial start to a sample in 

if isfield(eData.trial(1), 'eyeSyncTime')
    fileFormat = 0;
elseif isfield(eData.trial(1), 'eyeSyncStartTime')
    fileFormat = 1;
else
    fileFormat = 2;
end

nTrials = length(eData.trial);
%[eData.trial.eyePositionTimeUnits] = deal('s');
if isfield(eData.trial(1), 'eyeposition')
    eData.trial = rmfield(eData.trial, 'eyeposition'); %[eData.trial.eyeposition] = deal([]); %Empty unneeded field.
end
[eData.trial.nSaccades] = deal(0);
[eData.trial.saccades] = deal([]); %Will fill below.
[eData.trial.eyePosition] = deal([]); %Will fill below.

if isnan(eData.trial(end).stopTime)
    eData.trial(end).stopTime = max(eData.trial(end).flipScreen);
end
gaze_start_ms = [eData.trial.eyeSyncTime];
gaze_stop_ms = gaze_start_ms + int32(1000*[eData.trial.stopTime]);  %  + endBuffer;

fprintf('Finding saccades in %i trials. (.=50)', length(eData.trial));
for t_ix = 1:nTrials
    if mod(t_ix,1000)==1
        fprintf('\n');
    end
    if mod(t_ix,50)==0
        fprintf('.');
    end
    
    this = eData.trial(t_ix);
    
    % It might be desirable to not extract the whole eye data for the trial
    % and instead only include the data during the period in which the
    % saccades are important. For example, after this.flipScreen(2);
%     if numel(this.flipScreen) > 1
%         gaze_skip_ms = int32(1000*this.flipScreen(2));
%         this_gaze_bool = gazeData.sample >= (gaze_start_ms(t_ix) + gaze_skip_ms - preBuffer)...
        this_gaze_bool = gazeData.sample >= gaze_start_ms(t_ix) ...
            & gazeData.sample < gaze_stop_ms(t_ix);

        this.saccades = deriveSaccades(gazeData, this_gaze_bool);
        if ~isempty(this.saccades)
            [this.saccades.expTrial] = deal(this.expTrial);
        end
        this.nSaccades = length(this.saccades);
        
        % Get gazeData sample times in seconds since trial start
        
        thisGazeTime = double(gazeData.sample(this_gaze_bool) - gaze_start_ms(t_ix)) / 1000;

        % Save the gazeData to the trial. We're saving degrees but replacing
        % eyePosition.
        this.eyePosition = [thisGazeTime ...
            double(gazeData.degree(this_gaze_bool, :)) ...
            double(gazeData.state(this_gaze_bool, :))];
        
%         previewGazePerTrial(this);
%         [~] = waitforbuttonpress;
%     end
    eData.trial(t_ix) = this;
end

if isfield(eData.trial(1), 'response')
    eData.trial = rmfield(eData.trial, 'response');
end

fprintf(' Done.\n');