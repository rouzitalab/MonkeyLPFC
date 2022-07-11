function gazeData = correctForEyePosAdjustments(gazeData, adjustments)
%function gazeData = correctForEyePosAdjustments(gazeData, adjustments)
%adjustments columns:
%[xoffset xgain yoffset ygain xrotation yrotation tmarker]

%tmarker is the sample number found in gazeData.sample

% global DEF

% Eliminate adjustments that end in -1. I don't know what to do with them
% anyway.
if size(adjustments, 1) > 0
    adjustments = adjustments( adjustments(:, end) ~= -1 , : );
end
if size(adjustments,1) > 0
    % check for inconsistencies
    if any(diff(adjustments(:,end))==0)
        error('%s: Same sample appears twice in adjustments', mfilename);
    end
    if any(sum(adjustments(:,[5 6]))) % any changes in these columns
        warning('%s is not yet ready to cope with rotation adjustments', mfilename);
    end
    % Adjust
    for i = 1:size(adjustments,1)
        
        % Identify samples to be adjusted.
        block_bool = gazeData.sample >= adjustments(i, end);
        if i < size(adjustments, 1)
            block_bool = block_bool & gazeData.sample < adjustments(i+1, end);
        end
          
        % xValues
        gazeData.gaze(block_bool,1) = (gazeData.gaze(block_bool,1)*adjustments(i,2)) + adjustments(i,1);

        % yValues
        gazeData.gaze(block_bool,2) = (gazeData.gaze(block_bool,2)*adjustments(i,4)) + adjustments(i,3);
        
        %% for rotation: not yet implemented
        %     [x1b,Y0]=rotatePoint(x1a,y1a,caladj(5));
        %     [X0,y1b]=rotatePoint(x1a,y1a,caladj(6));
    end %i
end