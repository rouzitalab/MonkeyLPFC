function gazeFilterKernel = makeGazeFilterKernel(gazeFilter)
global DEF
switch gazeFilter.function
    %%%%%%%%%%%%%%%
    case 'gaussian'
        kernelSD = 3; %Standarddeviations
        kernelLength = 2*kernelSD*gazeFilter.param(1); 
        %%%aPars.smoothSpikeData.param(1) == sigma
        if gazeFilter.param(1) == 0 % no replacement: stay with the pattern
            gazeFilterKernel = 1;
        else
            x = -ceil(kernelLength/2):ceil(kernelLength/2);
            gazeFilterKernel = normpdf(x,0,gazeFilter.param(1));
        end
        
    %%%%%%%%%%%%%%%%%    
    case 'runavg'
        if gazeFilter.param(1) == 0 % no replacement: stay with the pattern
            gazeFilterKernel = 1;
        else
            gazeFilterKernel = ones(1,gazeFilter.param(1))/gazeFilter.param(1);
        end       
    %%%%%%%%%%%%%%%%    
%     case 'PSP'
%         kernelLength =  4*sum(aPars.smoothSpikeData.param);
%         x = 1:kernelLength;
%         % aPars.smoothSpikeData.param(1) and (2) == TAUg + TAUd
%         % (depolarisation/repolarisation)
%         if aPars.smoothSpikeData.param(1) == 0 % no replacement: stay with the pattern
%             sdfKernel = 1;
%         else
%             sdfKernel = (1-exppdf(x, aPars.smoothSpikeData.param(1))).*exppdf(x,aPars.smoothSpikeData.param(2));
%         end
%         
    %%%%%%%%%    
    otherwise
        error('%s: unknown filterGazeFunction',mfilename);
end

% cast,normalize to area of '1' (like the AUC of the delta-function) and
gazeFilterKernel = cast(gazeFilterKernel/sum(gazeFilterKernel),DEF.eyePosGazeFormat);