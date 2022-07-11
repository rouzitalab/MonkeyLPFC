% SdT_Arraymap 30/05/2013
%This script reads an array info file .cmp and finds the correct electrode placement.
% The wire is on the right side, these arrays were implanted in the left side.
%IF the array is split in blocks, block B starts at 33 and so on.

%% Read the .cmp map file
[map, path] = uigetfile('.cmp', 'Select Array map');
cd(path)
fid = fopen(map,'r');

L = fgetl(fid); % skip first line
lineCount = 1;
while ~feof(fid)
    lineCount = lineCount + 1;  % count the lines that have been read ii
    L = fgetl(fid);    % read currentLine
    if isempty(L)
        continue
    end
    firstWord = sscanf(L, '%s', 1);

    % search for keyword
    switch firstWord
        case 'Subject'
            continue
            
        case 'Hemisphere'
            arrayInfo.hemi        = sscanf(L, '%*s %s',1); 
            continue
            
        case 'WireOrientation'
            arrayInfo.wireOri     = sscanf(L, '%*s %d',1); 
            continue
            
        case 'ImplantOrientation'
            arrayInfo.implantOri = sscanf(L, '%*s %d',1); 
            continue
            
        case 'electrodeSpacing'
            arrayInfo.elSpace    = sscanf(L, '%*s %d',1); 
            continue
            
        case 'Cerebus' % whole line is Cerebus mapping 
            arrayInfo.map    = textscan(fid, '%f %f %s %f');
            arrayInfo.map{1} = arrayInfo.map{1}+1; % increase to index that starts with '1'
            arrayInfo.map{2} = arrayInfo.map{2}+1;
        otherwise
            error('Unknown keyword in %s',mfilenme)

    end
end
fclose(fid);


%% Output a matrix with electrode placement

array = zeros(10);

for k = 1:length(arrayInfo.map{1,3})
    
    if strcmp(arrayInfo.map{1,3}(k), 'A')
        array(arrayInfo.map{1,2}(k), arrayInfo.map{1,1}(k)) =  arrayInfo.map{1,4}(k);
        
    elseif strcmp(arrayInfo.map{1,3}(k), 'B')
        array(arrayInfo.map{1,2}(k), arrayInfo.map{1,1}(k)) =  arrayInfo.map{1,4}(k) + 32;
        
    elseif strcmp(arrayInfo.map{1,3}(k), 'C')
        array(arrayInfo.map{1,2}(k), arrayInfo.map{1,1}(k)) =  arrayInfo.map{1,4}(k) + 64;
        
    else
        array(arrayInfo.map{1,2}(k), arrayInfo.map{1,1}(k)) = NAN;
    end
    
end

array = flipud(array); % Flip array up down. This is confirmed the correct orientation.
