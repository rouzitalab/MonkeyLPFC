function gazeData = getEyeData(filename, eyeDataPath)

% gets the eyeposition Data collected with the SR_Research EyeLink in form
% of converted 'asc'-files (from'edf'). If there is no asc-file, this 
% function tries to convert the according edf-file,
% reads the produced asc-file and deletes it afterwards from the disc to save
% Disk space. One may toggle the deletion with the keyword 'deleteASCfile' [bool]
% in the analysis parameter file. The filled structure 'gazeData' contains fields
% for -samplenum, -pixels (horz-vert), -pupilsize, -velocities and resolutions.
% Thereafter the data may be corrected for recalibrations during the experiments.
% This will recalculate all the respective values.
% Depending on the kind of input data it
% performs a continuous analysis to determine the  actual 'state' of the
% eye (i.e., fixating,saccade smooth movement, etc.). In a third step this
% data is transferred to a trial basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global paths DEF
ascHeaderLines = 2;
validASCFileVersion = 1.0;

[~, name, ext] = fileparts(filename);

%Look for the ASC file.
ascPath = fullfile(paths.preprocessed, 'asc');
if ~isdir(ascPath)
    mkdir(ascPath);
end
ascFilename = [name '.asc'];
if ~exist(fullfile(ascPath, ascFilename), 'file')
    edfConversionOptions = ' -y -l -t -nst -nmsg -s -miss NaN -p '; % meaning: overwrite asc-file, leftEye, tab separated, no messages and eyeStates, plot NaN for missing values
    edf2ascCommandSting = [ '"' paths.edfConverterPathFile '"' edfConversionOptions '"' ascPath '" "'...
        fullfile(eyeDataPath, filename) '"'];
    [~] = system(edf2ascCommandSting);
end

% get the raw Gaze Positions
fprintf('Importing ASC-file... ');
fid = fopen(fullfile(ascPath, ascFilename));

firstLine = fgetl(fid);
if ~any(strcmpi({'fp_edf2asc' 'EDF2ASC'}, sscanf(firstLine,'%s',1)))
     error('%s: unable to read eyeFile: %s: wrong fileType',mfilename,filename);
end

if ~isequal(sscanf(firstLine,'%*s %*s %f',1),validASCFileVersion);
     error('%s: unable to read eyeFile: %s: wrong fileVersion',mfilename,filename);
end
[~] = fgetl(fid);

numDataCols  = [];
lineCount       = 0;
stopAfter       = 10;

% Find out how many columns the data have
while ~feof(fid)
    currLine = fgetl(fid);
    lineCount = lineCount + 1;
%     fprintf([currLine '\n']); %16086530	  472.5	   68.8	 4642.0	.
    if isempty(numDataCols)
        if ~isempty(currLine) && isnumeric(sscanf(currLine,'%d',1))
            numDataCols     = length(find(diff(isspace(currLine)) == -1));  % last character is always a dot
        end
        if lineCount > stopAfter
            error('%s: no first TimeStamp in file > %s',mfilename,filename);
            break;
        end
    else
        break;
    end
end

% start from the beginning
frewind(fid);

% jump over header
for i = 1:ascHeaderLines 
   [~] = fgetl(fid); 
end

% read the actual Matrix
% Keep consistent with DEF?
if numDataCols == 4
    fileData = textscan(fid, '%d %f32 %f32 %d16 %*s');
elseif numDataCols == 6
    fileData = textscan(fid, '%d %f32 %f32 %d16 %f32 %f32 %*s');
elseif numDataCols == 8
    fileData = textscan(fid, '%d %f32 %f32 %d16 %f32 %f32 %f32 %f32 %*s');
end
fclose(fid);

% Put the data into struture
if length(fileData)<4 % pre-allocate gazeData memory and check the basics of the input matrix
    error('No appropriate Data in asc-file: %s',filename);
end
gazeData.sample = fileData{1};  % Time stamps in msec.
gazeData.gaze   = [fileData{2},fileData{3}];
gazeData.pupil  = fileData{4};
if numDataCols > 4 
    gazeData.velocity   = [fileData{5},fileData{6}];
end
if numDataCols > 6
    gazeData.resolution = [fileData{7},fileData{8}];
end
%gazeData.sampleRate = DEF.eyeSampleRate;
gazeData.timeUnits = 'msec';
fprintf('Done.\n');