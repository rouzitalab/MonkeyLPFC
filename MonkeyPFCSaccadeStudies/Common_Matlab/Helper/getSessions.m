function sessionsOut = getSessions()
%TODO: try readtable
global paths
if isempty(paths)
    my_paths;
end

%% Load the file into a cell array for headers and a cell array for data.
xlfiles = dir(fullfile(paths.sessionInfo, '*.xls*'));
if ~isempty(xlfiles)
    %Ignore temp ~$ files.
    if length(xlfiles)>1
        isgood = false(size(xlfiles));
        for xx = 1:length(xlfiles)
            isgood(xx) = ~(xlfiles(xx).name(1) == '~');
        end
        xlfiles = xlfiles(isgood);
    end
    [~, txt, raw] = xlsread(fullfile(paths.sessionInfo, xlfiles(1).name));
    headers = txt(1,:);
    headers = normalize_headers(headers);
    data = raw(2:end,:);
    
    if ismac
        dateAdd = datenum('01-Jan-1904');
    else
        dateAdd = datenum('30-Dec-1899');
    end
    nRows = size(data, 1);
    for col_ix = 1:length(headers)
        if strcmpi(headers{col_ix}, 'Date')
            data(:, col_ix) = num2cell(cellfun(@plus, data(:, col_ix), repmat(num2cell(dateAdd), nRows, 1)));
        elseif strcmpi(headers{col_ix}, 'nevdir')
            for row_ix = 1:nRows
                data{row_ix, col_ix} = num2str(data{row_ix, col_ix});
                if length(data{row_ix, col_ix}) < 6
                    data{row_ix, col_ix} = ['0' data{row_ix, col_ix}];
                end
            end
        end
    end
else
    csvfiles = dir(fullfile(paths.sessionInfo, '*.csv'));
    fid = fopen(fullfile(paths.sessionInfo, csvfiles(1).name), 'r');
    headers = fscanf(fid, '%s', 1);
    headers = strsplit(headers, ',');
    headers = normalize_headers(headers);
    content = textscan(fid, '%s');
    fclose(fid);
    data = cell(length(content{1}), length(headers));
    for row_ix = 1:length(content{1})
        this_row = content{1}{row_ix};
        data(row_ix, :) = strsplit(this_row, ',');
    end
    %Convert the columns to the appropriate data type (date, str, num)
    for col_ix = 1:length(headers)
        temp = strcmpi(data(:, col_ix), 'none');
        if any(temp)
            data(temp, col_ix) = num2cell(nan(sum(temp), 1));
        end
        if strcmpi(headers{col_ix}, 'Date')
            data(:, col_ix) = num2cell(datenum(data(:, col_ix), 'yyyy-mm-dd'));
        elseif any(strcmpi(headers{col_ix}, {'IsGood','NTrials','NCentreOut','NCueColour','NDiagonal'}))
            data(:, col_ix) = num2cell(cellfun(@str2num, data(:, col_ix)));
        elseif strcmpi(headers{col_ix}, 'nevdir')
            for row_ix = 1:size(data, 1)
                if length(data{row_ix, col_ix}) < 6
                    data{row_ix, col_ix} = ['0' data{row_ix, col_ix}];
                end
            end
        end
    end
end
clear xlfiles xx isgood txt raw csvfiles fid content row_ix this_row
%% Convert headers and data to sessions struct array
for hh = 1:length(headers)
    sessionsOut.(headers{hh}) = [];
end
clear hh

%Prepare the output variable
nSessions = size(data,1);
sessionsOut = repmat(sessionsOut, nSessions, 1);
for sess_ix = 1:nSessions
    thisData = data(sess_ix, :);
    for hh = 1:length(headers)
        fn = headers{hh};
        cellVal = thisData{strcmpi(headers, fn)};
        if strcmpi(fn, 'isgood')
            sessionsOut(sess_ix).(fn) = logical(cellVal);
        else
            sessionsOut(sess_ix).(fn) = cellVal;
        end
    end
end


end

function headers = normalize_headers(headers)
    %Normalize header format
    for h_ix = 1:length(headers)
        this_header = headers{h_ix};
        this_header = strtrim(this_header);
        this_header = this_header(~isspace(this_header));
        headers{h_ix} = lower(this_header);
    end

end