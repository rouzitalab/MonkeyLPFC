function [sentences, events] = getInfoFromWords(wordMatrix, varargin)
%[sentences, events] = getInfoFromWords(words, 'DIO', DIO, 'sentenceID', sentenceID)

%% Defaults and keyword arguments
[DIO, sentenceID] = get_default_DIO_sentenceID();
params = struct('DIO', DIO, 'sentenceID', sentenceID);
params = varg2params(varargin, params, {'DIO', 'sentenceID'});

% Remove fields with Bit
for fn = fieldnames(params.DIO)'
    if isfield(params.DIO, fn) && ~isempty(cell2mat(strfind(lower(fn), 'bit')))
        params.DIO = rmfield(params.DIO, fn);
    end
end

%% sentences
sentences = repmat(struct('type', [], 'value', [], 'time', []), 0, 1);
sentenceStartIndex = find(wordMatrix(:,1) == params.DIO.startSentenceWord);
sentenceStopIndex = find(wordMatrix(:,1) == params.DIO.stopSentenceWord);
if ~isempty(sentenceStartIndex)
    
    while sentenceStopIndex(1) < sentenceStartIndex(1)
        sentenceStopIndex(1) = [];
    end
    
    sentenceTypes = [fieldnames(params.sentenceID) struct2cell(params.sentenceID)];
    
    nSentences = length(sentenceStartIndex);
    sentences = repmat(struct('type', [], 'value', [], 'time', []), nSentences, 1);
    for s_ix = 1:nSentences
        st_ix = ismember([sentenceTypes{:, 2}], wordMatrix(sentenceStartIndex(s_ix) + 1, 1));
        sentences(s_ix).type = sentenceTypes{st_ix, 1};
        sentences(s_ix).value = wordMatrix(sentenceStartIndex(s_ix) + 2: sentenceStopIndex(s_ix) - 1, 1);
        sentences(s_ix).time = wordMatrix(sentenceStartIndex(s_ix), 2);
        
        wordMatrix(sentenceStartIndex(s_ix):sentenceStopIndex(s_ix), :) = nan;
    end
    wordMatrix(isnan(wordMatrix(:, 1)), :) = [];
    
end

%% events
events = repmat(struct('type', [], 'time', []), 0, 1);
nWords = size(wordMatrix, 1);
if nWords > 0
    
    wordTypes = [fieldnames(params.DIO) struct2cell(params.DIO)];
    events = struct('type', [], 'time', num2cell(wordMatrix(:, 2)));
    for w_ix = 1:nWords
        events(w_ix).type = wordTypes{ismember([wordTypes{:,2}], wordMatrix(w_ix, 1)), 1};
    end
end