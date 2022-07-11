function params = varg2params(varg, params, varargin)
% params = varg2params(varg, params, varargin)
% This function will typically be called at the very beginning of a
% function to parse its varargin input if it is expected that the varargin
% input (appearing here as varg) will be either key-value-pairs or a
% structure.
% The second argument (params) must contain the structure onto which the
% input arguments will be placed. If empty then a new structure will be
% made.
% The optional third argument is a cell array of parameter names to
% include, parameters not included will be ignored.
debug = false;
if ~isempty(varg)
    dropped_fns = {};
    %varg may be a structure
    if isstruct(varg{1})
        fnames = fieldnames(varg{1});
        for fn = 1:length(fnames)
            if isempty(varargin) || any(strcmpi(fnames{fn}, varargin{1})) %No include list, or fn in include list
                params.(fnames{fn}) = varg{1}.(fnames{fn});
            else
                dropped_fns = cat(2, dropped_fns, fnames{fn});
            end
        end
        clear fnames fn
    
    %varg is a cell array of key,value pairs.
    else
        if mod(length(varg),2)~=0
            error('Specify additional arguments as a structure or name-value pairs.');
        end
        
        for pp = 1:length(varg)/2 %for each kvp
            % If there is no include list, or the key is in the include
            % list
            if isempty(varargin) || any(strcmpi(varg{2*pp-1}, varargin{1}))
                params.(varg{2*pp-1}) = varg{2*pp};
            else
                dropped_fns = cat(2, dropped_fns, varg{2*pp});
            end
        end
        clear pp
        
    end
    if ~isempty(dropped_fns) && debug
        warning(['Ignored parameters ' strjoin(dropped_fns, ',') '. Check spelling.\n']);
    end
elseif ~isempty(varargin)
    %If we did not receive inputs, return a struct anyway.
    if isempty(params)
        params = struct();
    end
    for f_ix = 1:length(varargin{1})
        if ~isfield(params, varargin{1}{f_ix})
            params.(varargin{1}{f_ix}) = [];
        end
    end
end