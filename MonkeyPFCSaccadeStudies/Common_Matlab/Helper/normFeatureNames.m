function [featureNames, xTickNames, re_ix] = normFeatureNames(featureNames)

%Replace names
remapFNames = {...
    'baseline,cue,delay,response' 'concat';
    'baseline,cue,delay' 'concat';
    'baseline' 'base';
    'target' 'targ';
    'visual' 'cue';
    'cue' 'cue';
    'fullDelay' 'LongDelay';
    'motor' 'resp';
    'response' 'resp';
    'fullTrial' 'AvgFRate';
    'smooth' 'ntraj';
    'pca' 'PCA';
    'fa' 'FA';
    'gpfa' 'GPFA';
    'canon_cue' 'CanonCue';
    'canon_delay' 'CanonDelay';
    'canon_resp' 'CanonResp'};
    
tempNames = featureNames;
for f_ix = 1:length(featureNames)
    if any(strcmpi(featureNames{f_ix}, remapFNames(:,1)))
        r_ix = strcmpi(featureNames{f_ix}, remapFNames(:,1));
        tempNames{f_ix} = remapFNames{r_ix, 2};
    end
end
clear f_ix r_ix
featureNames = tempNames;

%Reorder names
reorder_fnames = {'base', 'targ', 'cue', 'delay', 'LongDelay', 'resp', 'AvgFRate', 'concat',...
    'ntraj', 'PCA', 'FA', 'GPFA', 'dPCA',...
    'CanonCue', 'CanonDelay', 'CanonResp'};
re_ix = nan(1, length(reorder_fnames));
for rx = 1:length(reorder_fnames)
    rebool = strcmpi(featureNames, reorder_fnames{rx});
    if any(rebool)
        re_ix(rx) = find(rebool);
    end
end
re_ix(isnan(re_ix)) = [];
featureNames = featureNames(re_ix);

featToLabel = {...
    'resp' 'resp.';
    'AvgFRate' 'avg.\newline rate';
    'concat' 'concat.';
    'ntraj' 'full\newline traj.';
    'CanonCue' 'canon.\newline   cue';
    'CanonDelay' 'canon.\newline delay';
    'CanonResp' 'canon.\newline  resp.'};

xTickNames = featureNames;
for x_ix = 1:length(xTickNames)
    xbool = strcmpi(featToLabel(:, 1), xTickNames{x_ix});
    if any(xbool)
        xTickNames{x_ix} = featToLabel{xbool, 2};
    end
end

