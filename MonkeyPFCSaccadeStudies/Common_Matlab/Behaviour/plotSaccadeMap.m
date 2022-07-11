function plotSaccadeMap(ptb, varargin)
%function plotSaccadeMap(ptb)
% Plots a representation of the screen showing the start and stop point of
% the first saccade in all included trials. Each class gets its own colour.
% Trials with newOutcomeCode==0 are solid, with ==9 are dashed.

%TODO: Generalize to params.outcomeCodes instead of hard-coded 0 and 9.
params = varg2params(varargin,...
    struct(...
    'outcomeCodes', [0 9]));

exp_types = unique({ptb.eData.trial.expType});
n_exps = length(exp_types);
OCCs = params.outcomeCodes;
occ_ls = {'-', '--', '-.', '.'};

for exp_ix = 1:n_exps
    subplot(n_exps, 1, exp_ix)
    set(gca, 'Color', 'none');
    set(gca, 'FontSize', 24, 'LineWidth', 2)
    hold on
    
    exp_bool = strcmpi({ptb.eData.trial.expType}, exp_types{exp_ix});
    tr_class = cat(1, ptb.eData.trial(exp_bool).newClass);
    uq_classes = unique(tr_class);
    
    % Colour map for classes
    cmap = parula(length(uq_classes));
    
    for tr_ix = 1:length(ptb.eData.trial)
        this_trial = ptb.eData.trial(tr_ix);
        [~, class_id] = ismember(this_trial.newClass, uq_classes);
        [~, occ_id] = ismember(this_trial.newOutcomeCode, OCCs);
        
        sac = this_trial.saccades(1);
        plot(this_trial.eyePosition(sac.gazeSampleStart:sac.gazeSampleEnd, 2), ...
            this_trial.eyePosition(sac.gazeSampleStart:sac.gazeSampleEnd, 3), ...
            occ_ls{occ_id}, 'LineWidth', 3, 'Color', cmap(class_id, :));
    end

    hold off
    axis ij
    xlabel('X POS (deg)')
    ylabel('Y POS (deg)')
end