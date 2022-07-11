function plotWeights(W, varargin)
if nargin > 1
    W_ref = varargin{1};
else
    W_ref = W{end};
end

if nargin > 2
    x_vec = varargin{2};
else
    x_vec = 1:length(W);
end

if nargin > 3
    nFactors = varargin{3};
else
    nFactors = size(W_ref, 2);
end

w_empty = cellfun(@isempty, W);
last_entries = find(~w_empty(1:end) & [w_empty(2:end);true]);
ref_entries = last_entries;
for bl_ix = 1:length(ref_entries)
    while any(any(isnan(W{ref_entries(bl_ix)})))
        ref_entries(bl_ix) = ref_entries(bl_ix) - 1;
    end
end

nSteps = length(W);

temp = cellfun(@size, W, num2cell(ones(length(W), 1)), 'UniformOutput', false);
nUnits = max([temp{:}]);

temp = cellfun(@size, W, num2cell(2*ones(length(W), 1)), 'UniformOutput', false);
nAllFacs = max([temp{:}]);

weights_out = nan(nSteps, nFactors*nUnits);
D_norm = nan(1, nSteps);
full_weights_out = nan(nSteps, nUnits*nAllFacs);
for step_ix = 1:nSteps
    
    % Swap axes and invert other C_outs so they are as close as possible to
    % W_ref
    if ~isempty(W{step_ix})
        
        if ~isnan(W_ref)
            this_W_ref = W_ref;
        else
            this_W_ref = W{ref_entries(find(last_entries >= step_ix, 1, 'first'))};
        end
        ref_vec = ones(1, size(this_W_ref, 1)) * this_W_ref;
        
        constrainedX = minDistBySwapAndInvert(this_W_ref', W{step_ix}');
        ctest = W{step_ix} * constrainedX;
%         unittraj(step_ix, :) = ones(1, size(ctest,1))*ctest(:, 1:nFactors);
        temp = ctest;
        temp(any(isnan(temp), 2), :) = 0;
        D_norm(step_ix) = norm(ref_vec - ones(1, size(ctest, 1)) * temp);
        
        full_weights_out(step_ix, :) = ctest(:);
        
        ctest = ctest(:, 1:nFactors);
        weights_out(step_ix, :) = ctest(:);
    end
end

%%
new_xvec = min(x_vec):min(diff(x_vec)):max(x_vec);
new_weights_out = nan(length(new_xvec), size(weights_out, 2));
lia = ismember(new_xvec, x_vec);
new_weights_out(lia, :) = weights_out;
for m_ix = find(~lia)
    new_weights_out(m_ix, :) = new_weights_out(m_ix-1, :);
end
h = imagesc(new_xvec, 1:size(new_weights_out,2), new_weights_out');
alpha = ~isnan(new_weights_out');
set(h, 'AlphaData', alpha);
set(gca, 'color', 'none');
colLim = max(abs(get(gca, 'CLim')));
caxis([-colLim colLim])
axis ij
xlabel('Trial #')
ylabel('Factor #')
set(gca, 'FontSize', 24)
set(gca, 'YTick', (1:nFactors)*nUnits - nUnits/2)
ytls = cell(1, nFactors);
horiz_lines = nan(nFactors, 2);
for yt_ix = 1:nFactors
    ytls{yt_ix} = ['F' num2str(yt_ix)];
    horiz_lines(yt_ix, :) = [yt_ix*nUnits yt_ix*nUnits];
end
set(gca, 'YTickLabel', ytls)
set(gca, 'LineWidth', 3)
hold on
plot([x_vec(1) x_vec(end)], horiz_lines',...
    'k--', 'LineWidth', 2);
% bl_switch = trial_vec(diff(bl_ix)>0);
% 
% for sw_ix = 1:length(bl_switch)
%     plot([bl_switch(sw_ix) bl_switch(sw_ix)], [1 size(weights_out, 2)],...
%         'k--', 'LineWidth', 4)
% end
hold off

%%
% figure;
% unittraj = nan(nSteps, 3);
% % Plot movement of 3D unit vector across winsteps
% plot(trial_vec', unittraj - 5*repmat([0 1 2], size(unittraj,1), 1), 'LineWidth', 5);
% axis tight
% set(gca, 'Color', 'none')
% set(gca, 'FontSize', 24')
% set(gca, 'YTick', fliplr(-5*[0 1 2]))
% set(gca, 'YTickLabel', {'F3','F2','F1'})
% xlabel('Trial #')
% hold on
% ylims = get(gca, 'YLim');
% for sw_ix = 1:length(bl_switch)
%     plot([bl_switch(sw_ix) bl_switch(sw_ix)], ylims,...
%         'k--', 'LineWidth', 4)
% end
% set(gca, 'LineWidth', 3)
% hold off
% 
% figure;
% plot(trial_vec, -log(D_norm), 'LineWidth', 5)
% axis tight
% hold on
% ylims = get(gca, 'YLim');
% for sw_ix = 1:length(bl_switch)
%     plot([bl_switch(sw_ix) bl_switch(sw_ix)], ylims,...
%         'k--', 'LineWidth', 4)
% end
% set(gca, 'LineWidth', 3)
% hold off
