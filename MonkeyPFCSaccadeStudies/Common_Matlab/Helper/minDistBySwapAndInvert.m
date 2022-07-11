function constrainedX = minDistBySwapAndInvert(base, test)
%Use the solution to base = test'*X as a seed to find a
%constrainedX, where each column has only 1 entry {-1,0,1}

temp = test; temp(:, any(isnan(test))) = 0;
X = temp'\base';  %new_gpfa = tmp'*X;

constrainedX = zeros(size(X));

while any(sum(constrainedX, 1) == 0)
    % Columns that haven't already been locked in
    bRemainingCols = sum(constrainedX, 1) == 0;
    bRemainingRows = sum(constrainedX, 2) == 0;
    orig_cols = find(bRemainingCols);
    orig_rows = find(bRemainingRows);
    
    tmpX = X;
    tmpX(~bRemainingRows, :) = [];
    tmpX(:, ~bRemainingCols) = [];
    
    %Find the highest scoring row for each column
    [~, i] = max(abs(tmpX));
    
    %See if any columns share the best row
    [c, ia, ic] = unique(i);
    if length(c) < sum(bRemainingCols)
        
        %if any cols share their best row,
        %find the col that, when selected, maximizes the
        %total score.
        dup_ind = setdiff(1:length(i), ia);
        dup_row = i(dup_ind(1));
        rem_row = setdiff(1:size(tmpX, 1), dup_row);
        dup_cols = find(i == dup_row);
        
        this_score = nan(1, length(dup_cols));
        for col_ix = 1:length(dup_cols)
            this_col = dup_cols(col_ix);
            rem_cols = setdiff(dup_cols, this_col);
            [~, rem_rowix] = max(abs(tmpX(rem_row, rem_cols)));
            rem_inds = sub2ind(size(tmpX), rem_row(rem_rowix), rem_cols);
            this_score(col_ix) = sum(abs(...
                [tmpX(dup_row, this_col) tmpX(rem_inds)]));
        end
        
        % Lock the column with the max score into
        % constrainedX
        [~, col_ix] = max(this_score);
        full_col_ix = orig_cols(dup_cols(col_ix));
        full_row_ix = orig_rows(dup_row);
        constrainedX(full_row_ix, full_col_ix) = ...
            sign(X(full_row_ix, full_col_ix));
    else
        keep_inds = sub2ind(size(X), orig_rows(i), orig_cols');
        constrainedX(keep_inds) = sign(X(keep_inds));
    end
end

% figure;
% colLims = max(max(abs(base)));
% colLims = [-colLims colLims];
% subplot(3, 1, 1)
% imagesc(base')
% caxis(colLims);
% subplot(3, 1, 2)
% imagesc(test')
% caxis(colLims);
% subplot(3, 1, 3)
% imagesc(test' * constrainedX)
% caxis(colLims);