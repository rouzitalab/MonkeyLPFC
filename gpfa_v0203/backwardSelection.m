function lowDSpike = backwardSelection(spike, label)

%% Randomize trials
numberOfRows = size(spike, 1);
newRowOrder = randperm(numberOfRows);
spike_shuf = spike(newRowOrder, :, :);
label_shuf = label(newRowOrder);
%% Single channel performances
tr_size = floor(0.6 * numberOfRows);
for ch = 1:size(spike,3)
    x_tr = squeeze(spike_shuf(1:tr_size, :, ch));
    y_tr = label(1:tr_size);
    x_ts = squeeze(spike_shuf(tr_size:end,:,ch));
    y_ts = label(tr_size:end);
    %% Train ensemble model
    Mdl = fitcensemble(x_tr,y_tr);
    % view(Mdl.Trained{1}.CompactRegressionLearner,'Mode','graph');
    %% Prediction
    y_pr = predict(Mdl,x_ts);
    %% Accuracy
    acc = sum(y_pr'==y_ts)/length(y_ts);
    fprintf("Single channel accuracy for # %i is: %.3f\n",[ch acc]);
end
