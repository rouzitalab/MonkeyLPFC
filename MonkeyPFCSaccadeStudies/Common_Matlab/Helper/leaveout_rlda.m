function pred_y = leaveout_rlda(X, y, out_ind)
cvIdx = setdiff(1:length(y), out_ind);
trainX = X(cvIdx, :);
testX = X(out_ind, :);
trainY = y(cvIdx);

% %SVM
% SVMStruct = svmtrain(trainX, trainY, 'autoscale', false);
% blockY(cv_ix) = svmclassify(SVMStruct, testX);

% rLDA
model = ml_train({[ones(size(trainX, 1), 1) trainX], trainY},...
    {'lda',...
    'lambda', [],...
    'regularization', 'auto'});  %'auto' 'shrinkage' 'independence'
predictions = ml_predict([ones(size(testX, 1), 1) testX], model);
[~, max_ix] = max(predictions{2}, [], 2);
pred_y = predictions{3}(max_ix);
end