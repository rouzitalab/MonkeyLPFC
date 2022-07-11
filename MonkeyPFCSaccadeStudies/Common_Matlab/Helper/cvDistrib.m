function testBool = cvDistrib(classVec, ncv, varargin)
%function testBool = cvDistrib(classVec, ncv[, false])
%last argument is optional. Set true to match number of training trials per
%class.

nTrials = length(classVec);
uqY = unique(classVec);
[nPerClass, ~] = hist(classVec, uqY);
randTID = randperm(nTrials);
randY = classVec(randTID);
testBool = false(nTrials, ncv);
minTest = min(ceil(nPerClass./ncv));
maxTrain = min(nPerClass - minTest);

if ~isempty(varargin) && varargin{1}
    %Each class is represented equally. This has a reduced training size and
    %uneven test size. For larger class, some cv's may have repeat test
    %trials.
    for cl_ix = 1:length(uqY)
        cl_id = uqY(cl_ix);
        cl_tr_id = randTID(randY == cl_id);
        
        nTest = max(length(cl_tr_id) - maxTrain, minTest);
        nTrain = length(cl_tr_id) - nTest;
        
        temp = nan(nPerClass(cl_ix), ncv);
        for cv_ix = 1:ncv
            train_ix = (cv_ix-1)*maxTrain+1:cv_ix*maxTrain;
            test_ix = train_ix(end)+1:train_ix(end)+nTest;
            while any(train_ix > nPerClass(cl_ix))
                train_ix(train_ix > nPerClass(cl_ix)) = train_ix(train_ix > nPerClass(cl_ix)) - nPerClass(cl_ix);
            end
            while any(test_ix > nPerClass(cl_ix))
                test_ix(test_ix > nPerClass(cl_ix)) = test_ix(test_ix > nPerClass(cl_ix)) - nPerClass(cl_ix);
            end
            test_id = cl_tr_id(test_ix);
            testBool(test_id, cv_ix) = true;
            
            temp(train_ix, cv_ix) = true;
        end
        
    end
else
    %Each class is represented proportionally.
    for cl_ix = 1:length(uqY)
        cl_id = uqY(cl_ix);
        cl_tr_id = randTID(randY == cl_id);
        
        edges = 0:length(cl_tr_id)/ncv:length(cl_tr_id);
        temp_id = 1:length(cl_tr_id);
        for cv_ix = 1:ncv
            test_id = cl_tr_id(temp_id > edges(cv_ix) & temp_id <= edges(cv_ix + 1));
            testBool(test_id, cv_ix) = true;
        end
    end
end