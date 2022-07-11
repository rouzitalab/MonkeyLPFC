function plotCfn(cfnMat, clssDirVec)

if isempty(clssDirVec)
    clssDirVec = -90:45:225;
    clssDirVec(clssDirVec < 0) = clssDirVec(clssDirVec < 0) + 360;
end

nArr = length(clssDirVec);
imagesc(100*bsxfun(@rdivide, cfnMat, sum(cfnMat,2)))
set(gca, 'XTickLabel', {})
set(gca, 'YTickLabel', {})
xlabel('PREDICTED DIRECTION', 'FontSize', 14)
ylabel('TRUE DIRECTION', 'FontSize', 14)
caxis([0 100]);

colorbar


axPos = get(gca, 'Position');

arrLen = 0.035; %?

%Divide the axis limits into nArr points
for xy = 1:2
    start = axPos(xy);
    axLen = axPos(2+xy);
    blockSize = axLen/nArr;
    
    for dd = 1:nArr
        arr_start = start + blockSize*(dd-1) + blockSize/2;
        if xy==1
            [endx, endy] = pol2cart(deg2rad(clssDirVec(dd)), arrLen);
            arrX = [arr_start-endx/2 arr_start+endx/2];
            arrY = [0.05-endy/2 0.05+endy/2];
        else
            [endx, endy] = pol2cart(deg2rad(clssDirVec(end-dd+1)), arrLen);
            arrX = [0.05-endx/2 0.05+endx/2];
            arrY = [arr_start-endy/2 arr_start+endy/2];
        end
        annotation('arrow', arrX, arrY, 'LineWidth', 2);
    end
end