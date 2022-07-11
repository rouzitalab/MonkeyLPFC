function plotBehavyy(behav)
%function plotBehav(behav, eData)
%behav is a matrix with columns:
%[goodStep trIdx blockIdx classes targets distractors ...
%   newOutcomeCode outcomeCode tpOut fpOut fnOut tnOut];
pmarg = 0.01; %If probability approaches 0 or 1, fix by this amount.
arrRad = 0.05; %Arrow radius in normalized figure units.
arrPosD = [-0.5 -1.5]; %Arrow position rel to block middle in units of d'

triali = behav(:,2);
tp = behav(:,9); fp=behav(:,10); fn=behav(:,11); tn=behav(:,12);
%         cnfMat = [tp fp; ...
%                   fn tn];

%http://www.birmingham.ac.uk/Documents/college-les/psych/vision-laboratory/sdtintro.pdf
tpr = tp./(tp+fn); tpr(tpr<pmarg)=pmarg; tpr(tpr>1-pmarg)=1-pmarg;
fpr = fp./(fp+tn); fpr(fpr<pmarg)=pmarg; fpr(fpr>1-pmarg)=1-pmarg;
ztpr = norminv(tpr, 0, 1);
zfpr = norminv(fpr, 0, 1);
d = ztpr - zfpr;
c = -(ztpr+zfpr)/2;

%70% for both classes is d=1.0488.
%Highest possible is 6.93, but effectively 4.65 for 99%

%Other measures of performance.
% sens = tp ./ (tp+fp);
% spec = tn ./ (tn+fn);
% balAcc = (sens+spec)/2;
% informedness = sens+spec-1;

labelNames = {'d''' 'Bias'};
labelColor = {'b' 'r'};
ylims = {[-2.0 5], [-1.5 1.5]};
[yyax, h1, h2] = plotyy(triali,[d ones(size(d))*1.05],...
    triali, [c -ones(size(c)) ones(size(c))]);
box off
set(h1(1), 'LineWidth', 3)
set(h1(2), 'LineWidth', 1, 'LineStyle', '--', 'Color', 'b')
set(h2(1), 'LineWidth', 3, 'Color', 'r')
set(h2(2), 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
set(h2(3), 'LineWidth', 1, 'LineStyle', '--', 'Color', 'r')
for yi=1:2
    set(yyax(yi), 'LineWidth', 1.5, 'FontSize', 14)
    set(get(yyax(yi), 'YLabel'),...
        'String', labelNames{yi}, 'FontSize', 14, 'Color', labelColor{yi})
    set(yyax(yi), 'XLim', [triali(1) triali(end)])
    set(yyax(yi), 'YLim', ylims{yi})
end
xlabel('Trial Index');

%Plot each block's colour-direction rule.
goodStep = logical(behav(:,1));
blockIdx = behav(:,3);
classes = behav(:,4);
uqBlocks = unique(blockIdx(goodStep));

axes(yyax(1))
hold on
axes(yyax(2))
hold on
for bb=1:length(uqBlocks)
    blockStart = triali(find(blockIdx==uqBlocks(bb), 1, 'first'));
    plot(yyax(1), [blockStart blockStart], ylims{1}, 'k--');
end
axes(yyax(1))
hold off
axes(yyax(2))
hold off

for bb = 1:length(uqBlocks)
    bix = uqBlocks(bb);
    blockBool = blockIdx==bix;
    uqClasses = unique(classes(blockBool));
    [targDirDeg, cueCol] = class2DirCol(uqClasses);
    
    %Plot an arrow for each target/colour combination.
    %On top of a white box.
    %The arrow origin will be the middle of the block at d=-2 or -4
    arrPosIx = double(mod(bb,2)==0)+1;
    arrOrig = [nanmean(triali(blockBool)), arrPosD(arrPosIx)];
    [origxf, origyf] = ds2nfu(yyax(1), arrOrig(1), arrOrig(2));
    for arix = 1:length(targDirDeg)
        [endx, endy] = pol2cart(deg2rad(targDirDeg(arix)), arrRad);
        endxf = origxf + endx;
        endyf = origyf + endy;
        annotation('arrow', [origxf endxf], [origyf endyf], ...
            'Color', cueCol{arix}, 'LineWidth', 4);
    end
end

