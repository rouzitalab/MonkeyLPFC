clc;
filename = 'temporal.xlsx';
sheet = 1;
c = 'BCDEFGHIJKLMNOPQ';
for ch = 1 : 32
    fprintf('\n\nCH %d: \n', ch);
    if ~isnan(models{ch}.bestmodel)
        for cls = 1 : length(models{ch}.class)
            seq = cell2mat(models{ch}.class(cls));
            if seq(1) && ~seq(2) && ~seq(3)
                ll = mean(models{ch}.testFit{cls,1}(:,3));
                cel = [c(w*4 + 1), int2str(ch+2)];
                writematrix(ll,filename,'Sheet',sheet,'Range',cel);
            elseif ~seq(1) && seq(2) && ~seq(3)
                ll = mean(models{ch}.testFit{cls,1}(:,3));
                cel = [c(w*4 + 2), int2str(ch+2)];
                writematrix(ll,filename,'Sheet',sheet,'Range',cel);
            elseif ~seq(1) && ~seq(2) && seq(3)
                ll = mean(models{ch}.testFit{cls,1}(:,3));
                cel = [c(w*4 + 3), int2str(ch+2)];
                writematrix(ll,filename,'Sheet',sheet,'Range',cel);
            elseif ~seq(1) && seq(2) && seq(3)
                ll = mean(models{ch}.testFit{cls,1}(:,3));
                cel = [c(w*4 + 4), int2str(ch+2)];
                writematrix(ll,filename,'Sheet',sheet,'Range',cel);
            end
        end
    end
end