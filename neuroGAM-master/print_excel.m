clc;
filename = 'log.xlsx';
for ch = 1 : 32
    fprintf('\n\nCH %d: \n', ch);
    for cls = 1 : length(models{ch}.class)
        seq = cell2mat(models{ch}.class(cls));
        if seq(1) && ~seq(2) && ~seq(3)
            sac_ll = mean(models{ch}.testFit{cls,1}(:,3));
            sac_ve = mean(models{ch}.testFit{cls,1}(:,1));
            cel = ['D', int2str(ch+2)];
            writematrix(sac_ll,filename,'Sheet',1,'Range',cel);
            cel = ['J', int2str(ch+2)];
            writematrix(sac_ve,filename,'Sheet',1,'Range',cel);
            fprintf('Sac, LL: %d, VE: %d \n', sac_ll, sac_ve);
        elseif ~seq(1) && ~seq(2) && seq(3)
            rul_ll = mean(models{ch}.testFit{cls,1}(:,3));
            rul_ve = mean(models{ch}.testFit{cls,1}(:,1));
            cel = ['E', int2str(ch+2)];
            writematrix(rul_ll,filename,'Sheet',1,'Range',cel);
            cel = ['K', int2str(ch+2)];
            writematrix(rul_ve,filename,'Sheet',1,'Range',cel);
            fprintf('Rul, LL: %d, VE: %d \n', rul_ll, rul_ve);
        end
    end
end