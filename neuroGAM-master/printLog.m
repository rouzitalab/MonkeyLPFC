filename = 'statistics_l.xlsx';
% if sess < 5, sheet = 1; base = sess*5-10;
% elseif sess < 8, sheet = 2; base = sess*5-25;
% else, sheet = 3; base = sess*5-40;
% end
% base = (run-1)*4;
col_name2 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
% if run < 7, sheet = 1;
% elseif run < 13, sheet = 2;
% elseif run < 19, sheet = 3;
% elseif run < 25, sheet = 4;
% else, sheet = 5;
% end
% if kind < 1, sheet = 1; base=0;
% else, sheet = 2; base=0;
% end
for ch = 1 : 32
    if (~isnan(models{ch}.bestmodel))
        seq = cell2mat(models{ch}.class((models{ch}.bestmodel)));
        selected = '';
        if seq(1), selected = [selected, 'S']; end
        if seq(2), selected = [selected, 'C']; end
%         if seq(3), selected = [selected, 'T']; end
        if seq(3), selected = [selected, 'R']; end
        if (ch*4) > 25, col_name = [col_name2(fix((ch*4)/26)), col_name2(mod(ch*4,26)+1)];
        else, col_name = col_name2(mod(ch*4,26)+1);
        end
        cel = [col_name,int2str(run+2)];
%         cel = [col_name(base+4), int2str(2*ch+1)];
        writematrix(selected,filename,'Sheet',1,'Range',cel);
        writematrix(selected,filename,'Sheet',2,'Range',cel);
    end
    for cls = 1 : length(models{ch}.class)
        seq = cell2mat(models{ch}.class(cls));
        if seq(1) && ~seq(2) && ~seq(3), ofst = 1;
        elseif ~seq(1) && seq(2) && ~seq(3), ofst=2;
        elseif ~seq(1) && ~seq(2) && seq(3), ofst=3;
%         elseif ~seq(1) && ~seq(2) && ~seq(3) && seq(4), ofst=4;
        else, ofst=0;
        end
        if (ofst)
            if (ch*4 + ofst - 4) > 25, col_name = [col_name2(fix((ch*4 + ofst - 4)/26)), col_name2(mod(ch*4 + ofst - 4,26)+1)];
            else, col_name = col_name2(mod(ch*4 + ofst - 4,26)+1);
            end
            for fold = 1 : 10
                cel = [col_name, int2str(run+2+fold-1)];
                writematrix(mean(models{ch}.testFit{cls,1}(fold,3)),filename,'Sheet',1,'Range',cel);
                cel = [col_name, int2str(run+2+fold-1)];
                writematrix(mean(models{ch}.testFit{cls,1}(fold,1)),filename,'Sheet',2,'Range',cel);
            end
        end
    end
%     cl = [col(sess*3 - 2), int2str(2*ch+1)];
%     writematrix(mean(models{ch}.testFit{1,1}(:,3)),filename,'Sheet',1,'Range',cl)
%     cl = [col(sess*3 - 1), int2str(2*ch+1)];
%     writematrix(mean(models{ch}.testFit{3,1}(:,3)),filename,'Sheet',1,'Range',cl)
%     cl = [col(sess*3 - 2), int2str(2*ch+2)];
%     writematrix(mean(models{ch}.testFit{1,1}(:,1)),filename,'Sheet',1,'Range',cl)
%     cl = [col(sess*3 - 1), int2str(2*ch+2)];
%     writematrix(mean(models{ch}.testFit{3,1}(:,1)),filename,'Sheet',1,'Range',cl)
    fprintf('Channel %d Written!\n', ch);
end
