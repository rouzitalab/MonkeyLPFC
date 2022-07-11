filename = 'statistics2.xlsx';
% if sess < 5, sheet = 1; base = sess*5-10;
% elseif sess < 8, sheet = 2; base = sess*5-25;
% else, sheet = 3; base = sess*5-40;
% end
% base = (run-1)*4;
if kind == 0, s_c = 'B'; r_c = 'D';
else, s_c = 'C'; r_c='E';
end
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
    if kind==0, sel=int2str(10*(ch-1)+2);
    else, sel=int2str(10*(ch-1)+3);
    end
    if (~isnan(models{ch}.bestmodel))
        seq = cell2mat(models{ch}.class((models{ch}.bestmodel)));
        selected = '';
        if seq(1), selected = [selected, 'S']; end
        if seq(2), selected = [selected, 'T']; end
%         if seq(3), selected = [selected, 'T']; end
        if seq(3), selected = [selected, 'R']; end
        cel = ['H',sel];
%         cel = [col_name(base+4), int2str(2*ch+1)];
        writematrix(selected,filename,'Sheet',1,'Range',cel);
    else, cel = ['H',sel]; writematrix('None',filename,'Sheet',1,'Range',cel);
    end
    for fold = 1 : 10
        cel = [s_c, int2str(10*(ch-1)+fold+1)];
        writematrix(models{ch}.testFit{1,1}(fold,3)*models{ch}.testFit{1,1}(fold,1),filename,'Sheet',1,'Range',cel);
        cel = [r_c, int2str(10*(ch-1)+fold+1)];
        writematrix(models{ch}.testFit{3,1}(fold,3)*models{ch}.testFit{3,1}(fold,1),filename,'Sheet',1,'Range',cel);
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
