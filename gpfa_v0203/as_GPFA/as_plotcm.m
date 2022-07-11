% as_plotcm
%ajs 0911

sum_cm=sum(cm);
for i = 1:length(sum_cm)
    cm_norm(:,i)=cm(:,i)/sum_cm(i);
end
figure,imagesc(cm_norm)