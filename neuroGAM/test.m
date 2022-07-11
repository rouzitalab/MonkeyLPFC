for i=1:numel(models.class)
    for j=1:nvars
        other_vars = (j~=1:nvars & models.class{i});
        fprintf(['i = ', i, ' / j = ', j, '\n']);
        disp(size(cell2mat(x(other_vars))));
        disp(size(cell2mat(models.wts{i}(other_vars))'))
        isempty(models.wts{i}(other_vars)) || (strcmp(models.xtype{j},'event') && all(strcmp(models.xtype(other_vars),'event')))
    end
end
        
        