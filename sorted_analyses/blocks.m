clc
f = NaN;
s = NaN;
block = 1;
l = 0;
for i = 1:length(labels)
    if labels(i)~=f
        if isnan(f)
            f=labels(i);
        elseif labels(i)~=s
            if isnan(s)
                s=labels(i);
            else
                fprintf('Block #%i - Length: %i\n',block,i-l);
                fprintf('Rules: %i, %i\n\n',f,s);
                block = block + 1;
                f=labels(i);
                s=NaN;
                l=i;
            end
        end
    end
end
