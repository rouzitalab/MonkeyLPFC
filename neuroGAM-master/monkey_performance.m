function [learned, not_learned] = monkey_performance(data)
border = [98, 389, 465, 684, 757, 1146, 1255];
outcome = data.outcome;
flag = find(outcome~=-1);
outcome = outcome(flag);
monkey_perf = zeros(1,length(outcome));
cor = 0;
b = 1;
bin = 25;
for i=1:bin
if outcome(i)==0, cor=cor+1;
end
end
monkey_perf(1:bin) = 100 * (cor / bin);
i = bin + 1;
while i < (length(outcome) + 1)
if i == border(b)
cor = 0;
for j=1:bin
if outcome(i+j)==0, cor=cor+1;end
end
monkey_perf(i:i+bin) = 100 * (cor / bin);
i = i + bin;
b = rem(b,length(border))+1;
elseif outcome(i) == outcome(i-bin+1), monkey_perf(i)=monkey_perf(i-1);i=i+1;
elseif outcome(i) == 0, cor=cor+1; monkey_perf(i) = 100*(cor/bin); i=i+1;
else, cor = cor-1; monkey_perf(i)=100*(cor/bin); i=i+1;
end
end
learned = find(monkey_perf>70);
not_learned = find(monkey_perf<60);