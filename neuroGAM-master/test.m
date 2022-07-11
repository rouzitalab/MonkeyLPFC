cons_rate = zeros(trials, samples, channels);
for i = 1:trials
    for j = 1:channels
        cons_rate(i,:,j) = rate((i-1)*samples*channels + (j-1)*samples + 1 :(i-1)*samples*channels + (j-1)*samples + samples);
    end
end