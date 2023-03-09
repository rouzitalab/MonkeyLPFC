s = wf.trialInfo.trialStartTimes;
e = wf.trialInfo.trialEndTimes;
sps = 1000; %sample per second
tr = 306;
l = e(tr)-s(tr);
figure();
spy(sua(:,ceil(s(tr)):ceil(e(tr))), 'k.')
asp = floor(l/(nUnits*4));
daspect([asp,1,2])
ylabel("Neuron")
xlabel("Time (s)")
x=[0, floor(.7*sps):floor(.25*sps):l];
xt = [-.7, 0:.25:l/sps-0.25];
xline(floor(.7*sps), "-", LineWidth=3, color="#0072BD", DisplayName="Target")
xline(floor(.95*sps), "-.", LineWidth=3, color="#0072BD", DisplayName="Color")
xline(floor(1.95*sps), "--", LineWidth=3, color="#0072BD", DisplayName="Delay")
xline(floor(2.2*sps), ":", LineWidth=3, color="#0072BD", DisplayName="Saccade")
xticks(x)
xticklabels(xt)
legend()