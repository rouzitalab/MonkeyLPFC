s = wf.trialInfo.trialStartTimes;
e = wf.trialInfo.trialEndTimes;
sps = 1000; %sample per second
tr1=1004;
l = e(tr1)-s(tr1);
cc = zeros(62,1);
sc = zeros(62,1);
for i = 1:62
cc(i) = full(sum(sua(i,floor(s(tr1)+.95*sps):floor(s(tr1)+1.95*sps))));
sc(i) = full(sum(sua(i,floor(s(tr1)+2.2*sps):ceil(e(tr1)))));
end
[o_cc, i_cc] = sort(cc, 'descend');
[o_sc, i_sc] = sort(sc, 'descend');
sin_ua = full(sua(i_cc,ceil(s(tr1)):ceil(e(tr1))));
% fprintf('Trial %i Top Units for Color Cue: %i with activities: %i\n',tr,i_cc(end-5:end),o_cc(end-5:end));
% fprintf('Trial %i Top Units for Saccade: %i with activities: %i\n',tr,i_sc(end-5:end),o_sc(end-5:end));
fr = zeros(size(sin_ua,1),floor(size(sin_ua,2)/50));
for i = 1:size(fr,2)
    fr(:,i)=sum(sin_ua(:,((i-1)*50)+1:i*50),2);
end
fr = fr/max(fr,[],'all');
figure();
colormap(hot)
% imagesc(fr)
fr_sm = imgaussfilt(fr,2);
imagesc(fr_sm, 'XData', [0,size(fr,2)])
hold on;
% daspect([asp,1,2])
ylabel("Neuron")
xlabel("Time (s)")
% x=[0, floor(.7*sps/50):floor(.25*sps/50):l/50];
% xt = [-.7, 0:.25:50/sps-0.25];
x=[0.1, 700/50:250/50:size(fr,2)+250/50];
xt = [-.7, 0:.25:l/sps-0.25];
% xline(floor(.7*sps/50), "-", LineWidth=3, color="#0072BD", DisplayName="Target")
% xline(floor(.95*sps/50), "-.", LineWidth=3, color="#0072BD", DisplayName="Color")
% xline(floor(1.95*sps/50), "--", LineWidth=3, color="#0072BD", DisplayName="Delay")
% xline(floor(2.2*sps/50), ":", LineWidth=3, color="#0072BD", DisplayName="Saccade")
xline(floor(.7*sps/50), "--", LineWidth=3, color="w", DisplayName="Target")
xline(floor(.95*sps/50), "--", LineWidth=3, color="w", DisplayName="Color")
xline(floor(1.95*sps/50), "--", LineWidth=3, color="w", DisplayName="Delay")
xline(floor(2.2*sps/50), "--", LineWidth=3, color="w", DisplayName="Saccade")

xticks(x)
xticklabels(xt)
colorbar()
% legend()
tr2=1509
figure();
colormap(hot)
l = e(tr2)-s(tr2);
cc = zeros(62,1);
sc = zeros(62,1);
for i = 1:62
cc(i) = full(sum(sua(i,floor(s(tr1)+.95*sps):floor(s(tr1)+1.95*sps))));
sc(i) = full(sum(sua(i,floor(s(tr1)+2.2*sps):ceil(e(tr1)))));
end
[o_cc, i_cc] = sort(cc, 'descend');
[o_sc, i_sc] = sort(sc, 'descend');
i_cc2 = zeros(62,1);
i_cc2(21:50)=i_cc(1:30);
i_cc2(1:10)=i_cc(53:end);
i_cc2(10:20)=i_cc(30:40);
i_cc2(50:end)=i_cc(40:52);
sin_ua = full(sua(i_cc2,ceil(s(tr2)):ceil(e(tr2))));
fr = zeros(size(sin_ua,1),floor(size(sin_ua,2)/50));
for i = 1:size(fr,2)
    fr(:,i)=sum(sin_ua(:,((i-1)*50)+1:i*50),2);
end
% fprintf('Trial %i Top Units for Color Cue: %i with activities: %i\n',tr,i_cc(end-5:end),o_cc(end-5:end));
% fprintf('Trial %i Top Units for Saccade: %i with activities: %i\n',tr,i_sc(end-5:end),o_sc(end-5:end));
fr = fr/max(fr,[],'all');
fr_sm = imgaussfilt(fr,2);
imagesc(fr_sm, 'XData', [0,size(fr,2)])
hold on;
asp = floor(l/(nUnits*4));
% daspect([asp,1,2])
ylabel("Neuron")
xlabel("Time (s)")
% x=[0, floor(.7*sps/50):floor(.25*sps/50):l/50];
% xt = [-.7, 0:.25:50/sps-0.25];
x=[0.1, 700/50:250/50:size(fr,2)+250/50];
xt = [-.7, 0:.25:l/sps-0.25];
% xline(floor(.7*sps/50), "-", LineWidth=3, color="#0072BD", DisplayName="Target")
% xline(floor(.95*sps/50), "-.", LineWidth=3, color="#0072BD", DisplayName="Color")
% xline(floor(1.95*sps/50), "--", LineWidth=3, color="#0072BD", DisplayName="Delay")
% xline(floor(2.2*sps/50), ":", LineWidth=3, color="#0072BD", DisplayName="Saccade")
xline(floor(.7*sps/50), "--", LineWidth=3, color="w", DisplayName="Target")
xline(floor(.95*sps/50), "--", LineWidth=3, color="w", DisplayName="Color")
xline(floor(1.95*sps/50), "--", LineWidth=3, color="w", DisplayName="Delay")
xline(floor(2.2*sps/50), "--", LineWidth=3, color="w", DisplayName="Saccade")
xticks(x)
xticklabels(xt)
colorbar()