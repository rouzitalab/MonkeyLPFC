data_traj = load('traj.mat');
traj_l = data_traj.traj_l;
traj_ul = data_traj.traj_ul;
y_l = data_traj.y_l;
y_ul = data_traj.y_ul;
t = 1:9;
color = ["red", "green"];
samples = length(traj_ul);
figure();
for i = 1 : samples
plot3(traj_ul(i,:,1), traj_ul(i,:,2), t, '-.', 'MarkerSize', 5, color='#808080')
% plot3(traj_ul(i,2,1), traj_ul(i,2,2), 2, 'o', 'MarkerSize', 10, color=color(y_ul(i)+1))
plot3(traj_ul(i,6,1), traj_ul(i,6,2), 6, '*', 'MarkerSize', 10, color=color(y_ul(i)+1))
% plot3(traj_ul(i,8,1), traj_ul(i,8,2), 8, 'h', 'MarkerSize', 10, color=color(y_ul(i)+1))
hold on;
end
% plot3(traj_ul(1,6,1), traj_ul(1,6,2), 6, 'x', 'DisplayName', 'Color Period', 'MarkerSize', 10, color=color(y_l(i)+1))
% plot3(traj_ul(samples,6,1), traj_ul(samples,6,2), 6, 'x', 'DisplayName', 'Color Period', 'MarkerSize', 10, color=color(y_l(i)+1))
% legend();
% xlim([-200 200])
% ylim([-200 200])
xlabel("Latent 1");
ylabel("Latent 2");
zlabel("Time (s/sample)");

samples = length(traj_l);
gcf = figure();
for i = 1 : samples
plot3(traj_l(i,:,1), traj_l(i,:,2), t, '-.', 'MarkerSize', 5, color='#CCCCCC')
plot3(traj_l(i,2,1), traj_l(i,2,2), 2, '*', 'MarkerSize', 10, color=color(y_l(i)+1))
plot3(traj_l(i,6,1), traj_l(i,6,2), 6, '*', 'MarkerSize', 10, color=color(y_l(i)+1))
% plot3(traj_l(i,8,1), traj_l(i,8,2), 8, 'h', 'MarkerSize', 10, color=color(y_l(i)+1))
hold on;
end
% plot3(traj_l(1,6,1), traj_l(1,6,2), 6, 'x', 'DisplayName', 'Color Period', 'MarkerSize', 10, color=color(y_l(i)+1))
% plot3(traj_l(samples,6,1), traj_l(samples,6,2), 6, 'x', 'DisplayName', 'Color Period', 'MarkerSize', 10, color=color(y_l(i)+1))
% legend();
xlim([-800 800])
ylim([-200 200])
xlabel("Latent 1");
ylabel("Latent 2");
zlabel("Time (s/sample)");
x = [-500 -500 500 500];
y = [-200 200 200 -200];
z = [6 6 6 6];
patch(x, y, [2 2 2 2], 'cyan', 'facealpha',0.2)
patch(x, y, [6 6 6 6], 'cyan', 'facealpha',0.2)
print(gcf, '-dpdf', '-r1200', 'traj3d_2.pdf')