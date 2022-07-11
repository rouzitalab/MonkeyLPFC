x = [0.95, 1.05, 1.95, 2.05, 2.95, 3.05];
y = [6.534276, 6.074313, 34.472868, 24.140525, 239.11109, 137.67944];
x2 = [0.95, 1.95, 2.95];
x1 = [1.05, 2.05, 3.05];
y1 = [6.534276, 34.472868, 239.11109];
y2 = [6.074313, 24.140525, 137.67944];
err1 = [0.3738506,3.0307322,11.741554];
err2 = [0.9285516,2.261616,9.197304];
err = [0.3738506, 0.9285516, 3.0307322, 2.261616, 11.741554, 9.197304];
% y = [.085, .02, 0.035, 0.07, 0.07];
% err = [0.004, 0.007, 0.006, 0.005, 0.005];
% y = [.033, .03, 0.04, 0.038, 0.038];
% err = [0.0031, 0.0035, 0.002, 0.003, 0.003];
% y = [.001, .00075, 0.0017, 0.0016, 0.0016];
% err = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001];
errorbar(x1,y1,err1,'o', 'Color', [0 0.5 0], 'LineWidth', 2);
hold on;
errorbar(x2,y2,err2,'o', 'Color', 'red', 'LineWidth', 2);
% yline(0,'--');

% line([3],[0.04],'LineStyle', 'o');
hold on;
% plot(1, 0.01, '.', 'LineWidth', 10);
% plot(.1,.006,'.', 'Color', 'red', 'MarkerSize', 20);
% hold on;
% line([1 5], [0 0],'LineStyle','--','Color', 'black', 'LineWidth', 2);
ylim([-10 260]);
% xlim([0.075 .525]);
xticks([1,2,3]);
% yticks([0, 0.005, 0.01]);
xticklabels({'Saccade', 'Target', 'Rule'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',8)
ylabel('Normalized Likelihood of the Rule Model','FontSize',10);
xlabel('Channels','FontSize',10);
legend({'Learned', 'Unlearned'}, 'Location','northwest');
% daspect([13 1 1])
grid on;