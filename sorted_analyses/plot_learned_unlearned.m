unit=1:62;
figure();
for i=1:62
    plot(i-0.1,l(i),marker="^", color="#77AC30")
    hold on;
    plot(i+0.1,u(i),marker="v", color="#D95319")
    hold on;
    xline(i, ":", LineWidth=1, color="k")
    hold on;
end
xlim([0 63]);
xlabel("Units");
ylabel("Normalized Performance");
legend("HP", "LP");
title("Normalized Single Unit Rule Decoding Performance");