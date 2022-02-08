function plot_rmse_comparison(temperature,rmse,labels)

% rmse for one model is listed in one column
rmseSize = size(rmse);
for i = 1:rmseSize(2)
    scatter(temperature, rmse(:,i),100,'filled','^');
    hold on
end
hold off
legend(labels)
xlabel("Temperature (K)");
ylabel("RMSE")

fig = gcf;
if ~(fig.WindowStyle == "docked")
    set(fig,'position',[10,10,900,600]);
end

end