function plot_TmTf_comparison(temperatue,Tm,TmCoef,Tf,TfCoef)

xlabel('Temperature (K)');
yyaxis left
scatter(temperatue, Tm, 100, 'filled','o');
hold on
scatter(temperatue, Tf, 100, 'o');
ylabel("T_m & T_f (ns)")
yyaxis right
scatter(temperatue, TmCoef*100, 100, 'filled', 'd');
scatter(temperatue, TfCoef*100, 100, 'd');
ylabel('Percentage (%)')
legend('T_m', 'T_f', 'T_m percentage', 'T_f percentage', "Location","best");
hold off

fig = gcf;
if ~(fig.WindowStyle == "docked")
    set(fig,'position',[10,10,900,600]);
end

end