function plot_T1Tf_comparison(temperatue,T1,T1Coef,Tf,TfCoef)

xlabel('Temperature (K)');
yyaxis left
scatter(temperatue, T1, 100, 'filled','o');
hold on
scatter(temperatue, Tf, 100, 'o');
ylabel("T_1 & T_f (\mus)")
yyaxis right
scatter(temperatue, T1Coef*100, 100, 'filled', 'd');
scatter(temperatue, TfCoef*100, 100, 'd');
ylabel('Percentage (%)')
legend('T_1', 'T_f', 'T_1 percentage', 'T_f percentage', "Location","best");
hold off

fig = gcf;
if ~(fig.WindowStyle == "docked")
    set(fig,'position',[10,10,900,600]);
end

end



