function plot_T1_vs_temp(temp,T1,Color,Marker,MarkerSize,inverseT1)

arguments
    temp double
    T1 double
    Color = 'm'
    Marker {mustBeMember(Marker,{'o','+','*','.','x','-','|','s','d',...
        '^','v','>','<','p','h'})} = 'p'
    MarkerSize double = 100
    inverseT1 logical = false
end

if inverseT1
    y = 1./T1;
    ylab = "1/T_1 (1/\mus)";
    ttl = "1/T_1 versus Temperature";
else
    y = T1;
    ylab = "T_1 (\mus)";
    ttl = "T_1 versus Temperature";
end

scatter(temp, y, MarkerSize,Marker,Color,'filled')
xlabel("Temperature (K)");
ylabel(ylab);
title(ttl)

fig = gcf;
if ~(fig.WindowStyle == "docked")
    set(fig,'position',[10,10,900,600]);
end

end