function plot_Tm_vs_temp(temp,Tm,Color,Marker,MarkerSize,inverseT1)

arguments
    temp double
    Tm double
    Color = 'm'
    Marker {mustBeMember(Marker,{'o','+','*','.','x','-','|','s','d',...
        '^','v','>','<','p','h'})} = 'p'
    MarkerSize double = 100
    inverseT1 logical = false
end

if inverseT1
    y = 1./Tm;
    ylab = "1/T_m (1/ns)";
    ttl = "1/T_m versus Temperature";
else
    y = Tm;
    ylab = "T_m (ns)";
    ttl = "T_m versus Temperature";
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