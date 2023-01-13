%% Compare different relaxation mechanisms
clf; clc;

%%
clf; clc
figure(1)
t = linspace(3,300,300);

direct = @(A_dir, T) A_dir*T;

IntFunc = @(x)(x.^8.*exp(x))./((exp(x)-1).^2);
raman = @(A_ram, theta_D, T) ...
    arrayfun(@(xx)A_ram*xx.^9*integral(IntFunc, 0, 1/xx), ...
    T/theta_D);
local = @(A_loc, Delta_loc, T) ...
    A_loc*(exp(Delta_loc./(1*T)))./(exp(Delta_loc./(1*T))-1).^2;
orbach = @(A_orb, Delta_orb, T) ...
    A_orb*Delta_orb^3./(exp(Delta_orb./(1*T))-1);
orbach1 = @(A_orb, Delta_orb, T) ...
    A_orb*Delta_orb*exp(-Delta_orb./(1*T));
thermal = @(A_therm, tau_0, E_a, omega, T) ...
    A_therm*(2*tau_0*exp(E_a./(1*T)))./( ...
    1+omega^2*(tau_0*exp(E_a./(1*T))).^2);
cosech = @(A_tls, Delta_tls, T) ...
    A_tls*2./(exp(Delta_tls./T)-exp(-Delta_tls./T));
cross = @(A_cro, Delta_cro, T) ...
    A_cro./(1+exp(Delta_cro./T));
hold on

% remove boltzman constant, use K as unit

plot(t, direct(0.01, t), '--r', LineWidth=2);
plot(t, raman(20, 200, t), '-g', 'LineWidth', 2);
plot(t, local(1, 100, t), '--b', 'LineWidth',2)
plot(t, orbach(1e-6, 300, t), '--m', LineWidth=2)
plot(t, orbach1(0.12, 300, t), '-m')
plot(t, thermal(10,0.4,100,1, t), 'cyan', LineWidth=2)
plot(t, cosech(3,100, t), 'color','#800000', LineWidth=2)
plot(t, cross(10,10,t), 'color', 'r')
% plot(t, 7e-5*t.^2, 'k') % approximate Raman


legend('Direct', 'Raman', 'Local', 'Orbach', ...
    'Orbach approx', 'Thermal', 'Cosech', 'Cross', ...
    Location='best')

hold off
% set(gca, 'XScale', 'log', 'YScale', 'log')
%%
figure(2);
t = linspace(0,120,120);
d1 = direct(4, t);
r1 = raman(1.8e5,112, t);
t1i = d1+r1;
plot(t,t1i, 'k');
hold on
plot(t,d1, 'r--');
plot(t,r1,'g--');
set(gca,'xscale','log','yscale','log')

%%
t = linspace(0,300,100);
% y = direct(550,t) + cosech(3.9e5,35,t);
y = direct(21,t) + cosech(1.34e4,10,t);
figure(3);
plot(t,y)
set(gca, 'YScale', 'log')

%%

t = linspace(1,300,300);
y = thermal(6.8e10, 1.07e-12, 635, 9.5e9, t);
plot(t,y)
