% DEER
%% (pi/2)x-tau-(pi)x-tau, excite only A spin
%%
syms Omega_A Omega_B omega_ee t tau real
load("operators.mat");
%%
H = Omega_A*Sz + Omega_B*Iz + omega_ee*SzIz;
sigma0 = Sz+Iz;
sigma1 = spin_evolve(sigma0, Sx*pi/2);
repr(sigma1)
%%
sigma2 = spin_evolve(sigma1, H*tau);
repr(sigma2)
%%
sigma3 = spin_evolve(sigma2, Sx*pi);
repr(sigma3)
%%
sigma4 = spin_evolve(sigma3, H*tau);
repr(sigma4)
% Fully refocused
% Answer: Sy


%% (pi/2)x-tau-(pi)x-tau, excite both A and B spin
H = Omega_A*Sz + Omega_B*Iz + omega_ee*SzIz;
sigma0 = Sz+Iz;
sigma1 = spin_evolve(sigma0, (Sx+Ix)*pi/2);
repr(sigma1)
%%
sigma2 = spin_evolve(sigma1, H*tau);
repr(sigma2)
%%
sigma3 = spin_evolve(sigma2, (Sx+Ix)*pi);
repr(sigma3)
%%
sigma4 = spin_evolve(sigma3, H*tau);
repr(sigma4)
% Answer: Iy*cos(omega_ee*tau) + Sy*cos(omega_ee*tau) - 2*SxIz*sin(omega_ee*tau) - 2*SzIx*sin(omega_ee*tau)




%% DEER three pulse
H = Omega_A*Sz + Omega_B*Iz + omega_ee*SzIz;
sigma0 = Sz; % only observe A spin
%% pi/2 x A
sigma1 = spin_evolve(sigma0, Sx*pi/2);
repr(sigma1)
%% H, t
sigma2 = spin_evolve(sigma1, H*t);
repr(sigma2)
%% pi, x, B
sigma3 = spin_evolve(sigma2, Ix*pi);
repr(sigma3)
%% H, tau-t
sigma4 = spin_evolve(sigma3, H*(tau-t));
repr(sigma4)
%% pi, x, A
sigma5 = spin_evolve(sigma4, Sx*pi);
repr(sigma5)
%% H, tau
sigma6 = spin_evolve(sigma5, H*tau);
repr(sigma6)
% Answer: Sy*cos(omega_ee*t) - 2*SxIz*sin(omega_ee*t)

