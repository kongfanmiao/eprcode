clear;clc
load("operators.mat");

syms Omega_S omega_I A B real

H0 = Omega_S*Sz + omega_I*Iz + A*SzIz + B*SzIx;

eta_alpha = atan(-B/(A+2*omega_I));
eta_beta = atan(-B/(A-2*omega_I));
xi = (eta_alpha+eta_beta)/2;
eta = (eta_alpha-eta_beta)/2;
omega_12 = (omega_I+A/2)*cos(eta_alpha) - B/2*sin(eta_alpha);
omega_34 = (omega_I-A/2)*cos(eta_beta) + B/2*sin(eta_beta);
omega_plus = omega_12 + omega_34;
omega_minus = omega_12 - omega_34;

U = expm(-1i*(xi*Iy+eta*2*SzIy));
Hd = simplify(U*H0*U');
coef_Sz = simplify(trace(Sz*Hd));
coef_Iz = simplify(trace(Iz*Hd));
coef_SzIz = simplify(trace(SzIz*Hd));

syms a_iso t real
tmp = spin_evolve(Sy, Sz, Omega_S*t);
res = spin_evolve(tmp, SzIz, a_iso*t)
res_repr = repr(res)






