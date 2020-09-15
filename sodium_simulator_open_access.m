% Simulation of sodium tensor operator evolution in the rotating frame
% Version 1.0
% Author: Chengchuan Wu
% Last modified date: 15/9/2020

% simulation parameters
% fast-motion, isotropic parameter set
%{
J0         = 10.64; 
J1         = 9.24; 
J2         = 9.24; 
omegaQHz   = 0;
DeltaHz    = 0; 
%}
% slow-motion, anisotropic parameter set
%
J0         = 625;   % spectral densities, Hz
J1         = 50.8; 
J2         = 30.4; 
omegaQHz   = 21;    % static quadrupolar coupling, Hz
DeltaHz    = 0;     % off-resonance frequency, Hz
%}

noFrames        = 3e3;
duration        = 40e-3;  % seconds
tsim            = linspace(0.0, duration, noFrames);
M0              = zeros(15,1); M0(1) = 1; % Assume rho_eq = T10


% RF excitation
t1              = 0.0;    % seconds
tau             = 0.5e-3; % seconds
flipAngle       = 90;     % degrees
pulse           = @(t) 0.0 ...
    + (t >= t1) .* (t <= t1+tau) .* (deg2rad(flipAngle)/tau);
tr_phase        = @(t) 0.0 ...   
    + (t >= t1) .* (t <= t1+tau) .* deg2rad(180);

%---
% Differential Equation
sodium          = @(J0, J1, J2, Delta, omegaQ, omega1, phi, M, M0, t)[(-2/5).*J1+(-8/5).*J2,omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),0,0,0,0,0,(-4/5).*J1+(4/5).*J2,0,0,0,0,0,0;
    (-1).*omega1(t).*sin(-phi(t)),(-3/5).*J0+(-1).*J1+(-2/5).*J2,Delta*1i,0,0,1i.*(3/5).^(1/2).*omegaQ,0,0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,0,0;
    (1i*(-1)).*omega1(t).*cos(phi(t)),Delta*1i,(-3/5).*J0+(-1).*J1+(-2/5).*J2,0,1i.*(3/5).^(1/2).*omegaQ,0,0,0,0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,0;
    0,0,0,(-2).*J1+(-2).*J2,3.^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*3.^(1/2).*omega1(t).*cos(phi(t)),0,0,0,0,0,0,0,0,0;
    0,0,1i.*(3/5).^(1/2).*omegaQ,(-1).*3.^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J0+(-1).*J1+(-2).*J2,Delta*1i,(1i*(-1)).*omega1(t).*cos(phi(t)),omega1(t).*sin(-phi(t)),0,0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,0;
    0,1i.*(3/5).^(1/2).*omegaQ,0,(1i*(-1)).*3.^(1/2).*omega1(t).*cos(phi(t)),Delta*1i,(-1).*J0+(-1).*J1+(-2).*J2,omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,0,0;
    0,0,0,0,(1i*(-1)).*omega1(t).*cos(phi(t)),(-1).*omega1(t).*sin(-phi(t)),(-1).*J0+(-2).*J1+(-1).*J2,2*Delta*1i,0,0,0,0,1i.*omegaQ,0,0;
    0,0,0,0,(-1).*omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),2*Delta*1i,(-1).*J0+(-2).*J1+(-1).*J2,0,0,0,1i.*omegaQ,0,0,0;
    (-4/5).*J1+(4/5).*J2,0,0,0,0,0,0,0,(-8/5).*J1+(-2/5).*J2,6.^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*6.^(1/2).*omega1(t).*cos(phi(t)),0,0,0,0;
    0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,1i.*(2/5).^(1/2).*omegaQ,0,0,(-1).*6.^(1/2).*omega1(t).*sin(-phi(t)),(-2/5).*J0+(-1).*J1+(-3/5).*J2,Delta*1i,(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),(5/2).^(1/2).*omega1(t).*sin(-phi(t)),0,0;
    0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,(1i*(-1)).*6.^(1/2).*omega1(t).*cos(phi(t)),Delta*1i,(-2/5).*J0+(-1).*J1+(-3/5).*J2,(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),0,0;
    0,0,0,0,0,0,0,1i.*omegaQ,0,(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),(-1).*(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J0+(-1).*J2,2*Delta*1i,(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),(3/2).^(1/2).*omega1(t).*sin(-phi(t));
    0,0,0,0,0,0,1i.*omegaQ,0,0,(-1).*(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),2*Delta*1i,(-1).*J0+(-1).*J2,(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t));
    0,0,0,0,0,0,0,0,0,0,0,(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),(-1).*(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J1+(-1).*J2,3*Delta*1i;
    0,0,0,0,0,0,0,0,0,0,0,(-1).*(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),3*Delta*1i,(-1).*J1+(-1).*J2] * M...
    + [2/5.*J1+8/5.*J2;0;0;0;0;0;0;0;4/5.*J1-4/5.*J2;0;0;0;0;0;0];



% numerical integration of the Equations
options     = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,15));    % solver accuracy
[T, M]      = ode45(@(t, M)sodium(J0, J1, J2, 2*pi*DeltaHz, 2*pi*omegaQHz, pulse, tr_phase, M, M0, t), tsim, M0, options); 


T10         = M(:,1);
T11a        = M(:,2);
T11s        = M(:,3) * -1i; % take its imaginary value

%
T20         = M(:,4);
T21a        = M(:,5);
T21s        = M(:,6) * -1i;
T22a        = M(:,7) * -1i;
T22s        = M(:,8);

%
T30         = M(:,9);
T31a        = M(:,10);
T31s        = M(:,11) * -1i;
T32a        = M(:,12) * -1i;
T32s        = M(:,13);
T33a        = M(:,14);
T33s        = M(:,15) * -1i;

%
T           = T * 1e3; % ms

%

%
figure, set(gca,'fontsize',12);hold on;
plot(T, T10); 
plot(T, T11a); 
plot(T, T11s); 
legend('T10','T11a','Im\{T11s\}');
title('Rank-1 Tensor Operators')
xlabel('time (msec)'); ylabel('amplitude (a.u.)');
%
figure, set(gca,'fontsize',12);hold on;
plot(T, T20); 
plot(T, T21a); 
plot(T, T21s); 
plot(T, T22a); 
plot(T, T22s); 
legend('T20','T21a','Im\{T31s\}','Im\{T22a\}','T22s');
title('Rank-2 Tensor Operators')
xlabel('time (msec)'); ylabel('amplitude (a.u.)');
%
figure, set(gca,'fontsize',12);
subplot(211), hold on; set(gca,'fontsize',12);
title('Rank-3 Tensor Operators')
plot(T, T30); 
plot(T, T31a); 
plot(T, T31s);
legend('T30','T31a','Im\{T31s\}');
ylabel('amplitude (a.u.)');
subplot(212), hold on; set(gca,'fontsize',12);
plot(T, T32a); 
plot(T, T32s);
plot(T, T33a); 
plot(T, T33s); 
legend('Im\{T32a\}','T32s','T33a','Im\{T33s\}');
xlabel('time (msec)'); ylabel('amplitude (a.u.)');
