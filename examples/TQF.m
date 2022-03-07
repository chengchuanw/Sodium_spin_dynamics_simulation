% This script is an example demonstrating the use of SBE for the 6-step TQF simulation.
% The TQF sequence was originally published in 
% Ileana Hancu, Fernando E. Boada and Gary X. Shen, Three-Dimensional
% Triple-Quantum–Filtered 23Na Imaging of In Vivo Human Brain (1999)

clear; 

xanthan            = struct;
xanthan.J0         = 319; 
xanthan.J1         = 28.2; 
xanthan.J2         = 28.1; 
xanthan.T1         = 17.3e-3;
xanthan.T2         = 3.7e-3;
xanthan.deltaHz    = 0;
xanthan.omegaQHz   = 61.8;

% prepare pulse seq
tqf                = struct;
t1                 = 0.0;    % [s] 
tau                = 0.5e-3; 
t2                 = tau + 3e-3;
t3                 = t2 + tau + 0.4e-3;
flipAngle          = 90;    % [deg]
tqf.pulse              = @(t) (t >= t1) .* (t <= t1+tau) .* (deg2rad(flipAngle)/tau)...
                        + (t >= t2) .* (t <= t2+tau) .* (deg2rad(flipAngle)/tau)...
                        + (t >= t3) .* (t <= t3+tau) .* (deg2rad(flipAngle)/tau); % [rad/s]
tqf.grad           = @(t) 0*t;
ph                 = [30,90,150,-150,-90,-30]; % phase cycling steps

% prepare solver
sim                = struct;
sim.time           = [t1, t2, t3, 12e-3];
noFrames           = 300;
sim.teval          = linspace(0, 12e-3, noFrames);
sim.T0             = [1;zeros(14,1)];
sim.options        = odeset('RelTol',1e-3,'AbsTol',1e-5);
T10                = zeros(noFrames, 6);
T11a               = zeros(noFrames, 6);
T11s               = zeros(noFrames, 6);
T20                = zeros(noFrames, 6);
T21a               = zeros(noFrames, 6);
T21s               = zeros(noFrames, 6);
T22a               = zeros(noFrames, 6);
T22s               = zeros(noFrames, 6);
T30                = zeros(noFrames, 6);
T31a               = zeros(noFrames, 6);
T31s               = zeros(noFrames, 6);
T32a               = zeros(noFrames, 6);
T32s               = zeros(noFrames, 6);
T33a               = zeros(noFrames, 6);
T33s               = zeros(noFrames, 6);


for j = 1:6
    tqf.trPhase    = @(t) (t >= t1) .* (t <= t1+tau) .* deg2rad(  ph(j)   )...
                        + (t >= t2) .* (t <= t2+tau) .* deg2rad(  ph(j) + 90)...
                        + (t >= t3) .* (t <= t3+tau) .* deg2rad(0);   % rad
    ISTO           = SBE(xanthan, tqf, sim);
    % detector phase toggling
    T10(:,j)       = ISTO.T(1,:);                   
    T11a(:,j)      = (-1)^(j+1) * ISTO.T(2,:);
    T11s(:,j)      = (-1)^(j+1) * imag(ISTO.T(3,:));
    T20(:,j)       = ISTO.T(4,:);
    T21a(:,j)      = (-1)^(j+1) * ISTO.T(5,:);
    T21s(:,j)      = (-1)^(j+1) * imag(ISTO.T(6,:));
    T22a(:,j)      = (-1)^(j+1) * imag(ISTO.T(7,:));
    T22s(:,j)      = (-1)^(j+1) * ISTO.T(8,:);
    T30(:,j)       = ISTO.T(9,:);
    T31a(:,j)      = (-1)^(j+1) * ISTO.T(10,:);
    T31s(:,j)      = (-1)^(j+1) * imag(ISTO.T(11,:));
    T32a(:,j)      = (-1)^(j+1) * imag(ISTO.T(12,:));
    T32s(:,j)      = (-1)^(j+1) * ISTO.T(13,:);
    T33a(:,j)      = (-1)^(j+1) * ISTO.T(14,:);
    T33s(:,j)      = (-1)^(j+1) * imag(ISTO.T(15,:));
end
figure, subplot(321),plot(ISTO.tsim*1e3, T11a); grid on;
xlabel('Time (ms)'); ylabel('T11a'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(322),plot(ISTO.tsim*1e3, T11s); grid on;
xlabel('Time (ms)'); ylabel('imag(T11s)'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(323),plot(ISTO.tsim*1e3, T21a); grid on;
xlabel('Time (ms)'); ylabel('T21a'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(324),plot(ISTO.tsim*1e3, T21s); grid on;
xlabel('Time (ms)'); ylabel('imag(T21s)'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(325),plot(ISTO.tsim*1e3, T22a); grid on;
xlabel('Time (ms)'); ylabel('imag(T22a)'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(326),plot(ISTO.tsim*1e3, T22s); grid on;
xlabel('Time (ms)'); ylabel('T22s'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');

figure, subplot(321),plot(ISTO.tsim*1e3, T31a); grid on;
xlabel('Time (ms)'); ylabel('T31a'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(322),plot(ISTO.tsim*1e3, T31s); grid on;
xlabel('Time (ms)'); ylabel('imag(T31s)'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(323),plot(ISTO.tsim*1e3, T32a); grid on;
xlabel('Time (ms)'); ylabel('imag(T32a)'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(324),plot(ISTO.tsim*1e3, T32s); grid on;
xlabel('Time (ms)'); ylabel('T32s'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(325),plot(ISTO.tsim*1e3, T33a); grid on;
xlabel('Time (ms)'); ylabel('T31a'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
subplot(326),plot(ISTO.tsim*1e3, T33s); grid on;
xlabel('Time (ms)'); ylabel('imag(T31s)'); xlim([0, 12]);
legend('30^\circ','90^\circ','150^\circ','-150^\circ','-90^\circ','-30^\circ');
