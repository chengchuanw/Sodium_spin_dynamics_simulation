% This script is an example demonstrating the use of SBE for the SIRFLA simulation,
% which is a sequence proposed in Stobbe and Beaulieu 2005.

clear; 

saline             = struct;
saline.J0          = 8.9;  % spectral density [Hz]
saline.J1          = 8.9; 
saline.J2          = 8.9; 
saline.omegaQHz    = 0;    % residual quadrupolar interaction 
saline.deltaHz     = 0;

agar            = struct;       
agar.type       = '8% agar';
agar.J0         = 250;         
agar.J1         = 45.4; 
agar.J2         = 19.3;
agar.omegaQHz   = 0;            
agar.deltaHz    = 0;

xanthan         = struct;
xanthan.J0      = 319; 
xanthan.J1      = 28.2; 
xanthan.J2      = 28.1; 
xanthan.T1      = 17.3e-3;
xanthan.T2      = 3.7e-3;
xanthan.omegaQHz = 61.8; 
xanthan.deltaHz = 0;


% RF excitation
t1              = 0.0;     % seconds
flipAngle       = 180;     % degrees
tau             = 10e-3;   % seconds
sirfla.pulse           = @(t) 0.0 ...
     + (t >= t1) .* (t <= t1+tau) .* (deg2rad(flipAngle)/tau);
sirfla.grad        = @(t) 0.0 * t;
sirfla.trPhase        = @(t) 0.0 ...   
    + (t >= t1) .* (t <= t1+tau) .* deg2rad(0);

% numerical integration of the Equations
sim.T0          = zeros(15,1); sim.T0(1) = 1;
sim.time        = [t1, t1+tau];
sim.teval       = linspace(t1, t1+tau, 100);
sim.options     = odeset('RelTol',1e-2,'AbsTol',1e-3);    
ISTO            = SBE(saline, sirfla, sim);
longPulse       = figure; subplot(131), hold on; grid on; box on;
set(gca, 'LineWidth', 2, 'FontSize', 16);
plot(ISTO.tsim*1e3, ISTO.T(1,:), 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');
xlim([0,10]); ylim([-1,1]); title('Saline');

ISTO            = SBE(agar, sirfla, sim);
subplot(132); hold on; grid on; box on;
set(gca, 'LineWidth', 2, 'FontSize', 16);
plot(ISTO.tsim*1e3, ISTO.T(1,:), 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');
xlim([0,10]); ylim([-1,1]); title('Agar');

n = 100;
runTime = zeros(1,n);
for i = 1:n
    tic;
    ISTO = SBE(xanthan, sirfla, sim);
    runTime(i) = toc; 
end
SBERunTime = [mean(runTime), std(runTime)] * 1e3

subplot(133); hold on; grid on; box on;
set(gca, 'LineWidth', 2, 'FontSize', 16);
plot(ISTO.tsim*1e3, ISTO.T(1,:), 'LineWidth', 1);
xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');
xlim([0,10]); ylim([-1,1]); title('Xanthan');


