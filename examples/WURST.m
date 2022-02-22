% This script is an example demonstrating the use of SBE for the WURST IR simulation.
% The WURST IR for 23Na was published in Madelin et al 2010. 

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

wurst              = struct;
deltaHz            = 0;
t1                 = 0.0;    % [s] not put into a struct to reduce sim time
flipAngle          = 180;    % [deg]
tau                = 10e-3;  % [s]
ampHz              = 250;
freqSweepRangeHz   = 2e3;  
freqModulation     = 2*pi*freqSweepRangeHz/tau;
wurst.pulse        = @(t) (t >= t1) .* (t <= t1+tau).* 2*pi*ampHz.*(1 - abs(sin(pi/2*(2*(t-t1)-tau)/tau)).^20);
wurst.grad         = @(t) (t >= t1) .* (t <= t1+tau) .* (-freqModulation .* (t - (t1+tau)/2));
wurst.trPhase      = @(t) (t >= t1) .* (t <= t1+tau) .* deg2rad(180);

sim                 = struct;
sim.time            = [t1, t1+tau];
sim.teval           = linspace(t1, t1+tau, 300);
sim.T0              = [1;zeros(14,1)];
sim.M0              = [0;0;1;0;0;1];
sim.options         = odeset('RelTol',1e-2,'AbsTol',1e-4);
ISTO                = SBE(saline, wurst, sim);

figure('Position',[216 284 559 421]), hold on; grid on; box on;
set(gca, 'LineWidth', 2, 'FontSize', 14);
plot(ISTO.tsim*1e3, ISTO.T(1,:),'LineWidth',2);
xlabel('Time (ms)'); ylabel('Longitudinal Magnetisation (a.u.)');
xlim([0,10]); ylim([-1,1]);title('Saline');

ISTO = SBE(agar, wurst, sim);
figure('Position',[216 284 559 421]), hold on; grid on; box on;
set(gca, 'LineWidth', 2, 'FontSize', 14);
plot(ISTO.tsim*1e3, ISTO.T(1,:), 'LineWidth', 2);
xlabel('Time (ms)'); ylabel('Longitudinal Magnetisation (a.u.)');
xlim([0,10]); ylim([-1,1]);title('Agar');

n = 100;
runTime = zeros(1,n);
for i = 1:n
    tic;
    ISTO = SBE(xanthan, wurst, sim);
    runTime(i) = toc; 
end
SBERunTime = [mean(runTime), std(runTime)] * 1e3    % [ms]


figure('Position',[216 284 559 421]), hold on; grid on; box on;
set(gca, 'LineWidth', 2, 'FontSize', 14);
plot(ISTO.tsim*1e3, ISTO.T(1,:),'LineWidth', 2);
xlabel('Time (ms)'); ylabel('Longitudinal Magnetisation (a.u.)');
xlim([0,10]); ylim([-1,1]);title('Xanthan');


