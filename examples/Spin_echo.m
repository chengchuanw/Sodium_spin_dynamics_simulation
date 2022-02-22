% This script is an example demonstrating the CPMG SE sequence simulation
% for 23Na spin dynamics
clear; 

% prepare sample 
agar            = struct;       
agar.type       = '8% agar';
agar.J0         = 250;          % spectral density [Hz]
agar.J1         = 45.4; 
agar.J2         = 19.3;
agar.T1         = 21.7e-3;      % apparent T1,T2 [s]
agar.T2         = 5.95e-3;      
agar.omegaQHz   = 0;            % residual quadrupolar interaction [Hz]
agar.deltaHz    = 0;            % local field off resonance


% prepare pulse seq
spinEcho              = struct;
t1                    = 0.0;    
t2                    = 5e-3;
tau                   = 1e-3;

spinEcho.pulse        = @(t) (t >= t1) .* (t < t1+tau).* deg2rad(90) / tau...
                        + (t >= t2) .* (t < t2+tau).* deg2rad(180) / tau;
spinEcho.trPhase      = @(t) (t >= t1) .* (t < t1+tau).* deg2rad(180)...
                        + (t >= t2) .* (t < t2+tau).* deg2rad(0);
tGrad1                = t1+tau;
tGrad2                = t2+tau;
tGradDur              = 4e-3;
TR                    = 20e-3;

sim.options           = odeset('RelTol', 1e-2, 'AbsTol', 1e-3);
isochromats           = [];         % initialisation
sim.time              = [t1 t2 TR]; % specify the pulse action time to avoid ode45 mistakenly skipping them
sim.teval             = linspace(t1, TR, 200); % evaluate the solution at these given time points 
sim.T0                = [1; zeros(14,1)];

for j = -5:4    % 10 isochromats
        spinEcho.grad   = @(t) (t >= tGrad1) .* (t < tGrad1+tGradDur) * 4*pi / (10*tGradDur) * j...
            +(t >= tGrad2) .* (t < tGrad2+2*tGradDur) * 4*pi / (10*tGradDur) * j;
        ISTO        = SBE(agar, spinEcho, sim, 'default', 'full');
        isochromats = cat(1, isochromats, ISTO.Mxy);
end
isochromats = mean(isochromats,1);
figure, subplot(211), set(gca, 'FontSize', 14, 'LineWidth', 1),
plot(ISTO.tsim,abs(isochromats), 'LineWidth', 2); xlabel('Time [sec]'); ylabel('M_{xy} Magnitude (a.u.)');
grid on, box on;
subplot(212), set(gca, 'FontSize', 14, 'LineWidth', 1),
fplot(@spinEcho.pulse, [0,20e-3],'LineWidth', 2); ylim([-4000,4000]);
yticks([])
yyaxis right; fplot(@spinEcho.grad, [0,20e-3],'LineWidth', 2); ylim([0,8000]);
yticks([])
grid on, box on;
legend('Pulse', 'Gradient');
xlabel('Time [sec]');
