% This script is an example demonstrating the use of SBE for the SIRFLA off-resonance effect
% simulation, which was originally performed in Feldman et al 2013.

clear; 

% prepare sample 
fluid            = struct;     
fluid.J0         = 10;          % spectral density [Hz]
fluid.J1         = 11; 
fluid.J2         = 10;
fluid.omegaQHz   = 0;           % residual quadrupolar interaction [Hz]
fluid.deltaHz    = 0;           % local field off resonance

cartilage            = struct;     
cartilage.J0         = 1225;    
cartilage.J1         = 25; 
cartilage.J2         = 24;
cartilage.omegaQHz   = 0;       
cartilage.deltaHz    = 0;     

% prepare pulse seq
sirfla             = struct;
t1                 = 0;       % start time of the inversion pulse
tau                = 10e-3;   % [s]
tSp1               = 30.5e-3; % start time of Spoiler Gradient 1
tSp2               = 93e-3;   % start time of Spoiler Gradient 2  
tSpDur             = 1e-3;    % Spoiler gradients' duration
t2                 = 31.5e-3; 
tau2               = 0.25e-3;
TR                 = 94e-3;
GradAmp            = 100;     % [Hz] local residual Zeeman field produced by the gradient
sirfla.pulse       = @(t) (t >= t1) .* (t < t1+tau) .* (deg2rad(180)/tau)...
                       + (t >= t2) .* (t < t2+tau2) .* (deg2rad(90)/tau2);
sirfla.trPhase     = @(t) 0*t;
                  
sim.options        = odeset('RelTol', 1e-3, 'AbsTol', 1e-4);
FluidMxy           = [];
CartMxy            = [];

for j = -54:1:54
    fluid.deltaHz         = j; 
    cartilage.deltaHz     = j; 
    FluidIsochromats      = zeros(1,10);
    CartIsochromats       = zeros(1,10);
    for i = -5:4  % 10 isochromats
        sirfla.grad       = @(t)  (t >= tSp1) .* (t < tSp1+tSpDur) .* 2*pi*GradAmp * (i/5)...
                                 +(t >= tSp2) .* (t < tSp2+tSpDur) .* 2*pi*GradAmp * (i/5);
        sim.T0            = [1; zeros(14,1)];
        sim.time          = [t1, tSp1, t2, tSp2, TR];
        for rep = 1:5
            FluidISTO     = SBE(fluid, sirfla, sim);
            sim.T0        = (FluidISTO.T(:,end));
        end
        sim.time          = [t1, tSp1, t2, t2+tau2];
        FluidISTO         = SBE(fluid, sirfla, sim);
        FluidIsochromats(i+6) = FluidISTO.T(2,end)+FluidISTO.T(3,end);
        
        sim.T0            = [1; zeros(14,1)];
        for rep = 1:5
            CartISTO      = SBE(cartilage, sirfla, sim);
            sim.T0        = (CartISTO.T(:,end));
        end
        sim.time          = [t1, tSp1, t2, t2+tau2];
        CartISTO          = SBE(cartilage, sirfla, sim);
        CartIsochromats(i+6) = CartISTO.T(2,end)+CartISTO.T(3,end);
    end
    FluidMxy = cat(1, FluidMxy, abs(mean(FluidIsochromats)));
    CartMxy  = cat(1, CartMxy,  abs(mean(CartIsochromats)));
end


 
figure, plot(-54:54,FluidMxy); xlim([-54,54]);ylim([0,1]);grid on;
hold on, plot(-54:54,CartMxy);
xlabel('Off-resonance frequency (Hz)');
ylabel('Transverse Magnetisation (a.u.)');

