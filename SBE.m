function ISTO = SBE(sample, seq, sim, version, outputform)
%SBE returns a struct containing ISTO integrated from the Spin-3/2 Bloch Equation.
%Input:
%   sample (struct) specifies the sample's 
%                       - sample.J0 [Hz] - spectral densities 
%                       - sample.J1 [Hz]
%                       - sample.J2 [Hz]
%                       - sample.omegaQHz [Hz] - residual quadrupolar interaction 
%                       - sample.deltaHz [Hz] - the sample's off resonance 
%   seq (struct) specifies 
%                       - seq.pulse [rad/s] - the pulse envelop
%                       - seq.trPhase [rad/s] - the transmit phase  
%                       - seq.grad [rad/s] - the gradient-driven off resonance. The residual Zeeman field in the differential equation is
%                                               seq.grad + 2*pi*sample.deltaHz .
%   sim (struct) specifies
%                       - sim.options (struct) - the ode options
%                       - sim.T0 - the initial ISTO
%                       - sim.time (n-by-1 array) - the simulated time interval. If the sequence contains short pulses,
%                       also specify their execution times to make sure the ode solver does not miss the pulses.
%                       - sim.teval (n-by-1 array, optional) - if specified, the ISTO will be evaluated at the time
%                       points teval.
%   version (string, optional)
%                       if version =='default'
%                           seq.pulse, seq.trPhase and seq.grad are function handles
%                       if version == 'ConstParam'
%                           seq.pulse, seq.trPhase and seq.grad are variables
%
%   outputform (string, optional)
%                       if outputform =='less'
%                           the output ISTO contains only ISTO.tsim (the
%                           simulated time steps) and ISTO.T (the full
%                           matrix of ISTO evolution).
%                       if outputform == 'full'
%                           the output ISTO contains also each ISTO components 
%                           extracted from ISTO.T.
%Output: 
%   ISTO (struct)
%                       refer to outputform.


%Authur:        Chengchuan Wu (chengchuanw@student.unimelb.edu.au)
%Institution:   The University of Melbourne
%Reference:     Chengchuan Wu, Yasmin Blunck, Leigh A. Johnston, 
%               The ‘Spin-3/2 Bloch Equation’: System matrix formalism of excitation, relaxation and off-resonance effects in biological tissue, 2022


if nargin < 4
    version = 'default';
    outputform = 'less';
elseif nargin < 5
    outputform = 'less';
end

switch version
    case 'default'
    % Differential Equations
    sodium          = @(J0, J1, J2, omegaQ, delta, omega1, phi, T, t)[(-2/5).*J1+(-8/5).*J2,omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),0,0,0,0,0,(-4/5).*J1+(4/5).*J2,0,0,0,0,0,0;
        (-1).*omega1(t).*sin(-phi(t)),(-3/5).*J0+(-1).*J1+(-2/5).*J2,delta(t)*1i,0,0,1i.*(3/5).^(1/2).*omegaQ,0,0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,0,0;
        (1i*(-1)).*omega1(t).*cos(phi(t)),delta(t)*1i,(-3/5).*J0+(-1).*J1+(-2/5).*J2,0,1i.*(3/5).^(1/2).*omegaQ,0,0,0,0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,0;
        0,0,0,(-2).*J1+(-2).*J2,3.^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*3.^(1/2).*omega1(t).*cos(phi(t)),0,0,0,0,0,0,0,0,0;
        0,0,1i.*(3/5).^(1/2).*omegaQ,(-1).*3.^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J0+(-1).*J1+(-2).*J2,delta(t)*1i,(1i*(-1)).*omega1(t).*cos(phi(t)),omega1(t).*sin(-phi(t)),0,0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,0;
        0,1i.*(3/5).^(1/2).*omegaQ,0,(1i*(-1)).*3.^(1/2).*omega1(t).*cos(phi(t)),delta(t)*1i,(-1).*J0+(-1).*J1+(-2).*J2,omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,0,0;
        0,0,0,0,(1i*(-1)).*omega1(t).*cos(phi(t)),(-1).*omega1(t).*sin(-phi(t)),(-1).*J0+(-2).*J1+(-1).*J2,2*delta(t)*1i,0,0,0,0,1i.*omegaQ,0,0;
        0,0,0,0,(-1).*omega1(t).*sin(-phi(t)),(1i*(-1)).*omega1(t).*cos(phi(t)),2*delta(t)*1i,(-1).*J0+(-2).*J1+(-1).*J2,0,0,0,1i.*omegaQ,0,0,0;
        (-4/5).*J1+(4/5).*J2,0,0,0,0,0,0,0,(-8/5).*J1+(-2/5).*J2,6.^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*6.^(1/2).*omega1(t).*cos(phi(t)),0,0,0,0;
        0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,1i.*(2/5).^(1/2).*omegaQ,0,0,(-1).*6.^(1/2).*omega1(t).*sin(-phi(t)),(-2/5).*J0+(-1).*J1+(-3/5).*J2,delta(t)*1i,(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),(5/2).^(1/2).*omega1(t).*sin(-phi(t)),0,0;
        0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,(1i*(-1)).*6.^(1/2).*omega1(t).*cos(phi(t)),delta(t)*1i,(-2/5).*J0+(-1).*J1+(-3/5).*J2,(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),0,0;
        0,0,0,0,0,0,0,1i.*omegaQ,0,(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),(-1).*(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J0+(-1).*J2,2*delta(t)*1i,(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),(3/2).^(1/2).*omega1(t).*sin(-phi(t));
        0,0,0,0,0,0,1i.*omegaQ,0,0,(-1).*(5/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(5/2).^(1/2).*omega1(t).*cos(phi(t)),2*delta(t)*1i,(-1).*J0+(-1).*J2,(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t));
        0,0,0,0,0,0,0,0,0,0,0,(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),(-1).*(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(-1).*J1+(-1).*J2,3*delta(t)*1i;
        0,0,0,0,0,0,0,0,0,0,0,(-1).*(3/2).^(1/2).*omega1(t).*sin(-phi(t)),(1i*(-1)).*(3/2).^(1/2).*omega1(t).*cos(phi(t)),3*delta(t)*1i,(-1).*J1+(-1).*J2] * T...
        + [2/5.*J1+8/5.*J2;0;0;0;0;0;0;0;4/5.*J1-4/5.*J2;0;0;0;0;0;0];
    
    case 'ConstParam'
    sodium          = @(J0, J1, J2, omegaQ, delta, omega1, phi, T, t)[(-2/5).*J1+(-8/5).*J2,omega1.*sin(-phi),(1i*(-1)).*omega1.*cos(phi),0,0,0,0,0,(-4/5).*J1+(4/5).*J2,0,0,0,0,0,0;
        (-1).*omega1.*sin(-phi),(-3/5).*J0+(-1).*J1+(-2/5).*J2,delta*1i,0,0,1i.*(3/5).^(1/2).*omegaQ,0,0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,0,0;
        (1i*(-1)).*omega1.*cos(phi),delta*1i,(-3/5).*J0+(-1).*J1+(-2/5).*J2,0,1i.*(3/5).^(1/2).*omegaQ,0,0,0,0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,0;
        0,0,0,(-2).*J1+(-2).*J2,3.^(1/2).*omega1.*sin(-phi),(1i*(-1)).*3.^(1/2).*omega1.*cos(phi),0,0,0,0,0,0,0,0,0;
        0,0,1i.*(3/5).^(1/2).*omegaQ,(-1).*3.^(1/2).*omega1.*sin(-phi),(-1).*J0+(-1).*J1+(-2).*J2,delta*1i,(1i*(-1)).*omega1.*cos(phi),omega1.*sin(-phi),0,0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,0;
        0,1i.*(3/5).^(1/2).*omegaQ,0,(1i*(-1)).*3.^(1/2).*omega1.*cos(phi),delta*1i,(-1).*J0+(-1).*J1+(-2).*J2,omega1.*sin(-phi),(1i*(-1)).*omega1.*cos(phi),0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,0,0;
        0,0,0,0,(1i*(-1)).*omega1.*cos(phi),(-1).*omega1.*sin(-phi),(-1).*J0+(-2).*J1+(-1).*J2,2*delta*1i,0,0,0,0,1i.*omegaQ,0,0;
        0,0,0,0,(-1).*omega1.*sin(-phi),(1i*(-1)).*omega1.*cos(phi),2*delta*1i,(-1).*J0+(-2).*J1+(-1).*J2,0,0,0,1i.*omegaQ,0,0,0;
        (-4/5).*J1+(4/5).*J2,0,0,0,0,0,0,0,(-8/5).*J1+(-2/5).*J2,6.^(1/2).*omega1.*sin(-phi),(1i*(-1)).*6.^(1/2).*omega1.*cos(phi),0,0,0,0;
        0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,0,0,1i.*(2/5).^(1/2).*omegaQ,0,0,(-1).*6.^(1/2).*omega1.*sin(-phi),(-2/5).*J0+(-1).*J1+(-3/5).*J2,delta*1i,(1i*(-1)).*(5/2).^(1/2).*omega1.*cos(phi),(5/2).^(1/2).*omega1.*sin(-phi),0,0;
        0,0,(-1/5).*6.^(1/2).*(J0+(-1).*J2),0,1i.*(2/5).^(1/2).*omegaQ,0,0,0,(1i*(-1)).*6.^(1/2).*omega1.*cos(phi),delta*1i,(-2/5).*J0+(-1).*J1+(-3/5).*J2,(5/2).^(1/2).*omega1.*sin(-phi),(1i*(-1)).*(5/2).^(1/2).*omega1.*cos(phi),0,0;
        0,0,0,0,0,0,0,1i.*omegaQ,0,(1i*(-1)).*(5/2).^(1/2).*omega1.*cos(phi),(-1).*(5/2).^(1/2).*omega1.*sin(-phi),(-1).*J0+(-1).*J2,2*delta*1i,(1i*(-1)).*(3/2).^(1/2).*omega1.*cos(phi),(3/2).^(1/2).*omega1.*sin(-phi);
        0,0,0,0,0,0,1i.*omegaQ,0,0,(-1).*(5/2).^(1/2).*omega1.*sin(-phi),(1i*(-1)).*(5/2).^(1/2).*omega1.*cos(phi),2*delta*1i,(-1).*J0+(-1).*J2,(3/2).^(1/2).*omega1.*sin(-phi),(1i*(-1)).*(3/2).^(1/2).*omega1.*cos(phi);
        0,0,0,0,0,0,0,0,0,0,0,(1i*(-1)).*(3/2).^(1/2).*omega1.*cos(phi),(-1).*(3/2).^(1/2).*omega1.*sin(-phi),(-1).*J1+(-1).*J2,3*delta*1i;
        0,0,0,0,0,0,0,0,0,0,0,(-1).*(3/2).^(1/2).*omega1.*sin(-phi),(1i*(-1)).*(3/2).^(1/2).*omega1.*cos(phi),3*delta*1i,(-1).*J1+(-1).*J2] * T...
        + [2/5.*J1+8/5.*J2;0;0;0;0;0;0;0;4/5.*J1-4/5.*J2;0;0;0;0;0;0];
end
    
    J0          = sample.J0;
    J1          = sample.J1;
    J2          = sample.J2;
    deltaHz     = sample.deltaHz;
    omegaQHz    = sample.omegaQHz;
    if strcmp(version, 'ConstParam')
        Delta   = seq.grad + 2*pi*deltaHz;
    else
        Delta       = @(t) seq.grad(t) + 2*pi*deltaHz;
    end
    pulse       = seq.pulse;
    trPhase     = seq.trPhase;
    tsim        = sim.time;
    T0          = sim.T0;
    options     = sim.options;
    
    for j = 1:numel(tsim)-1
        if j == 1
            sol         = ode45(@(t, T)sodium(J0, J1, J2, 2*pi*omegaQHz, Delta, pulse, trPhase, T, t),...
                            [tsim(j), tsim(j+1)], T0, options);
        else
            sol = odextend(sol,[],tsim(j+1));
        end
    end
    ISTO        = struct; 
    if isfield(sim,'teval')
        ISTO.tsim = sim.teval;
        ISTO.T  = deval(sol, sim.teval);
    else
        ISTO.tsim   = sol.x;
        ISTO.T      = sol.y;
    end
    if strcmp(outputform,'full')
        ISTO.T10    = ISTO.T(1,:);
        ISTO.T11a   = ISTO.T(2,:);
        ISTO.T11s   = ISTO.T(3,:); % pure imaginary
        ISTO.T20    = ISTO.T(4,:);
        ISTO.T21a   = ISTO.T(5,:);
        ISTO.T21s   = ISTO.T(6,:); % pure imaginary
        ISTO.T22a   = ISTO.T(7,:); % pure imaginary
        ISTO.T22s   = ISTO.T(8,:);
        ISTO.T30    = ISTO.T(9,:);
        ISTO.T31a   = ISTO.T(10,:);
        ISTO.T31s   = ISTO.T(11,:); % pure imaginary
        ISTO.T32a   = ISTO.T(12,:); % pure imaginary
        ISTO.T32s   = ISTO.T(13,:);
        ISTO.T33a   = ISTO.T(14,:);
        ISTO.T33s   = ISTO.T(15,:); % pure imaginary
        ISTO.Mxy    = ISTO.T(2,:)+ISTO.T(3,:);
        ISTO.Mx     = ISTO.T11a;
        ISTO.My     = ISTO.T11s / -1i;
        
        
    end

end
