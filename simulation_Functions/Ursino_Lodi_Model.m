function [t,Pic,Ca,Cic] = Ursino_Lodi_Model(Tmax,fs)


%% Model parameter
%Tmax    = Tmax/60; % Time in seconds
h       = 1/(60*fs);  % step size
t       = (0:h:Tmax);
Nsteps  = numel(t);

Ro = 526.3; % CSF outflow resistance [mmHg*s/ml]
Ro = 12*Ro;

Rpv = 1.24; % Proximal venous resistance [mmHg*s/ml]
Rf  = 2.38*1e3; % CSF formation resistance [mmHg*s/ml]

DeltaCa1 = 0.75; % amplitude of sigmoidal curve [ml/mmHg]
DeltaCa2 = 0.075; % amplitude of sigmoidal curve [ml/mmHg]
Can      = 0.15; % basal arterial compliance [ml/mmHg]

kE  = 0.11; % elastance coefficient [1/ml]
kE  = 2.1*kE;

kR  = 4.91*1e4; % resistance coefficient [mmHg^3*s/ml]
tau = 20; % time constance [s]
qn  = 12.5; % ml/s
G   = 1.5; % Gain [ml/mmHg *100% CBF change^{-1}]
%G = 2*G;


%% Initial input quantities

mABP   = 100;
Pa     = simulateABP(mABP,Nsteps,fs);

Pic(1) = 9.5; % Intracranial pressure [mmHg]
Pc(1)  = 25; % Dural sinus pressure [mmHg]
Pvs    = 6; % mmHg
Ca(1)  = 0.15; % Arterial compliance [ml/mmHg]
Cic(1) = 1/(kE*Pic(1));

%% Input values

PaDot = 0;
Ii = 0;


%% Solve ODE via RK4 integration method

for k = 1:Nsteps

  % Compute filling volume arterial-arteriolar cerebrovascular bed
  Va(k) = Ca(k)*(Pa(k) - Pic(k));

  % Autoregulated resistance
  Ra = (kR*Can^2)/(Va(k)^2); % Hagen-Poiseuille

  % Compute capillary pressure (neglecting CSF production)
  Pc(k) = (Pa(k)*Rpv + Pic(k)*Ra)/(Rpv + Ra);

  % Compute CBF
  q(k) = (Pa(k)-Pc(k))/Ra;

  % Compute sigmoidal function for autoregulatory gain
  x(k) = (q(k)-qn)/qn;

  if(x(k)<0)
    DeltaCa = DeltaCa1;    
  else
    DeltaCa = DeltaCa2;        
  end

  kSigma = DeltaCa/4;

  sigmoidNumerator = (Can + DeltaCa/2) + (Can - DeltaCa/2)*exp(G*x(k)/kSigma);
  sigmoidDenominator = 1 + exp(G*x(k)/kSigma);
  sigmoid(k) = sigmoidNumerator/sigmoidDenominator;

  % Rate of change arterial-arteriolar compliance
  CaDot = @(t,y) (1/tau)*(-y+ sigmoid(k));
  Ca(k+1) = rk4_Wadehn(k,Ca(k),CaDot,h);

  % Change in ICP
  PicDot = @(t,y)(kE*y/(1+Ca(k)*kE*y))*(Ca(k)*PaDot + CaDot(t,Ca(k)) *(Pa(k)-y) + (Pc(k)-y)/Rf - (y-Pvs)/Ro + Ii);
  Pic(k+1) = rk4_Wadehn(k,Pic(k),PicDot,h);
  Cic(k+1) = 1/(kE*Pic(k+1));  % Ursino Eq(2)
end

Pic = Pic(1:end-1);
Ca  = Ca(1:end-1);
Cic = Cic(1:end-1);

%% Visualize
% ax(1) = subplot(4,1,1);
% plot(t,Pa,'LineWidth',2);
% xlabel('Time [s]');
% ylabel('ABP [mmHg]');
% 
% ax(2) = subplot(4,1,2);
% plot(t,Pic(1:end-1),'LineWidth',2);
% xlabel('Time [s]');
% ylabel('ICP [mmHg]');
% 
% ax(3) = subplot(4,1,3);
% plot(t,q,'LineWidth',2);
% xlabel('Time [s]');
% ylabel('CBFV [ml/s]');
% 
% ax(4) = subplot(4,1,4);
% plot(t,Va,'LineWidth',2);
% xlabel('Time [s]');
% ylabel('Arteriolar-Volume [ml]');
% 
% linkaxes(ax,'x');
% 
% keyboard;
end
