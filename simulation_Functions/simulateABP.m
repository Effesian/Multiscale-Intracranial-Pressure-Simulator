function abp = simulateABP(mABP,Nsteps,fs)
% Simulator for ABP traces
% Mean reverting process (Ornstein–Uhlenbeck process) + 
% Mayer waves

% Wikipedia: Mayer waves are cyclic changes or waves in arterial blood pressure
% brought about by oscillations in 
% baroreceptor and chemoreceptor reflex control systems

abp(1) = mABP;
theta = 0.1;
sig = 1;  

for k = 1:Nsteps
  abp(k+1) = abp(k)+ theta*(mABP-abp(k)) + sig*randn;    
end

abp = abp(1:end-1);
abp = movmean(abp,5*fs);
abp =  abp + 10*sin((1:Nsteps)*0.0001);
