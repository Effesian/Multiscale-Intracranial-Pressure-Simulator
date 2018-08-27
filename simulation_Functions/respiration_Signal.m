function [icpResp,foVec] = respiration_Signal(icpTrend,fs)

% Independent parameters
fo      = 10/60;   % Breaths per second
startT  = 1;
icpResp = zeros(size(icpTrend));
foVec   = fo*ones(size(icpTrend));
respCycle = 1;

% Simulate ICP
while(startT + round(fs/fo) < length(icpTrend))
    
  window = startT:(startT + round(fs/fo));
  amp    = 0.2*sqrt(abs(mean(icpTrend(window))));

  icpResp(window) = amp.*sin(2*pi*fo.*(1:length(window))/fs);
  icpResp(window) = icpResp(window) - mean(icpResp(window));
  foVec(window)   = ones(1,length(window))*fo;
 
  startT = window(end) + 1;
  respCycle = respCycle + 1;
end
end