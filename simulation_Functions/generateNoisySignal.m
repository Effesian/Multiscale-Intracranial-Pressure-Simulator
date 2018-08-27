function [noisyICP] = generateNoisySignal(icp,fs)
% Adds noise to the signal or distorts the signal
% 
% 1. Signal drop
% 2. Saturation due to catheter flushing (Cliffords function)
% 3. Blood clot on the catheter damping the signal 
% 4. High frequency noise (movement)
% 

% Common artifacts in ICP signal
% Coughing (temporary increase in ICP)
% Sensor drift
% Hydrostatic error
% Transport of patient


% Initialize noise process
noisyICP = icp;
len = length(noisyICP);

%% Opening of the drain low-pass filter
idx = 1*round(len/10); 
window = idx:idx+round(30*fs);
[b,a] = butter(3,0.3/fs);
noisyICP(window) = filtfilt(b,a,noisyICP(window));


%% Hard Saturation (signal loss) at 20% of signal length (10s duration)
idx =  2*round(len/10);
window = idx:idx+round(8*fs);
noisyICP(window) = 0;


%% Gradual Saturation at 30% of signal length (30s duration)
idx =  round(2.5*len/10);
window = idx:idx+round(15*fs);

eta = 0.05;
ICPmax = 25;
A = noisyICP(idx);
noisyICP(window) = tanh((pi*(window-idx))/fs*eta)*(ICPmax-A) + A;


%% Hard Saturation (10s duration)
idx = round(len/40);
window = idx:idx+round(10*fs);
noisyICP(window) = 50;

end

