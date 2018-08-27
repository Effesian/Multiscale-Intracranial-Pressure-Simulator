function [signal] = lowpassFilterPhysiologicalSignal(signal,fs,varargin)
% LOWPASSFILTERPHYSIOLOGICALSIGNAL Returns low-pass filtered physiological signal.
% If no cutoff frequency is defined, a 4-th order Butterworth filter with a cutoff frequency of 12 Hz 
% is used to filter the signl.
% This cut off seems reasonable for waveforms like ABP, CBFV and ICP (not for ECG though).

if(nargin < 3)
  fc = 12; 
else
  fc = varargin{1};
end

fNyquist = fs/2;
order    = 4;

[b,a] = butter(order,fc/fNyquist);
signal = filtfilt(b,a,signal);
end

