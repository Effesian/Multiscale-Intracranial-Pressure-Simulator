%% The Multiscale ICP Simulator
% simulates synthetic ICP signals using the model presented in Wadehn, et al. at the CinC 2018 conference.
% The slow ICP trends are obtained by simulating the Ursino-Lodi model and the
% respiratory-induced changes in ICP and the ICP pulses are statistical
% models, i.e., do not follow a physiological model, but are merely aimed
% at producing waveforms with similar shape to real recordings.
 
clear all;
close all;
clc;

addpath(genpath('simulation_Functions'));

% Global simulation parameters
Tmax     = 1500; % Total simulation time [s]
fs       = 125; % General sampling rate (lower sampled signals upsampled to this value)
t_Vector = (1:1/fs:Tmax);
fsLodi   = 1; % Sampling rate/simulation steps (1 Hz) of Ursino-Lodi Model 


%% Ursino Lodi Model for long term ICP trend

% Simulate ICP trend signal
[tTrend,icpTrend,Ca,Cic] = Ursino_Lodi_Model(Tmax,fsLodi);

% Upsample to 125 Hz
icpTrend = interp1(tTrend,icpTrend,t_Vector);
CaTrend  = interp1(tTrend,Ca,t_Vector);
CicTrend = interp1(tTrend,Cic,t_Vector);


%% Simulate respiration signal
[icpResp,foVec] = respiration_Signal(icpTrend,fs);


%% Simulate single beats
icpBeat = generateICPpulse(t_Vector,icpTrend + icpResp,fs,foVec);


%% Combine synthesized ICP signals
icp = icpTrend + icpBeat + icpResp;


%% Noise process simulation
[noisyICP] = generateNoisySignal(icp,fs);
icp = noisyICP;


%% Post-process signal to have less jumps
icp = movmean(icp,25);


%% Visualize simulation results

% Plot ICP
plot(t_Vector,icp,'LineWidth',2); hold on;
xlabel('Time [s]');
ylabel('ICP [mmHg]');

keyboard;
