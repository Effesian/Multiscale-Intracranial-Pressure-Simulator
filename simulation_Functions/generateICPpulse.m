function icpBeats = generateICPpulse(t_Vector,icpMean,fs,foVec)

% Initialization
startT   = 1;
icpBeats = zeros(size(icpMean));

% Arrhythmias
with_arrhythmia = false;

% Interbeat intervals
IBI = fs;
IBI_sig = ones(size(icpMean))*fs; % Baseline heart rate 60 BPM
IBI_sig = IBI_sig + 5*sin(2*pi*foVec.*t_Vector + pi); % RSA phase shift of pi

pulseCounter = 0;
  
while(startT + IBI <= length(t_Vector))
        
  % Interbeat interval (min 50 samples and max 200 samples)
  IBI = max(200,min(50,round(mean(IBI_sig(startT:startT + 50)))));
  
  % Temporal window of pulse  
  windowT = startT:startT + IBI - 1;  
      
  % Pulse amplitude
  
  % Use for linear relation (crude approximation)
  %amp = 0.5*mean(icpMean(windowT)) + 1; 
  
  % Exponential-pressure volume relation
  amp = 30*exp(mean(icpMean(windowT))/80)/exp(1); 
 
  
  %--- Generate Gamma pdf based beat ------------------------------------->
  % The position of P1-P2-P3 (and thus a1,b1 and a2,b2 and a3,b3
  % are a function of the meanICP and compliance from the Ursoino model
  x = 0:0.1:5; % Should be long enough for pdf to have decayed

  % First Gamma 
  aOne = 5;
  bOne = 0.2;

  % Second Gamma 
  aTwo     = 6;
  bTwo     = 0.2;
  shiftTwo = 1;
  wTwo   = (1.2-0.75)/(40)*(mean(icpMean(windowT))-10) + 0.75;


  % Third Gamma 
  aThree   = 8;
  bThree   = 0.3;
  shiftThree = 2;
  wThree = (1.1-0.75)/(40)*(mean(icpMean(windowT))-10) + 0.75;


  % Add a PVC every 10th beat ---------------------------------------------
  if(mod(pulseCounter,10)== 1 && with_arrhythmia)
    IBI = round(1.5*IBI);
    windowT = startT - round(IBI/4):startT + IBI - 1 - round(IBI/4); 
    amp = amp/2;
    wTwo = 0;
    x = 0:0.1:8;
  end
  
  pulseOneGamma   = gampdf(x,aOne,bOne);
  pulseTwoGamma   = wTwo*gampdf(x-shiftTwo,aTwo,bTwo);
  pulseThreeGamma = wThree*gampdf(x-shiftThree,aThree,bThree);
  pulseGamma      = pulseOneGamma + pulseTwoGamma + pulseThreeGamma;
  pulse           = amp*(pulseGamma/max(pulseGamma));

  % Scale beat to IBI  
  x     = ((1:length(x))/length(x))*IBI;
  pulse = interp1(x,pulse,1:IBI,'pchip');
  
  % Flatten at the beginning and end for concatenating beats smoothly
  if(~(mod(pulseCounter,10)== 1|| ~mod(pulseCounter,10)== 2) && with_arrhythmia)
    if(IBI < 50)
      winFnct          = 0.9:-0.1:0;
      pulse(1:10)      = pulse(1:10).*fliplr(winFnct);
      pulse(end-9:end) = pulse(end-9:end).*winFnct;
    else
      winFnct           = 0.9:-0.05:0;
      pulse(1:19)       = pulse(1:19).*fliplr(winFnct);
      pulse(end-18:end) = pulse(end-18:end).*winFnct;
    end
  end
  
  
  % Plug in pulse into ICP signal
  if(mod(pulseCounter,10)== 1 && with_arrhythmia &&  pulseCounter >1)
    previous_beat_end = icpBeats(windowT(1)-1);
    pulse = pulse - pulse(1) + previous_beat_end; % make PVC align with previous beat
    offsetPVC = -(0:IBI-1)*(pulse(end) - previous_beat_onset)/(IBI-1) ; % make PVC align with following beat
    
    icpBeats(windowT) = pulse + offsetPVC;  
  else
    icpBeats(windowT) = pulse - mean(pulse);
    previous_beat_onset = icpBeats(windowT(1));
  end
  
  % Go further in time
  startT = windowT(end);
  pulseCounter = pulseCounter + 1;
end


end



