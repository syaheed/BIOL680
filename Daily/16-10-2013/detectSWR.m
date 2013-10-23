function evt = detectSWR(csc,varargin)

ripple_band = [140 180]; % frequency band to use
threshold = 5; % number of SDs above mean to use for detection
tolerance = 2;
t0 = 6000;
t1 = 6002;
extract_varargin;

cscR = Restrict(csc,t0,t1);

x = Data(cscR);
tvec = Range(cscR);
info = getHeader(cscR);
Fs = info.SamplingFrequency;

Wp = ripple_band * 2 / Fs; % passband - between 50 and 100 Hz
Ws = [ripple_band(1)-tolerance ripple_band(2)+tolerance] * 2 / Fs; % stopband
[N,Wn] = cheb1ord( Wp, Ws, 3, 20); 
[b,a] = cheby1(N,0.5,Wn);

y = filtfilt(b,a,x);
SWR_power = y.^2;
SWR_power_filtered = medfilt1(SWR_power,101); % filter window is specified in samples, so this is ~50ms
plot(tvec,x,'b',tvec,SWR_power_filtered,'r');

env_m = nanmean(SWR_power_filtered);
env_s = nanstd(SWR_power_filtered);
z_SWR = (SWR_power_filtered-env_m)/env_s; 

evt.t = tvec(z_SWR > threshold);
evt.pwr = z_SWR(z_SWR > threshold);