clear all; clc;

%% Initialise everything

cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R016-2012-10-08');

run(FindFile('*keys.m'));
%fname = ExpKeys.goodSWR{:}; wrong name?

fname_hippo = 'R016-2012-10-08-CSC02b.ncs';
fname_vs = 'R016-2012-10-08-CSC04d.ncs';

csc_hippo = myLoadCSC(fname_hippo);
csc_vs = myLoadCSC(fname_vs);

wSize = 1024; np = 1024; dsf = 4;

%% Part 1: White noise random signal

Fs = 2000; % in samples per second (Hz)

t0 = 0; t1 = 300; % start and end times

tvec = t0:1./Fs:t1; % construct time axis
y = rand(1,length(tvec)); % white noise distribution

% psd white noise

figure(1);

subplot(2,1,1);
plot(tvec,y,'k'); xlabel('Time(s)'); ylabel('Amplitude');
xlim([0 0.1]); ylim([-1 2]); title('White-Noise signal');

subplot(2,1,2);
[Pxx,F] = pwelch(y,hamming(wSize),wSize/2,np,Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 300]); title('White-Noise PSD');

%% Part 1 notes

% 1. Compute the PSD of “white noise”, i.e. a random signal where each sample is drawn independently from the open interval (0,1) with equal probability. Is it 1/f? How would you describe its shape? Hint: use the MATLAB function rand().

% From figure 1, the PSD of white noise not 1/f, but more like a jaggy flat 
% line with a spike at the begining. Flatness makes sence given the points 
% are random. Consequntly, the power for any one frequency is not very large. 
% The bigger the sample, the smoother (flatter) the distribution will be. 
% Not sure why it has a high number of low frequency components though.

%% Part 2: VS vs. HC

%psd hippocampus

header_hippo = getHeader(csc_hippo);
Fs_hippo = header_hippo.SamplingFrequency;

csc_hippo_pre = Restrict(csc_hippo,0,ExpKeys.TimeOnTrack(1)-10);
csc_hippo_preR = Range(csc_hippo_pre);
csc_hippo_preD = Data(csc_hippo_pre);

csc_hippo_preD = decimate(csc_hippo_preD,dsf);
csc_hippo_preR = downsample(csc_hippo_preR,dsf);
Fs_hippo = Fs_hippo./dsf;

figure(2);
[Pxx,F] = pwelch(csc_hippo_preD,hamming(wSize),wSize/2,np,Fs_hippo);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 300]); ylim([-40 40]);
title('Hippocampus LFP');

% psd ventral striatum

header_vs = getHeader(csc_vs);
Fs_vs = header_vs.SamplingFrequency;

csc_vs_pre = Restrict(csc_vs,0,ExpKeys.TimeOnTrack(1)-10);
csc_vs_preR = Range(csc_vs_pre);
csc_vs_preD = Data(csc_vs_pre);

csc_vs_preD = decimate(csc_vs_preD,dsf);
csc_vs_preR = downsample(csc_vs_preR,dsf);
Fs_vs = Fs_vs./dsf;

figure(3);
[Pxx,F] = pwelch(csc_vs_preD,hamming(wSize),wSize/2,np,Fs_vs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 300]); ylim([-40 40]);
title('Ventral Striatum LFP');

%% Part 2 notes

% 2. Compute the PSD of a LFP, simultaneously recorded with the signal above but now from the hippocampus. The ExpKeys specify the filename of a verified hippocampal LFP (in the GoodSWR field, for “good sharp wave-ripple complexes”, characteristic of the hippocampus). How does it compare to the ventral striatal LFP? 

% From figure 2, the hippocampus (HC) LFP seems to have a similiar delta (~4Hz)& theta (~7-9Hz) 
% oscillation, but not the other components [e.g beta (~12-20Hz), and gamma(~40-100Hz oscillations]
% that the ventral striatum (VS) LFP (figure 3) seems to have. Hence, it's more 
% "1/f"-like in shape (at least at < 200 Hz), which makes the 60Hz line noise more pronounced.
% Something similiar occurs at 180Hz between HC and VS, but not sure what is the cause of
% this spike.

%% Part 3: Testing wSize

wSize_list = [32 64 256 512 1024];
figure(4);

for s = 1:length(wSize_list)
    
    wSize = wSize_list(s);
    
    subplot(2,length(wSize_list),s)
    [Pxx,F] = pwelch(csc_hippo_preD,hamming(wSize),wSize/2,np,Fs_hippo);
    plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    xlim([0 300]); ylim([-40 40]); 
    title_name = sprintf('HC, wSize = %d', wSize);
    title(title_name);
    
    subplot(2,length(wSize_list),s + length(wSize_list))
    [Pxx,F] = pwelch(csc_vs_preD,hamming(wSize),wSize/2,np,Fs_vs);
    plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    xlim([0 300]); ylim([-40 40]);
    title_name = sprintf('VS, wSize = %d', wSize);
    title(title_name);
    
end

%% Part 3 notes

%  3. For both LFPs, explore the effects of the window size parameter of the Welch power spectrum, making sure to vary it across an order of magnitude. What do you notice? 

% From figure 4, The smaller the window size, the 'smoother' the psd is. However, this
% also means that sharp spikes (limitied and specific ranges of frequencies 
% where there is high power), cannot be accurately gauged. For example, the
% 'bump' due to the 60Hz line noise gets spread out and reachs a lower peak
% as wSize decreases.

%% Part 4: Downsample vs decimate

figure(5);

shift = 0.5; % displace the two graphs a bit, to compare

subplot(2,1,1); hold on;
% using previous decimated data
[Pxx,F] = pwelch(csc_hippo_preD,hamming(wSize),wSize/2,np,Fs_hippo);
plot(F,10*log10(Pxx),'k');
% naive downsampling
csc_hippo_preDo = Data(csc_hippo_pre);
csc_hippo_preDo = csc_hippo_preDo(1:dsf:length(csc_hippo_preDo));  
[Pxx,F] = pwelch(csc_hippo_preDo,hamming(wSize),wSize/2,np,Fs_hippo);
plot(F,10*log10(Pxx)-shift,'r'); 

xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 300]); ylim([-40 40]); 
title('HC, decimate vs downsample');
legend('Decimated','Naive downsampling');


subplot(2,1,2); hold on;
% using previous decimated data
[Pxx,F] = pwelch(csc_vs_preD,hamming(wSize),wSize/2,np,Fs_vs);
plot(F,10*log10(Pxx),'k');
 % naive downsampling
csc_vs_preDo = Data(csc_vs_pre);
csc_vs_preDo = csc_vs_preDo(1:dsf:length(csc_vs_preDo)); 
[Pxx,F] = pwelch(csc_vs_preDo,hamming(wSize),wSize/2,np,Fs_vs);
plot(F,10*log10(Pxx)-shift,'r'); 

xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 300]); ylim([-40 40]); 
title('HC, decimate vs downsample');
legend('Decimated','Naive downsampling');

%% Part 4 notes

% 4. Compare the PSD following the use of decimate() as in the above example, to a PSD obtained from downsampling without using decimate(). Are there any differences?

% From figure 5, at lower frequencies, there does not seem to be much of a difference between
% using decimate and downsampling without decimating. But at higher
% frequncies, the difference is more drastic.