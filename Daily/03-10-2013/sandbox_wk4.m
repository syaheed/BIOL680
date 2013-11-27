% Week 4

%% plot a simple sinusoid
Fs = 100; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1./Fs:t1; % construct time axis
 
f = 2; % frequency of sine to plot
y = sin(2*pi*f*tvec); % note sin() expects arguments in radians, not degrees (see sind())
 
stem(tvec,y);

%% change phase and amplitude

phi = pi/2;
 
figure;
y = sin(2*pi*f*tvec + phi); % a phase shift
stem(tvec,y);
hold on;
plot(tvec,cos(2*pi*f*tvec),'LineWidth',2) % notice, cosine is simply phase shifted sine
legend('Sin (Phase Shifted)', 'Cos');
hold off;
 
a = 2;
 
figure;
y = a.*sin(2*pi*f*tvec + phi); % amplitude change
stem(tvec,y);

%% frequency modulation (FM) signal

f2 = 10;
m = 2;
 
subplot(311)
s1 = sin(2*pi*f*tvec);
plot(tvec,s1); title('message');
 
subplot(312);
s2 = sin(2*pi*f2*tvec);
plot(tvec,s2); title('carrier');
 
subplot(313);
s3 = sin(2*pi*f2*tvec + m.*sin(2*pi*f*tvec - pi/2));
plot(tvec,s3); title('FM signal');

%% Harmonic series example

mag = [0.1 0 1.3 0.5]; % magnitudes for each term
pha = [-pi/6 0 pi 2*pi/3]; % phases for each term
f = 2; % base frequency
 
signal_out = zeros(size(tvec));
for ii = 1:numel(mag) % note, the book chapter uses i, not best practice!
 
    this_signal = mag(ii)*cos(2*pi*f*ii*tvec + pha(ii));
    plot(tvec,this_signal,'r:'); hold on;
    signal_out = signal_out + this_signal; % build the sum
 
end
figure;
plot(tvec,signal_out,'LineWidth',2);

%% Fourier

x = round(rand(1,8)*10); % generate a length 8 vector of integers between 0 and 10
xlen = length(x);
 
% get magnitudes and phases of Fourier series

X = fft(x);
Xmag = abs(X); % magnitudes, a_n
Xphase = angle(X); % phases, phi_n
 
n = 0:xlen-1;
t = 0:0.05:xlen-1; % a finer timescale to show the smooth signal later
 
for iH = xlen-1:-1:0 % reconstruct each harmonic
    s(iH+1,:) = Xmag(iH+1)*cos(2*pi*iH*n/xlen + Xphase(iH+1))/xlen;
    sm(iH+1,:) = Xmag(iH+1)*cos(2*pi*iH*t/xlen + Xphase(iH+1))/xlen;
    % detail: xlen appears here because the fundamental frequency used by fft() depends on this
end
 
ssum = sum(s);
smsum = sum(sm);
 
figure;
plot(n, x, 'go', t, smsum, 'b', n, ssum, 'r*');
legend({'original','sum - all','sum - points only'});

%% fft()

Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1; % start and end times
tvec = t0:1/Fs:t1-(1/Fs); % construct time axis; generate exactly 20 samples
 
f = 2; % signal frequency
y = sin(2*pi*f*tvec); % construct signal, a 2Hz sine wave sampled at 20Hz for 1s
 
yfft = fft(y,length(y));
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
stem(yfft_mag)

Npoints = length(y);
F = [-Npoints/2:Npoints/2-1]./Npoints; % construct frequency axis
 
yfft_mag = fftshift(yfft_mag); % align output, see note below
stem(F,yfft_mag);
 
xlabel('Frequency (Fs^{-1})');

%%

tvec = t0:1/Fs:t1;
f = 2; % signal frequency
y = sin(2*pi*f*tvec); % construct signal, a 2Hz sine wave sampled at 20Hz for 1s
 
yfft = fft(y,length(y));
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
stem(yfft_mag)

Npoints = length(y);
F = [-Npoints/2:Npoints/2-1]./Npoints; % construct frequency axis
 
yfft_mag = fftshift(yfft_mag); % align output, see note below
stem(F,yfft_mag);
 
xlabel('Frequency (Fs^{-1})');

%%

tvec = t0:1/Fs:t1;
nPoints = [length(tvec) 64 256 1024];
 
for iP = 1:length(nPoints) % repeat fft with different numbers of points
 
    nP = nPoints(iP);
    subplot(2,2,iP);
 
    y = sin(2*pi*f*tvec);
    yfft = fft(y,nP);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
 
    title(sprintf('%d point FFT',nP));
    xlabel('Frequency (Fs^{-1})');
 
end

%%

tvec = t0:1/Fs:t1-(1/Fs);
nRepeats = [1 2 4 8];
 
nP =  1024;
 
for iP = 1:length(nRepeats)
 
    subplot(2,2,iP);
 
    y = sin(2*pi*f*tvec);
    y = repmat(y,[1 nRepeats(iP)]); % repeat the signal a number of times
 
    yfft = fft(y,nP);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
    F = [-nP/2:nP/2-1]./nP;
    yfft_mag = fftshift(yfft_mag);
    plot(F,yfft_mag,'kx',F,yfft_mag,'k');
 
    title(sprintf('%d repeats',nRepeats(iP)));
    xlabel('Frequency (Fs^{-1})');
 
end

%%

nP = 25;
nPFFT = 1024;
 
windows = {'rectwin','triang','hamming','hanning','blackman'};
cols = 'rgbcmyk';
 
for iW = 1:length(windows)
 
    eval(cat(2,'wn = ',windows{iW},'(nP);')); % make sure you understand this
    wn = wn./sum(wn);
 
    subplot(211); % plot the window
    plot(wn,cols(iW),'LineWidth',2); hold on;
 
    subplot(212);
    yfft = fft(wn,nPFFT);
    yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
    F = [-nPFFT/2:nPFFT/2-1]./nPFFT;
    yfft_mag = fftshift(yfft_mag);
 
    h(iW) = plot(F,yfft_mag,cols(iW),'LineWidth',2); hold on;
 
end
 
xlabel('Frequency (Fs^{-1})');
legend(h,windows);

%% Robust spectral estimation methods

[Pxx,F] = periodogram(y,[],nP,Fs);
plot(F,Pxx); xlabel('Frequency (Hz)');

hold on;
[Pxx,F] = periodogram(y,hanning(length(y)),nP,Fs);
plot(F,Pxx,'r');

%% pwelch()

Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1;
f = 2;
nRepeats = 4;
 
tvec = t0:1/Fs:t1-(1/Fs);
 
nP =  1024;
y = sin(2*pi*f*tvec);
y = repmat(y,[1 nRepeats]);
 
[Pxx,F] = periodogram(y,rectwin(length(y)),nP,Fs);
plot(F,Pxx);
 
hold on;
wSize = 40;
[Pxx,F] = pwelch(y,rectwin(wSize),wSize/2,nP,Fs);
plot(F,Pxx,'r'); xlabel('Frequency (Hz)');


%%

Fs = 20; % in samples per second (Hz)
t0 = 0; t1 = 1; f = 2;
nP =  1024;
gaps = [5 10 15]; % idx of samples to be removed
 
tvec = t0:1/Fs:t1;%-(1/Fs);
y = sin(2*pi*f*tvec);
 
subplot(211)
plot(tvec,y,'k*'); hold on;
 
yfft = fft(y,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);
 
subplot(212);
plot(F,yfft_mag,'kx',F,yfft_mag,'k'); hold on;
 
xlabel('Frequency (Fs^{-1})');
 
% signal with gaps
y = sin(2*pi*f*tvec);
y2(gaps) = []; tvec(gaps) = []; % remove
 
subplot(211);
plot(tvec,y2,'bo'); hold on;
 
yfft = fft(y2,nP);
yfft_mag = abs(yfft); yfft_ph = angle(yfft);
 
F = [-nP/2:nP/2-1]./nP;
yfft_mag = fftshift(yfft_mag);
 
subplot(212);
plot(F,yfft_mag,'bx',F,yfft_mag,'b');


%% Application to real data

cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R016-2012-10-08');
csc = myLoadCSC('R016-2012-10-08-CSC04d.ncs');
run(FindFile('*keys.m'));

% restrict to prerecord, leaving some time (10s) before rat actually goes on track
csc_pre = Restrict(csc,0,ExpKeys.TimeOnTrack(1)-10);

csc_preR = Range(csc_pre);
csc_preD = Data(csc_pre);
 
% check if sampling is ok
plot(diff(csc_preR)); % only minimal differences
Fs = 1./mean(diff(csc_preR));

% downsample
dsf = 4;
csc_preD = decimate(csc_preD,dsf);
csc_preR = downsample(csc_preR,dsf);
Fs = Fs./dsf;

wSize = 1024; nP = 1024;
[Pxx,F] = periodogram(csc_preD,hamming(length(csc_preD)),length(csc_preD),Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);

[Pxx,F] = pwelch(csc_preD,hamming(wSize),wSize/2,nP,Fs);
plot(F,10*log10(Pxx),'k'); xlabel('Frequency (Hz)'); ylabel('Power (dB)');
xlim([0 150]);
