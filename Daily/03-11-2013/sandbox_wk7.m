clear all; clc

x = randn(100,4);  % uncorrelated data - 4 signals of 100 points each
x(:,4) = sum(x,2);   % 4th signal is now correlated with first 3
[r,p] = corrcoef(x);  % compute sample correlation and p-values; notice signal 4 is correlated with the others
imagesc(r) % plot the correlation matrix -- note symmetry, diagonal, and r values of ~0.5 for signal 4

%%

cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R016-2012-10-03');
fname = 'R016-2012-10-03-CSC04a.Ncs';
csc = myLoadCSC(fname);

cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Daily/03-11-2013');


cscR = Restrict(csc,2700,3300); % risk session only
Fs = 2000;
 
[S,F,T,P] = spectrogram(Data(cscR),hanning(512),384,1:0.25:200,Fs); % spectrogram
 
[r,p] = corrcoef(10*log10(P')); % correlation matrix (across frequencies) of spectrogram
 
% plot
imagesc(F,F,r); 
caxis([-0.1 0.5]); axis xy; colorbar; grid on;
set(gca,'XLim',[0 150],'YLim',[0 150],'FontSize',14,'XTick',0:10:150,'YTick',0:10:150);

%%

wnoise = randn(1000,1);
[r,p] = corrcoef(wnoise(1:end-1),wnoise(2:end)); % same signal, offset by one sample

[acf,lags] = xcorr(wnoise,100,'coeff'); % compute correlation coefficients for lags up to 100
plot(lags,acf,'LineWidth',2); grid on;
set(gca,'FontSize',18); xlabel('time lag (samples)'); ylabel('correlation ({\itr})');

%%

Fs = 500; dt = 1./Fs;
tvec = 0:dt:2-dt;
 
wnoise = randn(size(tvec));
sgnl = sin(2*pi*10*tvec) + wnoise; % 10 Hz sine wave plus Gaussian white noise
 
subplot(211);
plot(tvec,sgnl); grid on;
set(gca,'FontSize',18); title('signal');
 
subplot(212);
[acf,lags] = xcorr(sgnl,100,'coeff');
lags = lags.*dt; % convert samples to time
plot(lags,acf); grid on;
set(gca,'FontSize',18); xlabel('time lag (s)'); ylabel('correlation ({\itr})');
title('autocorrelation');

%%

% spectrogram with better time resolution (slow!)
Fs = 500;
d = decimate(Data(cscR),4);
[S,F,T,P] = spectrogram(d,hanning(125),115,1:200,Fs);

F_idx = find(F > 70 & F < 100); % high-gamma power
F_idx2 = find(F > 50 & F < 65); % low-gamma power
 
pwr = mean(P(F_idx,:)); % average across frequencies
pwr2 = mean(P(F_idx2,:));
 
[ac,lags] = xcorr(pwr2-mean(pwr2),50,'coeff'); % remember to subtract the mean!
lags = lags.*mean(diff(T));
figure;
plot(lags,ac)
set(gca,'FontSize',18); xlabel('time lag (s)'); ylabel('correlation ({\itr})');
title('low-gamma autocorrelation');

%%

F_idx = find(F > 70 & F < 100); % high-gamma power
F_idx2 = find(F > 50 & F < 65); % low-gamma power
 
pwr = mean(P(F_idx,:)); % average across frequencies
pwr2 = mean(P(F_idx2,:));
 
[ac,lags] = xcorr(pwr-mean(pwr),pwr2-mean(pwr2),50,'coeff'); % note, two different input signals now
lags = lags.*mean(diff(T));
figure;
plot(lags,ac)
set(gca,'FontSize',18); xlabel('time lag (s)'); ylabel('correlation ({\itr})');
title('xcorr (high-gamma relative to low-gamma at t = 0)');

%%

Fs = 500; dt = 1./Fs;
tvec = 0:dt:1-dt;
 
wnoise = randn(size(tvec));
sgnl = sin(2*pi*4*tvec) + 0.1*wnoise;
 
h = hilbert(sgnl);
phi = angle(h);

[ax_h,h1,h2] = plotyy(tvec,sgnl,tvec,phi);
set(ax_h(2),'YColor',[1 0 0],'YLim',[-pi pi],'YTick',-pi:pi/2:pi,'YTickLabel',{'-p','-p/2','0','p/2','p'}, ...
   'fontname','symbol','XTick',[]);
set(h2,'Color',[1 0 0]);
set(get(ax_h(2),'Ylabel'),'String','phase (radians)') 
box off;

%%

Fs = 500; dt = 1./Fs;
tvec = 0:dt:1-dt;
 
f1 = 8;
s1 = sin(2*pi*f1*tvec)+0.05*randn(size(tvec));
s1_phi = angle(hilbert(s1)); % phase of slow signal
dphi = exp(-abs(s1_phi)/1.5); % envelope with peaks at phase 0, based on slow phase
 
f2 = 80;
s2 = 0.3.*sin(2*pi*f2*tvec).*dphi; % fast oscillation multiplied by envelope
 
s = s1+s2;
plot(tvec,s); xlim([0 0.5]);

%%

phi_edges = -pi:pi/8:pi;
pow_bin = averageXbyYbin(dphi,s1_phi,phi_edges);
 
phi_centers = phi_edges(1:end-1)+pi/16; % convert edges to centers
plot(phi_centers,pow_bin);

