clear all; clc;

basepath = pwd;

%% load the data
cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R042-2013-08-18');

fc = FindFiles('*.t');
S = LoadSpikes(fc);

cd(basepath)
%% plot
iC = 47;
t = [5801 5801.7];
 
spk_t = Data(Restrict(S{iC},t(1),t(2))); % get spike times
 
line([spk_t spk_t],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command

binsize = 0.001; % given in seconds
tbin_edges = t(1):binsize:t(2); % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)
 
spk_count = histc(spk_t,tbin_edges); % get spike counts for each bin
spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.

hold on;
h = bar(tbin_centers,spk_count);
set(h,'BarWidth',1,'EdgeColor','none','FaceColor',[1 0 0]); % reformat bar appearance
 
yt = get(gca,'YTick'); ytl = get(gca,'YTickLabel');
set(gca,'YLim',[-1.5 10],'YTick',yt(2:end),'YTickLabel',ytl(2:end,:)); xlabel('time (s)');

%% 

binsize = 0.001; % select a small bin size for good time resolution
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
 
spk_count = histc(spk_t,tbin_edges);
spk_count = spk_count(1:end-1);
 
rw_sdf = conv2(spk_count,rectwin(50),'same'); % convolve with rectangular window
plot(tbin_centers,rw_sdf,'b');
 
gau_sdf = conv2(spk_count,gausswin(50),'same'); % convolve with gaussian window
plot(tbin_centers,gau_sdf,'g');

%%

binsize = 0.001; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.02./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
plot(tbin_centers,gau_sdf,'g');

%%

iC = 47;
spk_t = Data(S{iC}); % spike times
isi = diff(spk_t); % interspike intervals
 
dt = 0.001; % in s, because spike times are in s
isi_edges = 0:dt:0.25; % bin edges for ISI histogram
isi_centers = isi_edges(1:end-1)+dt/2; % for plotting
 
isih = histc(isi,isi_edges);
 
bar(isi_centers,isih(1:end-1)); % remember to ignore the last bin of histc() output
set(gca,'FontSize',20,'XLim',[0 0.25]); xlabel('ISI (s)'); ylabel('count'); grid on;

isi_shift = [isi(2:end); 0];
plot(isi,isi_shift,'ok','MarkerFaceColor','k','Markersize',2)
grid on;
xlim([0 0.25]); ylim([0 0.25]);
xlabel('isi(n)'); ylabel('isi(n+1)')

%%

dt = 0.001;
t = [0 4500]; % time interval (length) of spike train to generate
tvec = t(1):dt:t(2);
 
pspike = 0.47/1000; % probability of generating a spike in bin
rng default; % reset random number generator to reproducible state
spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
spk_poiss_t = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
 
line([spk_poiss_t spk_poiss_t],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
axis([0 0.1 -1.5 5]);

%%

isi = diff(spk_poiss_t); % interspike intervals
isih = histc(isi,isi_edges);
 
bar(isi_centers,isih(1:end-1)); % remember to ignore the last bin of histc() output
set(gca,'FontSize',20,'XLim',[0 0.25]); xlabel('ISI (s)'); ylabel('count'); grid on;


isi_shift = [isi(2:end); 0];
plot(isi,isi_shift,'ok','MarkerFaceColor','k','Markersize',2)
grid on;
xlim([0 0.25]); ylim([0 0.25]);
xlabel('isi(n)'); ylabel('isi(n+1)')

%%

acorr = acf(ts(spk_poiss_t),0.01,1);
plot([-1:0.01:1]',acorr)

[acorr,xbin] = acf(S{47},0.01,1);
plot([-1:0.01:1]',acorr) % not flat

%%

cell1_id = 5; cell2_id = 42;
 
s1 = Restrict(S{cell1_id},3200,5650); % restrict to on-track times only
s2 = Restrict(S{cell2_id},3200,5650);
 
[xcorr,xbin] = ccf(s1,s2,0.01,1);
 
plot(xbin,xcorr);
set(gca,'FontSize',20); xlabel('lag (s)'); ylabel('xcorr');
title(sprintf('%d-%d',iX,iY));




