clear all; clc;
basepath = pwd;

%% load the data and return to basepath
cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R042-2013-08-18');

fc = FindFiles('*.t');
S = LoadSpikes(fc);

cd(basepath)
t = [3200 5650];

%% for each of the two neurons, restrict the data to [3200 5650] (the time interval when the rat was running on the track)
 
cell1_id = 5; cell2_id = 42;
 
s1 = Restrict(S{cell1_id},t(1),t(2)); % restrict to on-track times only
s2 = Restrict(S{cell2_id},t(1),t(2));

spk_t1 = Data(s1);
spk_t2 = Data(s2);

%% compute the spike density function for each, making sure that your tvec runs from 3200 to 5650 also, and that you have a 50ms SD for the Gaussian convolution kernel
 
binsize = 0.001; % select a small bin size for good time resolution
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
 
spk_count_1 = histc(spk_t1,tbin_edges);
spk_count_1 = spk_count_1(1:end-1);

spk_count_2 = histc(spk_t2,tbin_edges);
spk_count_2 = spk_count_2(1:end-1);

binsize = 0.001; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.05./binsize; % 0.05 seconds (50ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize

gau_sdf_s1 = conv2(spk_count_1,gk,'same'); % convolve with gaussian window
gau_sdf_s2 = conv2(spk_count_2,gk,'same'); % convolve with gaussian window

% figure(); hold on;
% plot(tbin_centers,gau_sdf_s1,'g');
% plot(tbin_centers,gau_sdf_s2,'r');
% xlim([t(1) t(2)])

%% to use these SDFs to generate Poisson spike trains, convert the firing rates given by the SDF to a probability of emitting a spike in a given bin. (As you did above for a 0.47 Hz constant firing rate.)
 
pspike_1 = gau_sdf_s1/1000; % create variable expectation of a spike per bin
pspike_2 = gau_sdf_s2/1000;

%% generate Poisson spike trains, making sure to use the same tvec
 
dt = 0.001;
tvec = t(1):dt:t(2);
tvec = tvec(1:end-1)';
 
%rng default; % reset random number generator to reproducible state

spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx_1 = find(spk_poiss < pspike_1); % index of bins with spike

spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
spk_poiss_idx_2 = find(spk_poiss < pspike_2); % index of bins with spike

spk_poiss_t_1 = tvec(spk_poiss_idx_1)'; % use idxs to get corresponding spike time
spk_poiss_t_2 = tvec(spk_poiss_idx_2)'; % use idxs to get corresponding spike time

% figure(); hold on;
% line([spk_poiss_t_1' spk_poiss_t_1'],[-1 -0.5],'Color',[1 0 0]); % note, plots all spikes in one command
% line([spk_poiss_t_2' spk_poiss_t_2'],[-0.5 0],'Color',[0 1 0]); % note, plots all spikes in one command

%% convert Poisson spike trains to ts objects and compute the ccf

ts1 = ts(spk_poiss_t_1);
ts2 = ts(spk_poiss_t_2);

[xcorr_g,xbin_g] = ccf(ts1,ts2,0.01,1);
[xcorr_o,xbin_o] = ccf(s1,s2,0.01,1);

figure(); set(gcf,'Color',[1,1,1]); hold on;
plot(xbin_o,xcorr_o,'-r'); % original
plot(xbin_g,xcorr_g,'-k','linewidth',2.0); % generated

set(gca,'FontSize',10); xlabel('lag (s)'); ylabel('xcorr');
title(sprintf('%d-%d',cell1_id,cell2_id));

legend('Original','Generated','Location', 'NorthEastOutside')

% comparisons
% Generated spike train seems less 'sinusoidal'? Reducing is theta
% sequences?