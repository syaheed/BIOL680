clear all; clc;

%% load the data
cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R016-2012-10-03');
fname = 'R016-2012-10-03-CSC04a.Ncs';
 
csc = myLoadCSC(fname);

hdr = getHeader(csc);
Fs = hdr.SamplingFrequency;
 
%% restrict data
cscR = Restrict(csc,3282,3286); % if you don't have this, implement it (it's one line of code!)
plot(cscR); % note there are various gamma oscillations present, as well as a large negative-going transient
 
%% construct and plot the spectrogram

figure();

subplot(2,1,1);
[S,F,T,P] = spectrogram(Data(cscR),hanning(512),384,1:0.25:200,Fs);
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');

hold on;
tvec = Range(cscR);
data = Data(cscR);
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = tvec - tvec(1); % align LFP with spectrogram
data = rescale(data,-lfp_minmax,lfp_minmax); data = data+lfp_cent;
 
lfp_h = plot(tvec0,data,'k');

subplot(2,1,2);
[S,F,T,P] = spectrogram(Data(cscR),hanning(1024),384,1:0.25:200,Fs);
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');

hold on;
tvec = Range(cscR);
data = Data(cscR);
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = tvec - tvec(1); % align LFP with spectrogram
data = rescale(data,-lfp_minmax,lfp_minmax); data = data+lfp_cent;
 
lfp_h = plot(tvec0,data,'k');

%%

cscR = Restrict(csc,3300,3340); % section of data with a gap
 
[S,F,T,P] = spectrogram(Data(cscR),rectwin(256),128,1:200,Fs);
 
imagesc(T,F,10*log10(P)); 
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
 
hold on;
tvec = Range(cscR); data = Data(cscR);
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = tvec - tvec(1); % align LFP with spectrogram
data = rescale(data,-lfp_minmax,lfp_minmax); data = data+lfp_cent;
 
lfp_h = plot(tvec0,data,'k');
xlim([tvec0(1) tvec0(end)]);

%%

[S,F,T,P] = spectrogram(Data(cscR),hanning(512),384,1:0.25:200,Fs-100);
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');


%%

fn = FindFile('*Events.nev');
[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

%%

%% remember to cd to your data folder
fc = {'R016-2012-10-03-CSC04a.ncs'};
data = ft_read_neuralynx_interp(fc);

cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5;
cfg.trialdef.post = 5;
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk'
cfg.trialdef.cue = {'c1','c3','c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'}
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);

plot(data_trl.time{47},data_trl.trial{47})

%%

cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 10:1:100; % frequencies of interest
cfg.t_ftimwin    = ones(size(cfg.foi)).*0.5;  % window size: fixed at 0.5s
cfg.toi          = -2:0.025:5; % times of interest
 
TFR = ft_freqanalysis(cfg, data_trl);
 
figure
cfg = []; cfg.channel = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%%

figure
cfg = [];
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.channel      = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%%

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100;
cfg.keeptrials   = 'yes'; % need this for stats later
cfg.t_ftimwin    = 20./cfg.foi;  % 20 cycles per time window
 
cfg.toi          = -2:0.05:0; % pre-nosepoke baseline
TFR_pre = ft_freqanalysis(cfg, data_trl);
 
cfg.toi          = 0:0.05:2; % post-nosepoke
TFR_post = ft_freqanalysis(cfg, data_trl);
 
TFR_pre.time = TFR_post.time; % time axes should be identical for comparison

%% t-test
cfg = [];
cfg.channel     = 'R016-2012-10-03-CSC04a';
cfg.latency     = 'all';
cfg.trials      = 'all';
cfg.frequency   = 'all';
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
cfg.parameter   = 'powspctrm';
cfg.method      = 'stats';
cfg.statistic   = 'ttest2';
cfg.alpha       = 0.05;
 
nTrials1 = size(TFR_pre.powspctrm,1); nTrials2 = size(TFR_post.powspctrm,1);
 
cfg.design = cat(2,ones(1,nTrials1),2*ones(1,nTrials2)); % two conditions
cfg.ivar = 1; % dimension of design var which contains the independent variable (group)
 
stat = ft_freqstatistics(cfg,TFR_post,TFR_pre);
 
cfg.parameter = 'stat';
ft_singleplotTFR(cfg,stat); % plot the t-statistic

%%

fc = FindFiles('*.ncs'); % get filenames of all LFPs recorded
data = ft_read_neuralynx_interp(fc); % load them all -- this will take a while
 
% define layout
cfg = [];
cfg.layout = 'ordered'; cfg.channel = data.label;
layout = ft_prepare_layout(cfg, data);

%%

cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5;
cfg.trialdef.post = 5;
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk'
cfg.trialdef.cue = {'c1','c3','c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'}
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);
 
%%
cfg              = [];
cfg.output       = 'pow';
%cfg.channel      = 'R016-2012-10-03-CSC04a';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100;
cfg.keeptrials = 'yes'; % should be needed for subsequent statistics...
cfg.t_ftimwin    = 20./cfg.foi;  % 20 cycles per time window
cfg.toi          = -1:0.05:4;
 
TFR = ft_freqanalysis(cfg, data_trl);

figure
cfg = [];
cfg.baseline     = [-2 0];
cfg.baselinetype = 'relative';
cfg.layout = layout;
 
ft_multiplotTFR(cfg, TFR); % note this is now multiplot rather than singleplot