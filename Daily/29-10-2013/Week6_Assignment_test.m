clear all; clc;

%% load the data
cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R016-2012-10-03');
fname = 'R016-2012-10-03-CSC04a.Ncs';
 
csc = myLoadCSC(fname);

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

%% Basic eventLFPplot parameters

time_window = [-1 3]; % how many seconds before ans after event times

event_times = event.ts*10^-6;
events = event_times(:); % to select events

eventLFPplot(csc,events,'t_window',time_window);

%% Fancier eventLFPplot parameters

decimate_ratio = 4; % to decimate or not, and by how much : ratio of 1 (default) does not decimate
colour_list = rand(length(events),3); % because I want random colours
FiltRange = [30 70]; % frequencies to KEEP

eventLFPplot(csc,events,'t_window',time_window, 'decimate_ratio', decimate_ratio, 'colour_list', colour_list, 'filter', 1, 'FiltRange', FiltRange);

%% Question & Answer

% Do you think the average is a reasonable charaterization of what is happening in the LFP for these events?

% Probably not. Not all events show the same pattern of LFPs. Averageing
% might dampen up the power spectrum observed, possibly leading to missed
% effects when doing comparisons.
