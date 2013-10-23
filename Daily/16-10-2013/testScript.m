clear all; clc;
BaseFolder = pwd;

%% Load data
cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R042-2013-08-18');
csc = myLoadCSC('R042-2013-08-18-CSC03a.ncs');
fc = FindFiles('*.t');
S = LoadSpikes(fc);
cd(BaseFolder)

%% Decide restriction

tvec = Range(csc);
t1 = tvec(end);
t0 = t1-(15*60);
cscR = Restrict(csc,t0,t1); % only look at lst 15 minutes

%% detectSWR implementation

evt = detectSWR(cscR);

%% Test with neuroplot

neuroplot(S,cscR, 'evt', evt.t)

%% Notes

% Optional strategies: Apply the procedure to LFPs recorded in other
% regions, since things like chewing might show up in multiple recordings.
% Filter out those common events to get hippocampal-specific neural activities
% like ripples.

% Z-scoring, looks at waves only at high amplitudes, more in-line with
% spike information?