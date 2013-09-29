% Script test modded LoadCSC function vs orig function

% samples that are 'zero' withing the area where the block's NumberOfValidSamples ~= 512 is deleted (by
% myLoadCSC, and highlighted by call to Invalid(mytsd), when
% 'handleMissing' is set to true, otherwise equivalent to original function

% seems to clean up the data, refer to the sample plot with the suggested
% xrange.

%% Initialise

% load week 1 data, change accordingly
cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R016-2012-10-08');
fname = 'R016-2012-10-08-CSC03b.ncs';

%% test original & modded LoadCSC methods

[csc_orig, csc_orig_info] = LoadCSC(fname); % original
csc_mod1 = myLoadCSC(fname, 'handleMissing', false); % my version
csc_mod2 = myLoadCSC(fname, 'handleMissing', true); % my version

%% Test Range, Data, and getHeader, and invalid (using NumberOfValidSamples) timestamps functions of new mytsd class

csc_data_mod = Data(csc_mod2);
csc_timestamp_mod = Range(csc_mod2);
csc_info_mod = getHeader(csc_mod2);
csc_invalid_mod = Invalid(csc_mod2);

%% Create comparison plot between original and my LoadCSC, also to test plot commands

disp('Plotting comparison graphs, hold on ...');

% set some interesting areas to test
y_plot_range = [-400 400];
x_plot_range = [1263.0 1263.5];
%x_plot_range = [1218.8 1219.2]; % cleans up a messy portion

figure(); % initialise new figure

subplot(3,1,1);
plot(csc_orig);  % original result
xlim(x_plot_range); ylim(y_plot_range);
title('original LoadCSC');
ylabel('LFP (mV)');

subplot(3,1,2); hold on; 
plot(csc_mod1); % myLoadCSC without deletion
xlim(x_plot_range); ylim(y_plot_range);
title('myLoadCSC (invalid samples not excluded)');
ylabel('LFP (mV)');

subplot(3,1,3); hold on; 
for n = 1:length(csc_invalid_mod) % to show that Invalid() function returns the correct 'area'
    plot([csc_invalid_mod(n) csc_invalid_mod(n)],[-csc_info_mod.InputRange csc_info_mod.InputRange], '-', 'color',[1 0.7 0.7],'linewidth', 0.001);
end
plot(csc_mod2); % myLoadCSC with deletion
xlim(x_plot_range); ylim(y_plot_range);
title('myLoadCSC (invalid samples excluded and highlighted)');
ylabel('LFP (mV)');
xlabel('Time (seconds)');
