BaseFolder = cd ('C:/Users/Syaheed/Dropbox/Waterloo/Data analysis for Neuroscience/Daily/18-09-2013');

%% first, cd to where the data you just grabbed is located
DataFolder = cd ('../../Data/R042-2013-08-18');

%% load the data (note, may need to unzip position data first)
fc = FindFiles('*.t');
S = LoadSpikes(fc);
 
[csc,csc_info] = LoadCSC('R042-2013-08-18-CSC03a.ncs');

[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );

cd(BaseFolder) % Jump back
%% Verify that your csc object indeed has timestamps and data of the same length.

if length(Data(csc)) == length(Range(csc))
    disp('csc has timestamps and data of the same length')
else
    disp('csc does not have timestamps and data of the same length')
end

%% csc_short variable with data restricted to between 5950 and 6050 s
%R = Restrict(X, t0, t1)
t0 = 5950; t1 = 6050;
csc_short = Restrict(csc,t0,t1);

%% tsd X and Y 
% tsa = tsd(t,data) , assumes secs
Timestamps_conv = Timestamps.*(10^-6);

tsa_X = tsd(Timestamps_conv,X');
tsa_Y = tsd(Timestamps_conv,Y');

% Verify that it worked by creating a X_short variable over the same range as csc_short

X_short = Data(Restrict(tsa_X,t0,t1))';

%% Plot the LFP (time against voltage) for the segment between 5950 and 6950 s

f_hd1 = figure(1);
hold on; box off;
set(gcf,'Color',[0 0 0],'InvertHardcopy','off');
set(gca,'Color',[0 0 0], 'XColor',[1 1 1], 'YColor',[1 1 1], 'XLim', [5989 5990], 'FontSize', 24);

csc_mean = nanmean(Data(csc));
xr = [t0 t1];

csc_short_plot = plot(Range(csc_short),Data(csc_short));
mean_hdl = plot(xr,[csc_mean csc_mean],'Color', [1 0 0], 'LineWidth', 2);

print(gcf,'-dpng','-r300','R042-2013-08-18-LFPsnippet.png');

%% Interactive plot
set(gca,'YLim', [-1000 1000]); % just to make it constant
key_pressed_shift(f_hd1)

%% Varargin
test_fun(1,'b',0,'c',0)

%% test

clf

lfp_color = [0 0 1];
color_set = repmat([1 0 0], length(S),1);

fighandle = figure(1); hold on; box off;
set(gcf,'Color',[1 1 1],'InvertHardcopy','off');
set(gca,'Color',[1 1 1], 'XColor',[0 0 0], 'YColor',[0 0 0], 'FontSize', 12);

xlabel('Time'); ylabel('Cell Number');
ncells = length(S);

for n = 1:ncells
    n
    time = Data(S{n});
    num_spikes = length(time);
        for i = 1:num_spikes
            plot([time(i) time(i)],[n-0.4 n+0.4], 'Color', color_set(n,:))
        end
end

csc_short_plot = plot(Range(csc_short),((Data(csc_short).*0.001)+ n + 3), 'Color', lfp_color);
mean_hdl = plot(xr,[csc_mean.*0.001 csc_mean.*0.001]+ n + 3,'Color', [1 0 0], 'LineWidth', 2);

axeshandle = get(fighandle,'CurrentAxes');
set(axeshandle,'yticklabel',[]);
xlim([t0 t0+2]);

%% Use neuroplot
%neuroplot(S,csc_short);