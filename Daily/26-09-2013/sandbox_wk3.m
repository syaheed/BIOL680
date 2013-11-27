BaseFolder = cd ('C:/Users/Syaheed/Dropbox/Waterloo/Data analysis for Neuroscience/Daily/18-09-2013');

%%

fs1 = 2000;
tvec = 0:1/fs1:4; % construct time axis, sampled at 2kHz
 
freq1 = 100;
y = sin(2*pi*freq1*tvec); % 100Hz signal
 
ax1 = subplot(211);
stem(tvec,y); title('original');

%%

subsample_factor = 10;
 
tvec2 = tvec(subsample_factor:subsample_factor:end); % take every 4th sample
y2 = y(subsample_factor:subsample_factor:end);
 
ax2 = subplot(212);
stem(tvec2,y2,'r'); title('subsampled');
xlabel('time (s)');

%%

xl = [1 1.04];
set(ax1,'XLim',xl); set(ax2,'XLim',xl);

%%

hold on;
 
y_interp = interp1(tvec2,y2,tvec,'linear');
p1 = plot(tvec,y_interp,'b');
 
y_interp2 = interp1(tvec2,y2,tvec,'spline');
p2 = plot(tvec,y_interp2,'g');
 
legend([p1 p2],{'linear','spline'},'Location','Northeast'); legend boxoff

%% 

sin(2*pi) == 0 % ...right? The answer might surprise you.
fprintf('Welcome to numerical computing!\n');

%%

figure(2)
subsample_factor = 10;
 
tvec2 = tvec(2:subsample_factor:end); % best case scenario -- can detect 100Hz signal
y2 = y(2:subsample_factor:end);
 
subplot(212)
stem(tvec2,y2,'r');
set(gca,'XLim',xl);

%% 2kHz Fs, 100Hz signal with 450Hz signal superimposed

figure(3)
 
fs1 = 2000;
tvec = 0:1/fs1:4;
 
freq1 = 100;
freq2 = 450; % note, above Nyquist frequency for our target subsampled Fs
 
y = sin(2*pi*freq1*tvec) + 0.5.*sin(2*pi*freq2*tvec);
 
subplot(211)
stem(tvec,y)
set(gca,'XLim',xl);

%% ss -- we don't care about the 450Hz signal, but...
subsample_factor = 4;
 
tvec2 = tvec(1:subsample_factor:end);
y2 = y(1:subsample_factor:end);
 
subplot(212)
stem(tvec2,y2,'r');
set(gca,'XLim',xl);

%% 
figure()
y2 = decimate(y,subsample_factor); 
stem(decimate(tvec,4),y2,'r');
set(gca,'XLim',xl);

%%
cd('C:/Users/Syaheed/Dropbox/Waterloo/Data analysis for Neuroscience/Data/R016-2012-10-08');

fname = 'R016-2012-10-08-CSC03b.Ncs';
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

%%
cell_split = strsplit(Header{strncmp(Header,'-ADMax',6)});
ADBitVolts = str2double(cell_split{2});

cell_split_2 = strsplit(Header{strncmp(Header,'-InputR',7)});
InputRange = str2double(cell_split_2{2});

cell_split_3 = strsplit(Header{strncmp(Header,'-Sampling',9)});
SamplingFrequency = str2double(cell_split_3{2});

Samples_conv = Samples*1500*2/(ADBitVolts*2+2);
nSamples = size(Samples,1)*size(Samples,2);
Timestamps_conv = Timestamps * 10^-6;

Samp_time = nSamples/SamplingFrequency;
Time_time = Timestamps_conv(end)-Timestamps_conv(1);

Time_diff = Time_time - Samp_time;
plot(diff(Timestamps))

%%

run(FindFile('*keys.m'))
val_start = ExpKeys.TimeOnTrack(1);
val_end = ExpKeys.TimeOffTrack(1);

a = Timestamps_conv >= val_start & Timestamps_conv <= val_end;
l_a = length(a(a == 1));

time_start_ind = find(a==1,1);
time_end_ind = time_start_ind + l_a;
TimestampsValue = Timestamps_conv(Timestamps_conv >= val_start & Timestamps_conv <= val_end);

SamplesValue = Samples_conv(:,time_start_ind:time_end_ind);

plot(diff(TimestampsValue))

%%
format long
a = 512.*(1/2000) == mode(diff(TimestampsValue));

actual_sampling_freq = 1/(mode(diff(TimestampsValue))/512);

%%

plot(diff(TimestampsValue))
hold on;
plot(NumberOfValidSamples == 512,'r')

%%

indices_crippled = find(NumberOfValidSamples ~= 512);
Samples(:,indices_crippled); % replaces with zeroes?

%%

which tsd %C:\Users\Syaheed\Documents\GitHub\vandermeerlab\util\MClust-3.5\Utilities\@tsd\tsd.m  % tsd constructor
which Range %C:\Users\Syaheed\Documents\GitHub\vandermeerlab\util\MClust-3.5\Utilities\@tsd\Range.m  % tsd method

%%
Xtsd = tsd(Timestamps,Samples(1,:));
plot(Xtsd);

%% create tsd

cd('C:/Users/Syaheed/Dropbox/Waterloo/Data analysis for Neuroscience/Data/R042-2013-08-18');
fname = 'R042-2013-08-18-CSC03a.ncs';
myLoadCSC(fname)
