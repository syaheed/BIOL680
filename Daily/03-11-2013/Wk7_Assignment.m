clear all; clc

%% load and restrict the data

cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Data/R016-2012-10-03');
fname = 'R016-2012-10-03-CSC04a.Ncs';
csc = myLoadCSC(fname);

cd('C:/Users/Syaheed/Documents/GitHub/BIOL680/Daily/03-11-2013');

cscR = Restrict(csc,2700,3300); % risk session only
Fs_orig = 2000;
Fs = 500;
d = decimate(Data(cscR),round(Fs_orig/Fs));

%% define frequency ranges of interest

Range1 = [3 4]; % delta
Range2 = [45 65]; % gamma

%% design filters for frequency ranges

Wp1 = Range1 * 2 / Fs; 
Ws1 = [Range1(1)-2 Range1(2)+2] * 2 / Fs;
[N,Wn] = cheb1ord( Wp1, Ws1, 3, 20); 
[b_c1,a_c1] = cheby1(N,0.5,Wn);

Wp2 = Range2 * 2 / Fs; 
Ws2 = [Range2(1)-2 Range2(2)+2] * 2 / Fs;
[N,Wn] = cheb1ord( Wp2, Ws2, 3, 20); 
[b_c2,a_c2] = cheby1(N,0.5,Wn);

%% filter the data (remember to use filtfilt!)

y1 = filtfilt(b_c1,a_c1,d);
y2 = filtfilt(b_c2,a_c2,d);

%% extract delta phase and low gamma power

s1_phi = angle(hilbert(y1)); % phase of delta sub-signal
s2_amp = abs(hilbert(y2));% fast oscillation amplitude

%% use averageXbyYbin to plot relationship (ideally with standard deviations)

phi_edges = -pi:pi/8:pi;
[pow_bin pow_sd pow_count] = averageXbyYbin(s2_amp,s1_phi,phi_edges);
 
phi_centers = phi_edges(1:end-1)+pi/16; % convert edges to centers
errorbar(phi_centers,pow_bin,(pow_sd./sqrt(pow_count)));

xlabel('Phase of Delta sub-signal (radians)')
ylabel('Average amplitude of Gamma sub-signal') 