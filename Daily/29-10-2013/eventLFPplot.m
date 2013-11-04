function eventLFPplot(csc,event_times,varargin)

    % function eventLFPplot(csc,event_times,varargin)
    %
    % INPUTS
    %
    % csc: [1 x 1] mytsd, LFP signal to be plotted
    % event_times: [nEvents x 1] double with event times to align LFP on
    %
    % varargins (with defaults):
    %
    % t_window: [2 x 1] double indicating time window to use, e.g. [-1 3] for 1 second before to 3 seconds after event times

    t_window = [-1 3];
    decimate_ratio = 1;
    colour_list = repmat([0 0 1],length(event_times),1);
    
    filter = 0;
    FiltRange = [3 100];
    
    extract_varargin;
    
    if size(colour_list,1) ~= length(event_times)
        warning('colour_list not proportional, setting to default ...');
        colour_list = repmat([0 0 1],length(event_times),1);
    end
    
    head = getHeader(csc);
    Fs = head.SamplingFrequency;
    fighandle = figure; hold on;
    
    csc_data = Data(csc);
    csc_range = Range(csc);
    
    if decimate_ratio ~= 1
            csc_data = decimate(csc_data,decimate_ratio);
            csc_range = downsample(csc_range,decimate_ratio);
            Fs = Fs/decimate_ratio;
    end
    
    if filter == 1
        Wp = FiltRange * 2 / Fs;
        Ws = [FiltRange(1)-2 FiltRange(2)+2] * 2 / Fs;
        [N,Wn] = cheb1ord( Wp, Ws, 3, 20);
        [b_c1,a_c1] = cheby1(N,0.5,Wn);
        csc_data = filtfilt(b_c1,a_c1,csc_data);
    end

    for nevent = 1:length(event_times)
        
        % *extract the corresponding piece of LFP
        idx = find(csc_range >= event_times(nevent)+ t_window(1) & csc_range <= event_times(nevent)+ t_window(2));
        x = csc_range(idx);
        y = csc_data(idx);
        
        % *replace the original time axis with a new one based on the time window asked for
        x = x - x(1) + t_window(1);
               
        %y-offset and rescaling
        y = (y*0.4) + 300*nevent;
        
        if nevent == 1
           ymin = min(y) -1 ; % just get the mininum y point for plotting purposes
           disp(ymin);
        end
        if nevent == length(event_times)
           ymax = max(y) + 1; % just get the maximum y point for plotting purposes
           disp(ymax);
        end
        
        % *plot the current piece of LFP
        plot(x,y, 'Color', colour_list(nevent,:));
        
    end
    
        % add a vertical line to the plot indicating time zero
        axeshandle = get(fighandle,'CurrentAxes');
        set(axeshandle,'yticklabel',[]);
        set(axeshandle,'ytick',[])
        
        plot([0 0],[ymin ymax],'-k')
        
        xlabel('Time-Window (Seconds from Event Times of Interest)');
        ylim([ymin max(y)]);
        title('Event-Based LFPs')

end