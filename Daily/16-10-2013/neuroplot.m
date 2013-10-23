function neuroplot(spikes,csc,varargin)
    
    cscColor = repmat([0 0 1], length(csc),1);
    spikeColor = repmat([1 0 0], length(spikes),1);
    t0 = 5950; t1 = t0+2;
    interactiveMode = true;
    evt = [t0 t1];
    
    evtColor= repmat([0 0 0], length(evt),1);
    
    extract_varargin;
    nEvents = length(evt);
    
    if size(spikeColor,1) ~= length(spikes)
       warning('spikeColor lists not proportional, setting to default ...');
       spikeColor = repmat([1 0 0], length(spikes),1);
    end
    
    if size(cscColor,1) ~= length(csc)
       warning('cscColor lists not proportional, setting to default ...');
       cscColor = repmat([0 0 1], length(csc),1);
    end
 
    if size(evtColor,1) ~= nEvents
       warning('evtColor lists not proportional, setting to default ...');
       evtColor= repmat([0 0 0], nEvents,1);
    end
    
    fighandle = figure(); hold on; box off;
    set(gcf,'Color',[1 1 1]);
    set(gca,'Color',[1 1 1], 'XColor',[0 0 0], 'YColor',[0 0 0], 'FontSize', 12);

    xlabel('Time'); %ylabel('Cell Number');
    
    ncells = length(spikes);
    
    fprintf('\nProcessing spike data. May take awhile depending on number of cells to plot.\n');
    
    for cell_num = 1:ncells
        fprintf('Processing cell %d of %d\n',cell_num, ncells);
        time = Data(spikes{cell_num});
        num_spikes = length(time);
            for i = 1:num_spikes
                plot([time(i) time(i)],[cell_num-0.4 cell_num+0.4], 'Color', spikeColor(cell_num,:))
            end
    end

    xlim([t0 t1]);
    ylim([0 length(spikes)+ 7 + length(csc)]);
    
    axeshandle = get(fighandle,'CurrentAxes');
    set(axeshandle,'yticklabel',[]);
    set(axeshandle,'ytick',[])
    
    x_orig = get(axeshandle, 'XLim');
    y_orig = get(axeshandle, 'YLim');
    
    for s = 1: length(csc)
        csc_mean = nanmean(Data(csc(s)));
        plot(Range(csc(s)),((Data(csc(s)).*0.003)+ ncells + 3 + s), 'Color', cscColor(s,:));
        plot([min(Range(csc(s))) max(Range(csc(s)))],[csc_mean.*0.003 csc_mean.*0.003]+ ncells + 3 + s,'Color', cscColor(s,:), 'LineStyle', '--','LineWidth', 0.5);
    end
    
    for e = 1:nEvents
        y_range = [min(Data(csc(1))).*0.0025 max(Data(csc(1))).*0.0025] + ncells + 4;
        plot([evt(e) evt(e)],y_range, 'Color', evtColor(e,:))
    end
    
    if interactiveMode == true        
        fprintf('\nInteractive Mode is on. Use arrow keys to navigate\n.');
        set(fighandle,'KeyPressFcn', @shifting);
    else
        fprintf('\nInteractive Mode is off\n.');
    end
