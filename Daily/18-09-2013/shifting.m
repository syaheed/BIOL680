function shifting(src,event)

    axeshandle = get(src, 'CurrentAxes');
    x_orig = get(axeshandle, 'XLim');
    x_range = x_orig(1)-x_orig(2);
    
    if strcmp(event.Key, 'leftarrow')==1
        x_new = x_orig + (x_range/2);
    elseif strcmp(event.Key, 'rightarrow')==1
        x_new = x_orig - (x_range/2);
    elseif strcmp(event.Key, 'uparrow')==1
        x_range_new = x_range ./2;
        x_new = [x_orig(1)+x_range_new x_orig(1)];
    elseif strcmp(event.Key, 'downarrow')==1
        x_range_new = x_range * 2;
        x_new = [x_orig(2) x_orig(2)-x_range_new];
    else
        x_new = x_orig;
    end
    
    set(axeshandle, 'XLim', x_new);