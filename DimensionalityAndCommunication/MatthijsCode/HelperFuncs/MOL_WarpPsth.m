function psth_out = MOL_WarpPsth(psth_in,edges,tstart,tend,t_events,t_align)

if isnan(t_align)
   t_align = nanmedian(t_events);
end

t_start_idx     = find(edges>tstart,1,'first');
t_end_idx       = find(edges>tend,1,'first');
t_align_idx     = find(edges>t_align,1,'first');

psth_out                    = NaN(size(psth_in)); %init psth out
psth_out(:,1:t_start_idx-1) = psth_in(:,1:t_start_idx-1);
psth_out(:,t_end_idx:end)   = psth_in(:,t_end_idx:end);

psth_out                    = psth_in; %init psth out
% psth_out(:,1:t_start_idx-1) = psth_in(:,1:t_start_idx-1);
% psth_out(:,t_end_idx:end)   = psth_in(:,t_end_idx:end);

for iTr = 1:size(psth_in,1)
    if t_events(iTr)>0.1e6 && t_events(iTr)<1e6
        t_ev_idx        = find(edges>t_events(iTr),1,'first');
        warpedbefore 	= interp1(linspace(0,1,t_ev_idx-t_start_idx),psth_in(iTr,t_start_idx:t_ev_idx-1),linspace(0,1,t_align_idx-t_start_idx),'linear');
        warpedafter     = interp1(linspace(0,1,t_end_idx+1-t_ev_idx),psth_in(iTr,t_ev_idx:t_end_idx),linspace(0,1,t_end_idx+1-t_align_idx),'linear');
        psth_out(iTr,t_start_idx:t_align_idx-1) = warpedbefore;
        psth_out(iTr,t_align_idx:t_end_idx) = warpedafter;
    end
        
end

end