function  MOL_PlotPSTH(sessionData,trialData,spikeData)

% eventofinterest = 'stimStart';
eventofinterest = 'stimChange';
% eventofinterest = 'rewardTime';
% eventofinterest = 'lick_firstofresponse';
% eventofinterest = 'lick_firstofitibout';

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);

    params = params_histresponse;
    
    for iNeuron = 1:length(tempspikeData.ts)
        if strfind(eventofinterest,'lick')
            events_ts           = MOL_getLicks(temptrialData,snkprm.AlignOn);
        else
            events_ts           = temptrialData.(eventofinterest);          %Get events
        end

        spikes_ts           = tempspikeData.ts{iNeuron};                    %Get spikes for this neuron
        [edges,hist_mat]    = calc_psth(events_ts,spikes_ts,params);        %Get histogram for this neuron
        [splits,colors]     = MOL_GUI_GetSplits(temptrialData);             %Get splits on the basis of the GUI
%         if OverrulewithITILicks
%             splits{1} = true(length(events_ts),1);
%         end
        MOL_Plot_Histmat(edges,hist_mat,splits,colors,params,tempspikeData.cell_ID{iNeuron})           %Plot the responses
    end
    
end

distFig('Not',1) %distribute figures over the screen
end
