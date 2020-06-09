function MOL_LickHistogram(sessionData,trialData,spikeData)

params          = params_histresponse();
showHistogram   = 1;

hist_mat_all    = [];
splits_all      = {[] []};

params.conv_twin            = 0.5e6;        %Window size for smoothing
params.conv_sigma           = 0.05e6;        %sd of gaussian window for smoothing

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [tempsessionData,temptrialData] = MOL_getTempPerSes(sesid,sessionData,trialData);
    
    events_ts       = temptrialData.stimChange;
%     events_ts       = temptrialData.stimStart;
    
    licktimes       = [temptrialData.lickTime{:}];
    if tempsessionData.VisualLeftCorrectSide
        licksvis      = licktimes(strfind([temptrialData.lickSide{:}],'L'));
        licksaud      = licktimes(strfind([temptrialData.lickSide{:}],'R'));
    else
        licksvis      = licktimes(strfind([temptrialData.lickSide{:}],'R'));
        licksaud      = licktimes(strfind([temptrialData.lickSide{:}],'L'));
    end
    
    if ~(sum([numel(licksvis) numel(licksaud)]) == numel(licktimes))
        error('no match between licktimes and licksides')
    end

    licks = licksvis;
%     licks = licksaud;
    [edges,hist_mat]    = calc_psth(events_ts,licks,params);        %Get histogram for this neuron
    [splits,colors]     = MOL_GUI_GetSplits(temptrialData);             %Get splits on the basis of the GUI

    for i = 1:length(splits)
        splits_all{i} = [splits_all{i}; splits{i}];
    end
    hist_mat_all = [hist_mat_all; hist_mat];
    
%     licksright = licksright * 1000000;
%     events_ts = events_ts * 1000000;
%     licksleft = licksleft * 1000000;
%     [resp,baseline,tmax] = calc_response(events_ts,licksright,params,showHistogram);
%     [resp,baseline,tmax] = calc_response(events_ts,licksleft,params,showHistogram);
%     [resp,baseline,tmax] = calc_response(events_ts,responses_ts,params,showHistogram);
%     eval(sprintf('%s(tempsessionData,temptrialData)',script));

end

% figure;
% set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
figtitle = 'Licks';
MOL_Plot_Histmat(edges,hist_mat,splits,colors,params,figtitle)           %Plot the responses
ylim([0 7])


end