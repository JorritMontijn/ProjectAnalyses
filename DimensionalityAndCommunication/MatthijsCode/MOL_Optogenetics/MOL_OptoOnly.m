function MOL_OptoOnly(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% Analysis of the effects of optical stimulation only:
% ToDo:
% -

%% Parameter settings for PSTH
params = params_hist_optoresponse(); % All time is in microseconds

% params.histmethod           = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.eventofinterest      = 'photostimStart';

meanbaseline            	= [];
meanopto                    = [];
sembaseline             	= [];
semopto                     = [];
signmodulated               = [];

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    if ~isempty(tempspikeData.ts)
        %init output histogram matrix:
        hist_mat = NaN(length(temptrialData.(params.eventofinterest)),(params.t_post-params.t_pre)/params.binsize,length(tempspikeData.ts));
        for iNeuron = 1:length(tempspikeData.ts)
            events_ts                   = temptrialData.(params.eventofinterest);       %Get events
            spikes_ts                   = tempspikeData.ts{iNeuron};                    %Get spikes for this neuron
            [edges,hist_mat(:,:,iNeuron)]    = calc_psth(events_ts,spikes_ts,params);   %Get histogram per neuron
        end
        
        %% Select only neurons that are present during the recording:
        if isfield(tempspikeData,'coverage')
            idx             = tempspikeData.coverage>0.5;
            hist_mat        = hist_mat(:,:,idx);
        end
        
        %% Calculate response:
        %Calculate baseline:
        baseline                = squeeze(mean(hist_mat(:,edges>params.twin_baseline_start & edges<=params.twin_baseline_stop,:),2));
        %Get response:
        opto                    = squeeze(mean(hist_mat(:,edges>params.twin_resp_start & edges<=params.twin_resp_stop,:),2));
        
        meanbaseline            = [meanbaseline nanmean(baseline,1)];
        meanopto                = [meanopto nanmean(opto,1)];
        %     sembaseline             = [sembaseline nanstd(baseline,1)/sqrt(numel(baseline)-1)];
        %     semopto                 = [semopto nanstd(opto,1)/sqrt(numel(opto)-1)];
        sembaseline             = [sembaseline nanstd(baseline,1)];
        semopto                 = [semopto nanstd(opto,1)];
        
        mod                     = zeros(1,size(baseline,2));
        mod(ttest(baseline,opto,'alpha',0.025,'tail','right')==1)            = -1;
        mod(ttest(baseline,opto,'alpha',0.025,'tail','left')==1)             = 1;
        signmodulated           = [signmodulated mod]; %#ok<*AGROW>
    end
end

%% Cap values at 0.01 Hz:
meanopto(meanopto<0.01) = 0.01;
meanbaseline(meanbaseline<0.01) = 0.01;

%% Convert to log scale:
sembaseline             = log10(meanbaseline+sembaseline)-log10(meanbaseline);  %Only used when showing errorbars
semopto                 = log10(meanopto+semopto)-log10(meanopto);              %Only used when showing errorbars

meanbaseline            = log10(meanbaseline);
meanopto                = log10(meanopto);

%% Make scatter figure
figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w'); hold on;
% title('Optical stimulation in V1','FontSize',25)

% scatter(meanbaseline,meanopto,40,'b','filled'); %Plot dots
scatter(meanbaseline(signmodulated==-1),meanopto(signmodulated==-1),40,'b','filled'); %Plot dots
scatter(meanbaseline(signmodulated==0),meanopto(signmodulated==0),40,'k','filled'); %Plot dots
scatter(meanbaseline(signmodulated==1),meanopto(signmodulated==1),40,'g','filled'); %Plot dots

% errorbar(meanbaseline, meanopto, semopto, 'k.');
% herrorbar(meanbaseline, meanopto, sembaseline,'k.');
plot([-5 5],[-5 5],'k','LineWidth',3); %Plot equity line

%Make up of figure:
grid on
ticks = [0.01 0.03 0.1 0.3 1 3 10 30 100];
set(gca,'XTick',log10(ticks),'XTickLabels',ticks,'FontSize', 20)
set(gca,'YTick',log10(ticks),'YTickLabels',ticks,'FontSize', 20)
set(gca,'linewidth',3)
xlim([min(log10(ticks)) max(log10(ticks))*1.1])
ylim([min(log10(ticks)) max(log10(ticks))*1.1])
xlabel('BASELINE (Hz)','FontSize', 20)
ylabel('OPTO (Hz)','FontSize', 20)

fprintf('n = %d neurons\n',numel(meanbaseline));

%% Check if significantly enhanced firing neurons are also classified as interneurons:
try %Get Interneurons and Pyramidal Neurons:
    [spikeData.celltype,spikeData.celltypelabel]      = calc_neurontype([spikeData.meanwf{:}]);
catch
    nNeurons = length(spikeData.meanwf);
    for iNeuron = 1:nNeurons
        spikeData.meanwf{iNeuron} = spikeData.meanwf{iNeuron}(1:32);
    end
    [spikeData.celltype,spikeData.celltypelabel]      = calc_neurontype([spikeData.meanwf{:}]);
end
uniqsignmod = unique(signmodulated);
uniqcelltype = unique(spikeData.celltype);
for iMod = 1:length(uniqsignmod)
    for iCType = 1:length(uniqcelltype)
        barmat(iMod,iCType) = sum(spikeData.celltype(signmodulated==uniqsignmod(iMod))==uniqcelltype(iCType));
    end
end
figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w'); hold on;
bar(barmat);
ylabel('Number of neurons', 'FontSize',20)

end

% %% To laser stim onset/offset:
%
% npulses = 20; dur = 0.05;
% laserstimonset = repmat(trialData.photostimonset,1,npulses);
% laserstimonset = laserstimonset + repmat(linspace(0,(npulses-1)*dur,npulses),length(trialData.photostimonset),1)*1e6;
% % [1 2 4 5 8];
% showFig = 1; close all;
% for iCell = 1:length(spikeData.cell_ID)
%     params              = params_optoresponse; %Get parameters
%     [Temp_Resp(iCell),Temp_Baseline(iCell),Temp_Onset(iCell)]    = calc_response(laserstimonset(:)',spikeData.ts{iCell,1},params,showFig); %#ok<*SAGROW> %Get response
% end
% distFig('Screen','Ext'); %Distribute figures over screen
%
% %% Make scatter figure
% figure; set(gcf,'color','w'); hold on;
% title('Optical inhibition of AL/PM')
% scatter(Temp_Baseline,Temp_Resp,15,'b','filled');
% maximum = max([Temp_Baseline Temp_Resp]);
% plot([0 maximum],[0 maximum],'k','LineWidth',2)
% xlabel('BASELINE')
% ylabel('OPTO')
%
% %% Show neuron type:
% figure; colors  = 'brg';
% for iCell = 1:length(spikeData.cell_ID)
%     if any(ismember([1 2 4 5 8],iCell))
%         plot(linspace(0,2,64), mean(spikeData.wf{iCell,1},1),'b');        hold on
%     else
%         plot(linspace(0,2,64), mean(spikeData.wf{iCell,1},1),'r');        hold on
%     end
% end

