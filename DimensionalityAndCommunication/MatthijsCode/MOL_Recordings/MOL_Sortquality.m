function MOL_Sortquality(varargin)

% sessionData     = varargin{1};
% trialData       = varargin{2};
spikeData       = varargin{3};

%% Parameters
nExNeurons          = 25; %number of example neurons which is plotted
params.minID        = 10; %minimum isolation distance
params.maxISI_FA    = 0.015; %maximum fraction of spikes within ISI of 1.5? ms
params.maxIDdivISI  = 20; %Minimum IsoDis divided by ISIFA

%% get selection based on criteria
nTotalNeurons       = length(spikeData.cell_ID);

idx_IsoDist     = spikeData.QM_IsolationDistance>params.minID;
idx_ISI_FA      = spikeData.QM_ISI_FA<params.maxISI_FA;

TotalNeurons    = 1:nTotalNeurons;

selectedNeurons     = randsample(TotalNeurons(idx_IsoDist&idx_ISI_FA),nExNeurons); %returns a vector of k values sampled uniformly from vector population
       

%% ISI figure:
figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')

edges = [0:0.5:50]*1e3; %#ok<*NBRAK>
for iNeuron = 1:nExNeurons
    subplot(5,5,iNeuron);
    isihist = histcounts(diff(spikeData.ts{selectedNeurons(iNeuron)}),edges,'Normalization','probability');
    bar(edges(1:end-1),isihist);
    xlim([-2 50]*1e3);
    set(gca,'XTick',[0 10 20 30 40] *1e3, 'XTickLabels', [0 10 20 30 40],'YTick',[], 'YTickLabels', [])
    xlabel('ISI in ms');
    ylim([0 max(isihist)*1.1])
    patch([-2 2.5 2.5 -2]*1e3,[0 0 max(isihist)*1.1 max(isihist)*1.1],'k','FaceAlpha',.1,'EdgeColor','none')
end

%% Parameters
figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')

ymax        = 60;
xmax        = 0.03;

scatter(spikeData.QM_ISI_FA,spikeData.QM_IsolationDistance,250,'k.' )
% patch([0 params.maxISI_FA params.maxISI_FA 0],[params.minID params.minID ymax ymax],'g','FaceAlpha',.3,'EdgeColor','none')
% fprintf('Percent rejected: %2.2f%%\n',sum(spikeData.QM_IsolationDistance<params.minID | spikeData.QM_ISI_FA>params.maxISI_FA)/length(spikeData.QM_IsolationDistance)*100)

patch([0 xmax xmax 0],[0 params.maxIDdivISI*xmax*100 params.maxIDdivISI*xmax*100 params.maxIDdivISI*xmax*100 ],'g','FaceAlpha',.3,'EdgeColor','none')
fprintf('Percent rejected: %2.2f%%\n',sum((spikeData.QM_IsolationDistance./(spikeData.QM_ISI_FA*100))<params.maxIDdivISI)/length(spikeData.QM_IsolationDistance)*100)

ylim([0 ymax])
xlim([0 xmax])

set(gca,'XTick',0:0.005:xmax, 'XTickLabels', [0:0.005:xmax]*100,'YTick',0:5:ymax,'FontSize',22)
grid on

ylabel('Isolation Distance')
xlabel('% Spikes within 1.5 ms')

end
