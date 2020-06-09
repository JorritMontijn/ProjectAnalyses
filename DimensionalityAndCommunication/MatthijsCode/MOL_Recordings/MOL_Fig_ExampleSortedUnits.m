function MOL_Fig_ExampleSortedUnits(varargin)
%% Example session 1
rootdir             = 'F:\Data\MODDISCR\RawData\2.14\2017-09-29_13-29-18\klustaSort\probe1shank1\';
rootfile            = '2017-09-29_13-29-18_probe1_shank1';
nChans              = 32; %nchannels on probe
ProbeType           = 'A1x32-Poly2-5mm-50s-177-CM32';
%% Example session 2
rootdir             = 'F:\Data\CHDET\RawData\2003\2018-02-07_14-28-57\klustaSort\probe2shank1\';
rootfile            = '2018-02-07_14-28-57_probe2_shank1';
nChans              = 32; %nchannels on probe
ProbeType           = 'A1x32-Poly2-5mm-50s-177-CM32';
%% Example session 3
rootdir             = 'F:\Data\CHDET\RawData\1008\2019-03-06_12-04-15\klustaSort\probe1shank1\';
rootfile            = '2019-03-06_12-04-15_probe1_shank1';
nChans              = 64; %nchannels on probe
ProbeType           = 'A1x64-Poly2-6mm-23s-160';

%% Parameters:
params              = setParamMatthijs;

baseFilename        = fullfile(rootdir,rootfile);
kwikFilename        = [baseFilename '.kwik'];
datFilename         = [baseFilename '.dat'];

nSamples            = 48;
nExNeurons          = 10; %number of example neurons which is plotted
nWaves              = 200; %number of waves that are loaded to estimate mean waveform
TipDepth            = 800; %doesn't matter
[ChannelX,ChannelY] = MOL_GetChannelPosition(ProbeType,TipDepth); %get configuration of channels

%% Get info about clusters and make subselection:
fprintf(1, 'loading data from kwik and kwx...\n');
clusters            = double(hdf5read(kwikFilename, '/channel_groups/0/spikes/clusters/main'));
uClusters           = unique(clusters)';
cluster_group       = NaN(size(uClusters));
time_samples        = hdf5read(kwikFilename, '/channel_groups/0/spikes/time_samples');

for iClus = 1:length(uClusters)
    cluster_group(iClus) = hdf5read(kwikFilename, ['/channel_groups/0/clusters/main/' num2str(uClusters(iClus))], 'cluster_group');
end

goodClusters        = uClusters(cluster_group == 2);
selectedClusters    = randsample(goodClusters,nExNeurons); %returns a vector of k values sampled uniformly from vector population

FileInfo            = dir(datFilename);
Source              = memmapfile(datFilename, 'Format', {'int16', [nChans, (FileInfo.bytes/nChans/2)], 'x'});

%% Get the mean waveforms of selected clusters and make the figure:
figure; set(gcf,'units','normalized','Position',[0.05 0.2 0.75 0.5],'color','w')
colors              = rand(nExNeurons,3)*0.8; %get some random colors, a little darker than random

for iSub = 1:nExNeurons %loop over selected neurons and thus subplots:
    iClus           = selectedClusters(iSub);
    spikeSamples    = time_samples(clusters==iClus); % ts in samples, not in seconds
    waveforms       = NaN(nChans,nSamples,nWaves); %initialize waveforms
    spikeSamples    = randsample(spikeSamples,nWaves); %take rand subsample for meanwave
    
    for i=1:nWaves %loop over subselection
        waveforms(:,:,i) = Source.Data.x(1:nChans,spikeSamples(i)-params.extract_s_before:spikeSamples(i)+params.extract_s_after);
        waveforms(:,:,i) = waveforms(:,:,i) - repmat(mean(waveforms(:,1:12,i),2),1,nSamples); %subtract baseline
    end
    
    meanWaveforms(:,:) = mean(waveforms(:,:,:),3); %get mean wave
    subplot('Position',[(iSub-1)*1/nExNeurons 1/nExNeurons 0.1 0.8]); hold all; %make subplot
    for iCh = 1:nChans %plot for every channel:
        plot([1:nSamples]+ChannelX(iCh)*1.6,meanWaveforms(iCh,:)+ChannelY(iCh)*15,'Color',colors(iSub,:),'LineWidth',2)
    end
    %some figure cosmetics:
    set(gca,'XTick',[], 'XTickLabels', [],'YTick',[], 'YTickLabels', [])
    xlim([-40 90]);
    ylim([-10600 12000]);
    box on
end

end
