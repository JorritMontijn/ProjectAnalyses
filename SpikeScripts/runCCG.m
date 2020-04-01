%pre-process data
%for intMouse=1:8
intMouse=1;
close all
clearvars -except intMouse
%get block data
if ~exist('cellMultiSes','var')
	if intMouse == 1
		strSes = '20140207';
	elseif intMouse == 2
		strSes = '20140314';
	elseif intMouse == 3
		strSes = '20140425';
	elseif intMouse == -1
		%return; %exclude?; bad behavior, weird signals
		strSes = '20140430';
	elseif intMouse == 4
		strSes = '20140507';
	elseif intMouse == 5
		strSes = '20140530';
	elseif intMouse == 6
		strSes = '20140604';
	elseif intMouse == 7
		strSes = '20140711';
	elseif intMouse == 8
		strSes = '20140715';
	end
	fprintf('Loading pre-processed data for %s [%s]\n',strSes,getTime);
	load(['D:\Data\Results\spikeAnalysis\dataPreProAggregate' strSes '.mat']);
end

%vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
vecBlockTypes = unique(vecBlock);
intNumBlocks = length(vecBlockTypes);
%vecNeuronNum = zeros(1,intNumBlocks);
%cellKeepList = cell(1,intNumBlocks);
%#ok<*ASGLU>
%#ok<*AGROW>
clear sLoad sSesAggregate ses
sParams.strFigDir = ['D:\Data\Results\spikeAnalysis' filesep strSes filesep];
sParams.boolSavePlots = true;
sParams.boolSaveData = true;
strOldDir = cd(sParams.strFigDir);

%change name for no split
strSes = ['AE' strSes];

%for intPopulation = vecBlockTypes
intPopulation=2;
%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
fprintf('Removing untuned neurons and reordering data for %s pop %d [%s]\n',strSes,intPopulation,getTime);

%take opposite directions as the same
cellMultiSes{intPopulation}.structStim.Orientation = mod(cellMultiSes{intPopulation}.structStim.Orientation,180);
vecOrientations = unique(cellMultiSes{intPopulation}.structStim.Orientation);
intNeurons = numel(cellMultiSes{intPopulation}.neuron);

%get neuronal responses per trial
[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
intContrasts = length(cellSelectContrasts);

%get orientation-based trial selection vectors
sTypesOri = getStimulusTypes(cellMultiSes{intPopulation},{'Orientation'});
cellSelectOri = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesOri);

%make normalized activity matrix per contrast
matTrialNormResponse = zeros(size(matTrialResponse));
for intContrast=1:intContrasts
	matTrialNormResponse(:,cellSelectContrasts{intContrast}) = imnorm(matTrialResponse(:,cellSelectContrasts{intContrast}));
end

%calculate preferred stimulus orientation for all neurons
vecOriPrefContrasts = 4:6; % 8% - 100%
matPref = nan(length(vecOriPrefContrasts),intNeurons);
intCounter = 0;
for intContrast = vecOriPrefContrasts
	intCounter = intCounter + 1;
	matResp = matTrialNormResponse(:,cellSelectContrasts{intContrast});
	structStimC{intContrast} = cellMultiSes{intPopulation}.structStim;
	cellFields = fieldnames(cellMultiSes{intPopulation}.structStim);
	for intField=1:length(cellFields)
		strField = cellFields{intField};
		structStimC{intContrast}.(strField) = structStimC{intContrast}.(strField)(cellSelectContrasts{intContrast});
	end
	cellSelect = getSelectionVectors(structStimC{intContrast},sTypesOri);
	sTuning{intContrast} = calcTuningRespMat(matResp,cellSelect,vecOrientations);
	matPref(intCounter,:) = sTuning{intContrast}.vecPrefIndex;
end
vecNeuronPrefStim = nan(1,intNeurons);
for intOri=1:length(vecOrientations);
	vecNeuronPrefStim(sum(matPref == intOri,1) > 1) = intOri;
end
%vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = randi(length(vecOrientations),[1 sum(isnan(vecNeuronPrefStim))]);

%remove non-tuned neurons
%{
cellMultiSes{intPopulation}.neuron(isnan(vecNeuronPrefStim)) = [];
vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = [];
%}
intNeurons = numel(cellMultiSes{intPopulation}.neuron);
intTrials = length(cellMultiSes{intPopulation}.structStim.Orientation);

%get updated response matrices
%get neuronal responses per trial
[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);

%% pre-pro to spikes
fprintf('Transforming dF/F0 data to activation events [%s]\n',getTime);
[cellSpikes,structStimSpikes] = getSpikeArrays(cellMultiSes{intPopulation});
ses = cellMultiSes{intPopulation};
fprintf('Transformation complete; retrieving population spiking data [%s]\n',getTime);

% get mean population activity
dblStep = 1/ses.samplingFreq;
dblTotTime = structStimSpikes.TimeOff(end)+5;
matSpikeCounts = getSpikeCounts(cellSpikes,0:dblStep:dblTotTime,dblStep);
fprintf('Population spiking data retrieved [%s]\n',getTime);


%% do analysis
%single-neuron CCGs
%{
for intNeuron=1:intNeurons
	figure
	
	vecCCG_bins = -1:dblStep:1;
	
	vecSpikeTimes = cellSpikes{intNeuron};
	
	intNeurons = length(cellSpikes);
	matCCG = zeros(intNeurons,length(vecCCG_bins));
	for intSpike=1:length(vecSpikeTimes)
		vecStart = vecCCG_bins + vecSpikeTimes(intSpike);
		matCCG = matCCG + getSpikeCounts(cellSpikes,vecStart,dblStep);
	end
	
	subplot(2,2,1)
	matCCG2 = conv2(matCCG,vecBlur,'valid');
	matCCG2 = imnorm(matCCG2,1);
	imagesc(matCCG2);
	colormap('hot')
	
	subplot(2,2,2)
	matCCG3 = matCCG2./repmat(sum(matCCG2,2),[1 size(matCCG2,2)]);
	vecNeurons=1:intNeurons;
	vecNeurons(intNeuron) = [];
	dblMaxWithoutAuto = max(max(matCCG3(vecNeurons,:)));
	%[dummy,vecReorder] = sort(sum(matCCG3.*repmat(1:size(matCCG3,2),[intNeurons 1]),2),'ascend');
	[dummy,vecReorder] = sort(sum(cumsum(matCCG3,2)>0.5,2),'descend');
	
	imagesc(cumsum(matCCG3(vecReorder,:),2))
	colormap('hot')
	
	subplot(2,2,3)
	imagesc(matCCG2(vecReorder,:))
	colormap('hot')
	
	subplot(2,2,4)
	imagesc(matCCG3(vecReorder,:),[0 dblMaxWithoutAuto])
	colormap('hot')
	
	pause
	close
end
%}

%pop-act based CCGs
vecCCG_bins = -1:dblStep:1;
matCCG = zeros(intNeurons,length(vecCCG_bins));
intNeurons = length(cellSpikes);

vecCCG_bins = -1:dblStep:1;
for intNeuron=1:intNeurons
	fprintf('Now at neuron %d/%d\n',intNeuron,intNeurons);
	vecSpikeTimes = cellSpikes{intNeuron};
	for intSpike=1:length(vecSpikeTimes)
		vecStart = vecCCG_bins + vecSpikeTimes(intSpike);
		matThisCCG = getSpikeCounts(cellSpikes,vecStart,dblStep);
		matThisCCG(intNeuron,:)=0;
		matCCG = matCCG + matThisCCG;
	end
end
matCCG(end+1,:) = sum(matCCG,1);

%plot mean pop activity
figure
subplot(2,2,1)
vecPopAct = sum(matSpikeCounts,1);
vecBlur = normpdf(-4:4,0,2);
vecBlur = vecBlur./sum(vecBlur);
vecPopAct = conv(vecPopAct,vecBlur,'same');
plot(vecPopAct)
xlim([0 dblTotTime/dblStep])

%try taking only events where pop act is above shuffled spiking prob
subplot(2,2,2)
vecBlur = 1;
matCCG2 = conv2(matCCG,vecBlur,'valid');
matCCG2 = imnorm(matCCG2,1);
imagesc(matCCG2);
colormap('hot')

subplot(2,2,3)
matCCG3 = matCCG2./repmat(sum(matCCG2,2),[1 size(matCCG2,2)]);

vecMuCC = (sum((matCCG.*repmat(1:size(matCCG,2),[size(matCCG,1) 1])),2) ./ sum(matCCG,2) - ceil(length(vecCCG_bins)/2)) * dblStep;
[dummy,vecReorder] = sort(vecMuCC,'ascend');
%[dummy,vecReorder] = sort(sum(cumsum(matCCG3,2)>0.5,2),'descend');

imagesc(cumsum(matCCG3(vecReorder,:),2))
colormap('hot')

subplot(2,2,4)
imagesc(matCCG2(vecReorder,:))
colormap('hot')


figure
subplot(2,2,1)
dblMaxWithoutAuto = max(max(matCCG3));
imagesc(matCCG3(vecReorder,:),[0 dblMaxWithoutAuto]);
colormap('hot');
freezeColors;

subplot(2,2,2)
matCorr = corr(matCCG3');
imagesc(matCorr,[-1 1]);
colormap('redblue');
freezeColors;

%take only packets

