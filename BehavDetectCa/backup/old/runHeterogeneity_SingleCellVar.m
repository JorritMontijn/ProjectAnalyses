
%starting message
clc;
fprintf('\n		WELCOME TO NAT MOV ANAL :: THE NATURAL MOVIE ANALYSER\n\n')
fprintf('Process started at [%s] %s. Please be patient.\n\n',getTime,date);

%load session file database
strDataDir = 'D:\Data\Guido\';
strOutputDir = 'D:\Data\Results\correlationAnalysis\';
sLoad = load([strDataDir 'SessionDatabase.mat']);
cellRecGratings = sLoad.cellRecGratings;
cellRecNatMovs = sLoad.cellRecNatMovs;

%get which mouse/mice to process
strPath = mfilename('fullpath');
strPath = strPath(1:(end-length(mfilename)));

boolSavePlots = true;
for intMouse=91:95
	%mouse name
	strMouse = num2str(intMouse);
	fprintf('Retrieving data for mouse %s [%s]\n',strMouse,getTime);
	
	[sSesAggregate,indInclude] = getNatMovAnalPrep(strDataDir, cellRecNatMovs{intMouse});
	
	%get info
	sTypes = getStimulusTypes(sSesAggregate);
	cellSelect = getSelectionVectors(sSesAggregate.structStim,sTypes);
	intTrials = length(sSesAggregate.structStim.Scene);
	intNeurons = numel(sSesAggregate.neuron);
	intCombinations = (intNeurons *(intNeurons -1))/2;
	intStimTypes = length(cellSelect);
	intRepetitions = intTrials/intStimTypes;
	
	%get response matrix
	matAllResp = cat(1,sSesAggregate.neuron(:).dFoF);
	vecTrialStarts = sSesAggregate.structStim.FrameOn(sSesAggregate.structStim.Scene==1);
	vecTrialStops = sSesAggregate.structStim.FrameOff(sSesAggregate.structStim.Scene==1);
	
	
	[vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matAllResp);
	
	%% make PSTHs; rep x T x N
	vecSelectT = 11;
	vecN = 3;
	intTimePoints = unique(vecTrialStops - vecTrialStarts + 1);
	matHetPSTH = nan(intRepetitions,intTimePoints);
	matRespPSTH  = nan(intRepetitions,intTimePoints,intNeurons);
	for intRep = 1:intRepetitions
		matHetPSTH(intRep,:) = vecHeterogeneity(vecTrialStarts(intRep):vecTrialStops(intRep));
		matRespPSTH(intRep,:,:) = matAllResp(:,vecTrialStarts(intRep):vecTrialStops(intRep))';
	end
	%matRespPSTH_zScored = bsxfun(@rdivide,abs(bsxfun(@rdivide,bsxfun(@minus,matRespPSTH,mean(matRespPSTH,1)),std(matRespPSTH,[],1))),max(matRespPSTH,[],3));
	matRespPSTH_zScored = bsxfun(@rdivide,bsxfun(@minus,matRespPSTH,mean(matRespPSTH,1)),std(matRespPSTH,[],1));
	matHetPSTH = matHetPSTH(:,vecSelectT);
	matRespPSTH = matRespPSTH(:,vecSelectT,vecN);
	matRespPSTH_zScored = matRespPSTH_zScored(:,vecSelectT,vecN);
	matSDPSTH = mean(matRespPSTH_zScored,3);
	scatter(matHetPSTH(:),matSDPSTH(:),'xk')
	
	dblStep = 0.1;
	vecBin = 0.6:dblStep:2;
	vecBinPlot = (0.6+dblStep/2):dblStep:(2-dblStep/2);
	[nVec,meanVec,stdVec] = makeBins(matHetPSTH(:),matSDPSTH(:),vecBin);
	
	hold on
	errorbar(vecBinPlot,meanVec,stdVec,'r')
	hold off
end