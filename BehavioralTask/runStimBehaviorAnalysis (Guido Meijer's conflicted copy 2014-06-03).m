%% define variables
clear all;
close all;
intMouse = 1;
if intMouse == 1
	strSes = '20140207';
	strMasterPath = 'D:\Data\Processed\imagingdata';
	vecRecordings = 1:8;
	vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
elseif intMouse == 2
	strSes = '20140314';
	strMasterPath = 'D:\Data\Processed\imagingdata';
	vecRecordings = 1:7;
	vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
elseif intMouse == 3
	strSes = '20140425';
	strMasterPath = 'D:\Data\Processed\imagingdata';
	vecRecordings = 1:8;
	vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
elseif intMouse == 4
	strSes = '20140430';
	strMasterPath = 'D:\Data\Processed\imagingdata';
	vecRecordings = 2:8;
	vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
elseif intMouse == 5
	strSes = '20140507';
	strMasterPath = 'D:\Data\Processed\imagingdata';
	vecRecordings = 1:3;
	vecBlock = ones(size(vecRecordings)); %define whether neurons are in same population or not
end

%put in structure
sIn.strSes = strSes;
sIn.strMasterPath = strMasterPath;
sIn.vecRecordings = vecRecordings;

%define parameters
intAgg = 0; %will calculate everything over aggregate session if ~=0
sParams = struct;
sParams.strFigDir = ['D:\Data\Results\stimdetection' filesep strSes filesep];
strOldDir = cd(sParams.strFigDir);
sParams.boolSplitData = strcmp(strSes(end-4:end),'split');
sParams.boolSavePlots = false;

%which analyses to perform
sParams.boolDoPsychoCurve = true;
sParams.boolDoTraceStimPEP = true;
sParams.boolDoTraceRespPEP = true;
sParams.boolDoStimFrames = true;
sParams.boolDoStimMinusPreStimFrames = true;
sParams.boolDoPreStimFrames = true;
sParams.boolDoActivityDistributions = true; %memory intensive

close

sParams.boolDoPrecedingAlphaPower = false;
sParams.boolDoStimFramesDecoding = false;
%calculate block definitions
intStartSes = 1;
intStopSes = length(vecRecordings);
vecBlockTypes = unique(vecBlock);
vecBiDirSes = [nan 9];
intNumBlocks = length(vecBlockTypes);
vecNeuronNum = zeros(1,intNumBlocks);
cellKeepList = cell(1,intNumBlocks);
vecFirstBlock = length(vecBlockTypes);
vecLastBlock = length(vecBlockTypes);
for intBlockType = vecBlockTypes
	vecFirstBlock(intBlockType) = find(vecBlock == intBlockType,1,'first');
	vecLastBlock(intBlockType) = find(vecBlock == intBlockType,1,'last');
end
if sParams.boolSplitData,vecBlock = 1:intStopSes;end

%% get stim aggregate
fprintf('Building aggregate stimulus structure...\n')
sStimAggregate = buildStimAggregate(sIn);

%% calculate psychometric curve on % responded and RT
fprintf('Calculating psychometric curve...\n')
sPC = struct;
sPC = calcPsychometricCurve(sStimAggregate);
%matChiSquareP = sPC.matP;
cellSaveBehavDetect = shiftdim(sPC.cellResponses,1); %behavioral detection [6 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
cellSaveBehavRT = cell(5,1); %behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
for intC = 1:5
	vecRT = sPC.matRTs(:,intC+1);
	cellSaveBehavRT{intC} = vecRT(~isnan(vecRT));
end

%% get stim detection data per session
%create links to files
intCounter = 0;
cellAggregate = cell(1,length(vecRecordings));
for intRec=vecRecordings
	intCounter = intCounter + 1;
	cellAggregate{intCounter} = sprintf('%s%s%s%sxyt%02d%s%sxyt%02d_ses',strMasterPath,filesep,strSes,filesep,intRec,filesep,strSes,intRec);
end

%pre-allocate
cellSes = cell(1,numel(cellAggregate));
sContrastDataStim = struct;
sMetaDataStim = struct;
sContrastDataSB = struct;
sMetaDataSB = struct;
sContrastDataResp = struct;
sMetaDataResp = struct;
sContrastDataSF = struct;
sMetaDataSF = struct;
sContrastDataPS = struct;
sMetaDataPS = struct;
sContrastDataSD = struct;
sMetaDataSD = struct;
sContrastDataAP = struct;
sMetaDataAP = struct;

%aggregate all behavioral sessions
fprintf('Loading all data and creating neuron inclusion list..\n')
for intThisBlock=1:intNumBlocks
	cellPresenceList{intThisBlock} = [];
	cellResponseList{intThisBlock} = [];
end
for intFile=1:numel(cellAggregate)
	%load data
	sLoad = load([cellAggregate{intFile} '.mat']);
	ses = sLoad.ses;
	cellSes{intFile} = ses;
	
	%create inclusion list
	intThisBlock = vecBlock(intFile);
	vecNeuronNum(intThisBlock) = numel(ses.neuron);
	indPresenceList = cellPresenceList{intThisBlock};
	indResponseList = cellResponseList{intThisBlock};
	if isempty(indPresenceList),indPresenceList = true(1,numel(ses.neuron));end
	if isempty(indResponseList),indResponseList = false(1,numel(ses.neuron));end
	for intNeuron=1:numel(ses.neuron)
		if strcmp(ses.neuron(intNeuron).strPresence,'absent')
			indPresenceList(intNeuron) = false;
		end
		if ~strcmp(ses.neuron(intNeuron).strRespType,'silent')
			indResponseList(intNeuron) = true;
		end
	end
	cellPresenceList{intThisBlock} = indPresenceList;
	cellResponseList{intThisBlock} = indResponseList;
	
	%msg
	fprintf('Processed [%s]\n',cellAggregate{intFile})
end
for intThisBlock=1:intNumBlocks
	cellKeepList{intThisBlock} = cellPresenceList{intThisBlock} & cellResponseList{intThisBlock};
end

if ~sParams.boolSplitData
	%add list to parameter structure
	for intThisBlock=vecBlockTypes
		intBlockStart = find(vecBlock == intThisBlock,1,'first');
		intBlockStop = find(vecBlock == intThisBlock,1,'last');
		indKeepList = cellKeepList{intThisBlock};
		intNeuronNum = vecNeuronNum(intThisBlock);
		intNeuronsIncluded = sum(indKeepList);
		vecNeuronNum(intThisBlock) = intNeuronsIncluded;
		fprintf('>> Will include [%d/%d] neurons for block %d [%d - %d]\n',intNeuronsIncluded,intNeuronNum,intThisBlock,intBlockStart,intBlockStop)
	end
else
	cellKeepList = cell(1,intNumBlocks);
end

%% activity distributions
if sParams.boolDoActivityDistributions || intAgg ~= 0
	sParamsAD.intStartSes = intStartSes;
	sParamsAD.intStopSes = intStopSes;
	
	%build multi ses aggregate structure combination
	fprintf('\nBuilding multi-recording aggregate per block...\n')
	cellMultiSes = {};
	if ~sParams.boolSplitData
		cellMultiSes = cell(1,intNumBlocks);
		for intBlockType = vecBlockTypes
			sSesAggregate = [];
			for intSes=vecFirstBlock(intBlockType):vecLastBlock(intBlockType)
				%load session file
				ses = cellSes{intSes};
				
				%remove all excluded neurons
				indKeepList = cellKeepList{intBlockType};
				ses = doRecalcdFoF(ses,6,indKeepList);
				
				%create aggregate
				sSesAggregate = buildMultiSesAggregate(ses,sSesAggregate);
			end
			cellMultiSes{intBlockType} = sSesAggregate;
		end
	end
	fprintf('\b    Done! Time is %s\n',getTime)
	
	%close all;
	%doPlotActivityDistributions%(cellContrastDataAD,sMetaDataAD,sPC,sParamsAD);
end

%perform analysis
if intAgg ~= 0
	intStartSes=1;
	intStopSes=1;
	vecFirstBlock = 1;
	vecLastBlock = 1;
end
for intSes=intStartSes:intStopSes
	if intAgg ~= 0
		ses = cellMultiSes{intAgg};
		intThisBlock = 1;
		indKeepList = true(1,numel(ses.neuron));
	else
		%load data
		ses = cellSes{intSes};

		%get inclusion list
		intThisBlock = vecBlock(intSes);
		indKeepList = cellKeepList{intThisBlock};
	end
	
	%get actual data
	sParams.indKeepList = indKeepList;
	sProcData = doStimBehaviorAnalysis(ses,sParams,sStimAggregate);
	
	%aggregate for different analyses
	%stim-locked
	if sProcData.sParams.boolDoTraceStimPEP
		sMetaDataStim = sProcData.sTS.sMetaData;
		cellSubFields = fieldnames(sProcData.sTS.sContrast);
		for intContrastIndex=1:length(sMetaDataStim.cellSelectC)
			for intField=1:length(cellSubFields)
				strField = cellSubFields{intField};
				if intContrastIndex > numel(sContrastDataStim) || ~isfield(sContrastDataStim(intContrastIndex),strField)
					sContrastDataStim(intContrastIndex).(strField) = [];
				end
				sContrastDataStim(intContrastIndex).(strField) = [sContrastDataStim(intContrastIndex).(strField); sProcData.sTS.sContrast(intContrastIndex).(strField)];
			end
		end
	end
	
	%resp-locked
	if sProcData.sParams.boolDoTraceRespPEP
		sMetaDataResp = sProcData.sTR.sMetaData;
		cellSubFields = fieldnames(sProcData.sTR.sContrast);
		for intContrastIndex=1:length(sMetaDataResp.cellSelectC)
			for intField=1:length(cellSubFields)
				strField = cellSubFields{intField};
				if intContrastIndex > numel(sContrastDataResp) || ~isfield(sContrastDataResp(intContrastIndex),strField)
					sContrastDataResp(intContrastIndex).(strField) = [];
				end
				sContrastDataResp(intContrastIndex).(strField) = [sContrastDataResp(intContrastIndex).(strField); sProcData.sTR.sContrast(intContrastIndex).(strField)];
			end
		end
	end
	
	%activity stim frames
	if sProcData.sParams.boolDoStimFrames
		sMetaDataSF = sProcData.sSF.sMetaData;
		cellSubFields = fieldnames(sProcData.sSF.sContrast);
		for intContrastIndex=1:length(sMetaDataSF.cellSelectC)
			for intField=1:length(cellSubFields)
				strField = cellSubFields{intField};
				if intContrastIndex > numel(sContrastDataSF) || ~isfield(sContrastDataSF(intContrastIndex),strField)
					sContrastDataSF(intContrastIndex).(strField) = [];
				end
				sContrastDataSF(intContrastIndex).(strField) = [sContrastDataSF(intContrastIndex).(strField) sProcData.sSF.sContrast(intContrastIndex).(strField)];
			end
		end
	end
	
	%activity stim frames minus baseline frames
	if sProcData.sParams.boolDoStimMinusPreStimFrames
		sMetaDataSB = sProcData.sSB.sMetaData;
		cellSubFields = fieldnames(sProcData.sSB.sContrast);
		for intContrastIndex=1:length(sMetaDataSB.cellSelectC)
			for intField=1:length(cellSubFields)
				strField = cellSubFields{intField};
				if intContrastIndex > numel(sContrastDataSB) || ~isfield(sContrastDataSB(intContrastIndex),strField)
					sContrastDataSB(intContrastIndex).(strField) = [];
				end
				sContrastDataSB(intContrastIndex).(strField) = [sContrastDataSB(intContrastIndex).(strField) sProcData.sSB.sContrast(intContrastIndex).(strField)];
			end
		end
	end
	
	%decoding stim frames
	if sProcData.sParams.boolDoStimFramesDecoding
		sMetaDataSD = sProcData.sSD.sMetaData;
		cellSubFields = fieldnames(sProcData.sSD.sContrast);
		for intContrastIndex=1:length(sMetaDataSD.cellSelectC)
			for intField=1:length(cellSubFields)
				strField = cellSubFields{intField};
				if intContrastIndex > numel(sContrastDataSD) || ~isfield(sContrastDataSD(intContrastIndex),strField)
					sContrastDataSD(intContrastIndex).(strField) = [];
				end
				sContrastDataSD(intContrastIndex).(strField) = [sContrastDataSD(intContrastIndex).(strField) sProcData.sSD.sContrast(intContrastIndex).(strField)];
			end
		end
	end
	
	%activity during preceding baseline frames
	if sProcData.sParams.boolDoPreStimFrames
		sMetaDataPS = sProcData.sPS.sMetaData;
		cellSubFields = fieldnames(sProcData.sPS.sContrast);
		for intContrastIndex=1:length(sMetaDataPS.cellSelectC)
			for intField=1:length(cellSubFields)
				strField = cellSubFields{intField};
				if intContrastIndex > numel(sContrastDataPS) || ~isfield(sContrastDataPS(intContrastIndex),strField)
					sContrastDataPS(intContrastIndex).(strField) = [];
				end
				sContrastDataPS(intContrastIndex).(strField) = [sContrastDataPS(intContrastIndex).(strField) sProcData.sPS.sContrast(intContrastIndex).(strField)];
			end
		end
	end
	
	%activity distributions
	if sParams.boolDoActivityDistributions
		sMetaDataAD{intSes} = sProcData.sAD.sMetaData;
		cellContrastDataAD{intSes} = sProcData.sAD.sContrast;
	end
	
	%preceding alpha power
	if sParams.boolDoPrecedingAlphaPower
		sMetaDataAP = sProcData.sAP.sMetaData;
		cellSubFields = fieldnames(sProcData.sAP.sContrast);
		for intContrastIndex=1:length(sMetaDataAP.cellSelectC)
			for intField=1:length(cellSubFields)
				strField = cellSubFields{intField};
				if intContrastIndex > numel(sContrastDataAP) || ~isfield(sContrastDataAP(intContrastIndex),strField)
					sContrastDataAP(intContrastIndex).(strField) = [];
				end
				sContrastDataAP(intContrastIndex).(strField) = [sContrastDataAP(intContrastIndex).(strField); sProcData.sAP.sContrast(intContrastIndex).(strField)];
			end
		end
	end
end


%% plot different data types
%% plot psychometric curve
if sParams.boolDoPsychoCurve
	doPlotPsychoMetricCurve(sStimAggregate,sPC,1,sParams.boolSavePlots,strSes);
end

%% stim-locked
if sProcData.sParams.boolDoTraceStimPEP
	vecWindow = sMetaDataStim.vecWindow;
	intWindowFrames = sMetaDataStim.intWindowFrames;
	
	cellSaveTraceAct = cell(2,2,6); %dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values

	for intContrastIndex=1:length(sMetaDataStim.cellSelectC)
		
		%figure
		handleFig = figure;
		set(handleFig,'Color',[1 1 1]);
		figure(handleFig);
		
		%pre-compute variables
		vecWindowSelect = vecWindow(1):vecWindow(end);
		intWL = length(vecWindowSelect);
		vecWindowInv = intWL:-1:1;
		vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
		vecX = [vecWindowSelect vecWindowPlotInv]/ses.samplingFreq;
		vecLineX = vecWindowSelect/ses.samplingFreq;
		dblContrast = sMetaDataStim.vecContrasts(intContrastIndex);
		
		for intPlot=1:4
			subplot(2,2,intPlot)
			if intPlot == 1
				matTemp = sContrastDataStim(intContrastIndex).matPR;
				cellSaveTraceAct{1,1,intContrastIndex} = matTemp;
				strTitle = sprintf('Contrast %.3f; Pref / Resp',dblContrast);
			elseif intPlot == 2
				matTemp = sContrastDataStim(intContrastIndex).matNPR;
				cellSaveTraceAct{2,1,intContrastIndex} = matTemp;
				strTitle = sprintf('Contrast %.3f; Non-Pref / Resp',dblContrast);
			elseif intPlot == 3
				matTemp = sContrastDataStim(intContrastIndex).matPN;
				cellSaveTraceAct{1,2,intContrastIndex} = matTemp;
				strTitle = sprintf('Contrast %.3f; Pref / No Resp',dblContrast);
			elseif intPlot == 4
				matTemp = sContrastDataStim(intContrastIndex).matNPN;
				cellSaveTraceAct{2,2,intContrastIndex} = matTemp;
				strTitle = sprintf('Contrast %.3f; Non-Pref / No Resp',dblContrast);
			end
			
			%get data
			vecMeanTrace = mean(matTemp,1);
			vecSE = std(matTemp,[],1)/sqrt(size(matTemp,1));
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			fill(vecX,vecY,[0.7 0.7 1],'EdgeColor','none');
			hold on
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',[0 0 1]);
			hold off
			title(strTitle)
			grid on
			xlabel('Time after stimulus onset (s)')
			ylabel('dF/F')
			xlim([-3 5])
			ylim([-0.02 0.06])
		end
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_stimPEP_c%.1f_raw',strSes,dblContrast*100);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
	end
end

%% resp-locked
if sProcData.sParams.boolDoTraceRespPEP
	%pre-compute variables
	vecWindow = sMetaDataResp.vecWindow;
	intWindowFrames = sMetaDataResp.intWindowFrames;
	
	vecWindowSelect = vecWindow(1):vecWindow(end);
	intWL = length(vecWindowSelect);
	vecWindowInv = intWL:-1:1;
	vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
	vecX = [vecWindowSelect vecWindowPlotInv]/ses.samplingFreq;
	vecLineX = vecWindowSelect/ses.samplingFreq;
	
	
	for intContrastIndex=1:length(sMetaDataResp.cellSelectC)
		
		%get contrast
		dblContrast = sMetaDataResp.vecContrasts(intContrastIndex);
		if dblContrast == 0,continue;end
		
		%figure
		handleFig = figure;
		set(handleFig,'Color',[1 1 1]);
		figure(handleFig);
		
		for intPlot=1:3
			subplot(2,2,intPlot)
			if intPlot == 1
				matTemp = sContrastDataResp(intContrastIndex).matPR;
				strTitle = sprintf('Contrast %.3f; Pref / Resp',dblContrast);
			elseif intPlot == 2
				matTemp = sContrastDataResp(intContrastIndex).matNPR;
				strTitle = sprintf('Contrast %.3f; Non-Pref / Resp',dblContrast);
			elseif intPlot == 3
				matTemp = sContrastDataResp(intContrastIndex).matNS;
				strTitle = sprintf('No stim / Resp');
			end
			
			%get data
			vecMeanTrace = mean(matTemp,1);
			vecSE = std(matTemp,[],1)/sqrt(size(matTemp,1));
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			fill(vecX,vecY,[0.7 0.7 1],'EdgeColor','none');
			hold on
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',[0 0 1]);
			hold off
			title(strTitle)
			grid on
			xlabel('Time after response (s)')
			ylabel('dF/F')
			xlim([-3 0])
			ylim([-0.01 0.06])
		end
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_respPEP_c%.1f_raw',strSes,dblContrast*100);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
	end
end

%% act during only stim frames
if sParams.boolDoStimFrames
	cellSaveDuringStim = cell(5,2); %dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values

	%pre-allocate
	intContrasts = length(sMetaDataSF.cellSelectC);
	vecResp = [];
	vecNoResp = [];
	vecMetaRespY = nan(1,intContrasts);
	vecMetaRespE = nan(1,intContrasts);
	vecMetaNoRespY = nan(1,intContrasts);
	vecMetaNoRespE = nan(1,intContrasts);
	intContrastCounter = 0;
	for intContrastIndex=1:length(sMetaDataSF.cellSelectC)
		%get contrast
		dblContrast = sMetaDataSF.vecContrasts(intContrastIndex);
		if dblContrast == 0,continue;end
		intContrastCounter = intContrastCounter + 1;
		
		%figure
		handleFig = figure;
		set(handleFig,'Color',[1 1 1]);
		figure(handleFig);
		
		intPR = length(sContrastDataSF(intContrastIndex).vecPR);
		intNPR = length(sContrastDataSF(intContrastIndex).vecNPR);
		intPN = length(sContrastDataSF(intContrastIndex).vecPN);
		intNPN = length(sContrastDataSF(intContrastIndex).vecNPN);
		intDiffR = length(sContrastDataSF(intContrastIndex).vecDiffR);
		intDiffN = length(sContrastDataSF(intContrastIndex).vecDiffN);
		
		%{
		% preferred
		vecGroupsP = [ones(intPR,1); 2*ones(intPN,1)];
		vecDataP = [sContrastDataSF(intContrastIndex).matPR; sContrastDataSF(intContrastIndex).matPN];
		
		subplot(1,2,1)
		boxplot(vecDataP,vecGroupsP,'notch','on','labels',{'Resp','No resp'})
		ylabel('dF/F during stimulus frames')
		
		strTitle = sprintf('Contrast %.3f Pref',dblContrast);
		title(strTitle);
		%ylim([0 0.11])
		
		%non preferred
		vecDataNP = [sContrastDataSF(intContrastIndex).matNPR; sContrastDataSF(intContrastIndex).matNPN];
		vecGroupsNP = [ones(intNPR,1); 2*ones(intNPN,1)];
		
		subplot(1,2,2)
		boxplot(vecDataNP,vecGroupsNP,'notch','on','labels',{'Resp','No resp'})
		ylabel('dF/F during stimulus frames')
		
		strTitle = sprintf('Contrast %.3f Non-Pref',dblContrast);
		title(strTitle);
		%ylim([0 0.11])
		%}
		
		%difference
		%vecDataD = [sContrastDataSF(intContrastIndex).vecDiffR; sContrastDataSF(intContrastIndex).vecDiffN];
		%vecGroupsD = [ones(intDiffR,1); 2*ones(intDiffN,1)];
		%boxplot(vecDataD,vecGroupsD,'labels',{'Resp','No resp'})
		
		%perform analyses per contrast
		[h,p,ci] = ttest2(sContrastDataSF(intContrastIndex).vecPR,sContrastDataSF(intContrastIndex).vecPN);
		vecY = [mean(sContrastDataSF(intContrastIndex).vecPR) mean(sContrastDataSF(intContrastIndex).vecPN)];
		vecE = [std(sContrastDataSF(intContrastIndex).vecPR)/sqrt(length(sContrastDataSF(intContrastIndex).vecPR)) std(sContrastDataSF(intContrastIndex).vecPN)/sqrt(length(sContrastDataSF(intContrastIndex).vecPN))];
		errorbar(vecY,vecE,'bx')
		set(gca,'XTick',[1 2],'XTickLabel',{'Resp','No resp'})
		ylabel('dF/F during stimulus frames')
		
		strTitle = sprintf('Contrast %.3f Response Pref Pop [p=%.3f]',dblContrast,p);
		title(strTitle);
		%ylim([0 0.11])
		drawnow;
		
		%put in output
		cellSaveDuringStim{intContrastCounter,1} = sContrastDataSF(intContrastIndex).vecPR;
		cellSaveDuringStim{intContrastCounter,2} = sContrastDataSF(intContrastIndex).vecPN;
		
		%put in meta vector
		vecMetaRespY(intContrastIndex) = mean(sContrastDataSF(intContrastIndex).vecPR);
		vecMetaRespE(intContrastIndex) = std(sContrastDataSF(intContrastIndex).vecPR)/sqrt(length(sContrastDataSF(intContrastIndex).vecPR));
		vecMetaNoRespY(intContrastIndex) = mean(sContrastDataSF(intContrastIndex).vecPN);
		vecMetaNoRespE(intContrastIndex) = std(sContrastDataSF(intContrastIndex).vecPN)/sqrt(length(sContrastDataSF(intContrastIndex).vecPN));
		
		%normalize responses per contrast
		if dblContrast ~= 0 && dblContrast ~= 1
			vecNormalized = imnorm([sContrastDataSF(intContrastIndex).vecPR sContrastDataSF(intContrastIndex).vecPN]);
			vecRespThisContrast = vecNormalized(1:length(sContrastDataSF(intContrastIndex).vecPR));
			vecNoRespThisContrast = vecNormalized((length(sContrastDataSF(intContrastIndex).vecPR) + 1):end);
			
			vecResp = [vecResp vecRespThisContrast];
			vecNoResp = [vecNoResp vecNoRespThisContrast];
		end
	end
	
	%group all responses (except highest & lowest contrasts (0 & 100%))
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);
	
	[h,p,ci] = ttest2(vecResp,vecNoResp);
	vecY = [mean(vecResp) mean(vecNoResp)];
	vecE = [std(vecResp)/sqrt(length(vecResp)) std(vecNoResp)/sqrt(length(vecNoResp))];
	errorbar(vecY,vecE,'bx')
	set(gca,'XTick',[1 2],'XTickLabel',{'Resp','No resp'})
	ylabel('Normalized activity during stimulus frames')
	
	strTitle = sprintf('Normalized responses over contrasts for Pref Pop [p=%.3f]',p);
	title(strTitle);
	%ylim([0 0.11])
	
	%plot over contrasts
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);
	
	%pre-compute variables
	vecWindow = [2 6];
	vecWindowSelect = vecWindow(1):vecWindow(end);
	intWL = length(vecWindowSelect);
	vecWindowInv = intWL:-1:1;
	vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
	vecX = [sMetaDataSF.vecContrasts(vecWindowSelect) sMetaDataSF.vecContrasts(vecWindowPlotInv)];
	vecLineX = sMetaDataSF.vecContrasts(vecWindowSelect);
	
	%get data
	for intResp=[0 1]
		if intResp == 1
			vecMeanTrace = vecMetaRespY(vecWindowSelect);
			vecSE = vecMetaRespE(vecWindowSelect);
			vecColorFill = [0.7 0.7 1];
			vecColorLine = [0 0 1];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = vecMetaNoRespY(vecWindowSelect);
			vecSE = vecMetaNoRespE(vecWindowSelect);
		end
		
		
		vecMinTrace = vecMeanTrace-vecSE;
		vecMaxTrace = vecMeanTrace+vecSE;
		vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
		
		%plot
		hold on
		fill(vecX,vecY,vecColorFill,'EdgeColor','none');
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		hold off
		
		
	end
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
	title('Mean +/- st err of population dF/F during stimulus frames')
	grid on
	xlabel('Contrast')
	ylabel('Mean dF/F during stimulus frames')
	xlim(sMetaDataSF.vecContrasts(vecWindow))
	%ylim([-0.01 0.06])
	legend({'StErr','No Resp','StErr','Resp'},'Location','Best')
	
	drawnow;
	if sParams.boolSavePlots
		drawnow;
		strFig = sprintf('%s_stimact_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%% decoding during only stim frames
if sParams.boolDoStimFramesDecoding
	%pre-allocate
	intContrasts = length(sMetaDataSD.cellSelectC);
	vecMetaRespY = nan(1,intContrasts);
	vecMetaNoRespY = nan(1,intContrasts);
	for intContrastIndex=1:length(sMetaDataSD.cellSelectC)
		%get contrast
		dblContrast = sMetaDataSD.vecContrasts(intContrastIndex);
		if dblContrast == 0,continue;end
		
		%perform analyses per contrast
		%chi-square test
		
		%put in meta vector
		vecMetaRespY(intContrastIndex) = nanmean(sContrastDataSD(intContrastIndex).vecDecodedStimTypeR == sContrastDataSD(intContrastIndex).vecStimTypeR);
		vecMetaNoRespY(intContrastIndex) = nanmean(sContrastDataSD(intContrastIndex).vecDecodedStimTypeN == sContrastDataSD(intContrastIndex).vecStimTypeN);
		
		%{
		%normalize responses per contrast
		if dblContrast ~= 0 && dblContrast ~= 1
			vecNormalized = imnorm([sContrastDataPS(intContrastIndex).vecPR sContrastDataPS(intContrastIndex).vecPN]);
			vecRespThisContrast = vecNormalized(1:length(sContrastDataPS(intContrastIndex).vecPR));
			vecNoRespThisContrast = vecNormalized((length(sContrastDataPS(intContrastIndex).vecPR) + 1):end);
			
			vecResp = [vecResp vecRespThisContrast];
			vecNoResp = [vecNoResp vecNoRespThisContrast];
		end
		%}
	end
	
	%{
	%group all responses (except highest & lowest contrasts (0 & 100%))
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);

	[h,p,ci] = ttest2(vecResp,vecNoResp);
	vecY = [mean(vecResp) mean(vecNoResp)];
	vecE = [std(vecResp)/sqrt(length(vecResp)) std(vecNoResp)/sqrt(length(vecNoResp))];
	errorbar(vecY,vecE,'bx')
	set(gca,'XTick',[1 2],'XTickLabel',{'Resp','No resp'})
	ylabel('Normalized activity during stimulus frames')
	
	strTitle = sprintf('Normalized responses over contrasts for Pref Pop [p=%.3f]',p);
	title(strTitle);
	%ylim([0 0.11])
	%}
	
	%plot over contrasts
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);
	
	%pre-compute variables
	vecWindow = [2 6];
	vecWindowSelect = vecWindow(1):vecWindow(end);
	intWL = length(vecWindowSelect);
	vecWindowInv = intWL:-1:1;
	vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
	vecX = [sMetaDataSD.vecContrasts(vecWindowSelect) sMetaDataSD.vecContrasts(vecWindowPlotInv)];
	vecLineX = sMetaDataSD.vecContrasts(vecWindowSelect);
	
	%get data
	for intResp=[0 1]
		if intResp == 1
			vecMeanTrace = vecMetaRespY(vecWindowSelect);
			vecColorFill = [0.7 0.7 1];
			vecColorLine = [0 0 1];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = vecMetaNoRespY(vecWindowSelect);
		end
		
		%plot
		hold on
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		hold off
	end
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
	title('Fraction correct decoded')
	grid on
	xlabel('Contrast')
	ylabel('Fraction correct decoded')
	xlim(sMetaDataSD.vecContrasts(vecWindow))
	%ylim([-0.01 0.06])
	legend({'No Resp','Resp'},'Location','Best')
	
	drawnow;
	if sParams.boolSavePlots
		drawnow;
		strFig = sprintf('%s_decoding_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%% act preceding stim frames
if sParams.boolDoPreStimFrames
	cellSavePreStim = cell(5,2);%dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
	
	%pre-allocate
	intContrasts = length(sMetaDataPS.cellSelectC);
	vecResp = [];
	vecNoResp = [];
	vecMetaRespY = nan(1,intContrasts);
	vecMetaRespE = nan(1,intContrasts);
	vecMetaNoRespY = nan(1,intContrasts);
	vecMetaNoRespE = nan(1,intContrasts);
	intContrastCounter = 0;
	for intContrastIndex=1:length(sMetaDataPS.cellSelectC)
		%get contrast
		dblContrast = sMetaDataPS.vecContrasts(intContrastIndex);
		if dblContrast == 0,continue;end
		intContrastCounter = intContrastCounter + 1;
		
		%figure
		handleFig = figure;
		set(handleFig,'Color',[1 1 1]);
		figure(handleFig);
		
		intPR = length(sContrastDataPS(intContrastIndex).vecPR);
		intNPR = length(sContrastDataPS(intContrastIndex).vecNPR);
		intPN = length(sContrastDataPS(intContrastIndex).vecPN);
		intNPN = length(sContrastDataPS(intContrastIndex).vecNPN);
		
		%perform analyses per contrast
		[h,p,ci] = ttest2(sContrastDataPS(intContrastIndex).vecPR,sContrastDataPS(intContrastIndex).vecPN);
		vecY = [mean(sContrastDataPS(intContrastIndex).vecPR) mean(sContrastDataPS(intContrastIndex).vecPN)];
		vecE = [std(sContrastDataPS(intContrastIndex).vecPR)/sqrt(length(sContrastDataPS(intContrastIndex).vecPR)) std(sContrastDataPS(intContrastIndex).vecPN)/sqrt(length(sContrastDataPS(intContrastIndex).vecPN))];
		errorbar(vecY,vecE,'bx')
		set(gca,'XTick',[1 2],'XTickLabel',{'Resp','No resp'})
		ylabel('dF/F preceding stimulus frames')
		
		strTitle = sprintf('Contrast %.3f Response Pref Pop preceding stimulus [p=%.3f]',dblContrast,p);
		title(strTitle);
		%ylim([0 0.11])
		drawnow;
		
		%put in output
		cellSavePreStim{intContrastCounter,1} = sContrastDataPS(intContrastIndex).vecPR;
		cellSavePreStim{intContrastCounter,2} = sContrastDataPS(intContrastIndex).vecPN;
		
		%put in meta vector
		vecMetaRespY(intContrastIndex) = mean(sContrastDataPS(intContrastIndex).vecPR);
		vecMetaRespE(intContrastIndex) = std(sContrastDataPS(intContrastIndex).vecPR)/sqrt(length(sContrastDataPS(intContrastIndex).vecPR));
		vecMetaNoRespY(intContrastIndex) = mean(sContrastDataPS(intContrastIndex).vecPN);
		vecMetaNoRespE(intContrastIndex) = std(sContrastDataPS(intContrastIndex).vecPN)/sqrt(length(sContrastDataPS(intContrastIndex).vecPN));
		
		%normalize responses per contrast
		if dblContrast ~= 0 && dblContrast ~= 1
			vecNormalized = imnorm([sContrastDataPS(intContrastIndex).vecPR sContrastDataPS(intContrastIndex).vecPN]);
			vecRespThisContrast = vecNormalized(1:length(sContrastDataPS(intContrastIndex).vecPR));
			vecNoRespThisContrast = vecNormalized((length(sContrastDataPS(intContrastIndex).vecPR) + 1):end);
			
			vecResp = [vecResp vecRespThisContrast];
			vecNoResp = [vecNoResp vecNoRespThisContrast];
		end
	end
	
	%group all responses (except highest & lowest contrasts (0 & 100%))
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);
	
	[h,p,ci] = ttest2(vecResp,vecNoResp);
	vecY = [mean(vecResp) mean(vecNoResp)];
	vecE = [std(vecResp)/sqrt(length(vecResp)) std(vecNoResp)/sqrt(length(vecNoResp))];
	errorbar(vecY,vecE,'bx')
	set(gca,'XTick',[1 2],'XTickLabel',{'Resp','No resp'})
	ylabel('Normalized activity preceding stimulus frames')
	
	strTitle = sprintf('Normalized responses over contrasts for Pref Pop preceding stimulus [p=%.3f]',p);
	title(strTitle);
	%ylim([0 0.11])
	
	%plot over contrasts
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);
	
	%pre-compute variables
	vecWindow = [2 6];
	vecWindowSelect = vecWindow(1):vecWindow(end);
	intWL = length(vecWindowSelect);
	vecWindowInv = intWL:-1:1;
	vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
	vecX = [sMetaDataPS.vecContrasts(vecWindowSelect) sMetaDataPS.vecContrasts(vecWindowPlotInv)];
	vecLineX = sMetaDataPS.vecContrasts(vecWindowSelect);
	
	%get data
	for intResp=[0 1]
		if intResp == 1
			vecMeanTrace = vecMetaRespY(vecWindowSelect);
			vecSE = vecMetaRespE(vecWindowSelect);
			vecColorFill = [0.7 0.7 1];
			vecColorLine = [0 0 1];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = vecMetaNoRespY(vecWindowSelect);
			vecSE = vecMetaNoRespE(vecWindowSelect);
		end
		
		
		vecMinTrace = vecMeanTrace-vecSE;
		vecMaxTrace = vecMeanTrace+vecSE;
		vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
		
		%plot
		hold on
		fill(vecX,vecY,vecColorFill,'EdgeColor','none');
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		hold off
		
		
	end
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
	title('Mean +/- st err of population dF/F preceding stimulus frames')
	grid on
	xlabel('Contrast')
	ylabel('Mean dF/F preceding stimulus frames')
	xlim(sMetaDataPS.vecContrasts(vecWindow))
	%ylim([-0.01 0.06])
	legend({'StErr','No Resp','StErr','Resp'},'Location','Best')
	
	drawnow;
	if sParams.boolSavePlots
		drawnow;
		strFig = sprintf('%s_actprecedingstim_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%% alpha power preceding stimulus presentation
if sParams.boolDoPrecedingAlphaPower
	%{
	for intContrastIndex=1:length(sContrastDataAP)
		%plot over contrasts
		handleFig = figure;
		set(handleFig,'Color',[1 1 1]);
		figure(handleFig);
		for intSubPop=1:4
			if intSubPop == 1
				matData = sContrastDataAP(intContrastIndex).matPR;
				strTitle = 'Resp / Pref';
			elseif intSubPop == 2
				matData = sContrastDataAP(intContrastIndex).matNPR;
				strTitle = 'Resp / Non Pref';
			elseif intSubPop == 3
				matData = sContrastDataAP(intContrastIndex).matPN;
				strTitle = 'No Resp / Pref';
			elseif intSubPop == 4
				matData = sContrastDataAP(intContrastIndex).matNPN;
				strTitle = 'No Resp / Non Pref';
			end
			subplot(2,2,intSubPop)
			
			% Plot single-sided amplitude spectrum.
			%pre-compute variables
			dblSamplingFreq = ses.samplingFreq;
			intNFFT = size(matData,2);
			vecLineX = dblSamplingFreq/2*linspace(0,1,intNFFT/2);
			vecInvert = length(vecLineX):-1:1;
			vecX = [vecLineX vecLineX(vecInvert)];
			
			%get data
			matY = 2*abs(matData(:,1:intNFFT/2));
			vecMeanTrace = mean(matY,1);
			vecSE = std(matY,[],1)/sqrt(size(matY,1));
			
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecColorFill = [0.7 0.7 1];
			vecColorLine = [0 0 1];
			
			%plot
			fill(vecX,[vecMinTrace vecMaxTrace(vecInvert)],vecColorFill,'EdgeColor','none');
			hold on
			plot(vecLineX,vecMeanTrace,'LineWidth',2,'Color',vecColorLine)
			hold off
			title(sprintf('Single-Sided Amplitude Spectrum of y(t) for %s',strTitle))
			xlabel('Frequency (Hz)')
			ylabel('|Y(f)|')
			
		end
	end
	%}
	matPower = [];
	for intContrastIndex=1:length(sContrastDataAP)
		matResp = [sContrastDataAP(intContrastIndex).matPR; sContrastDataAP(intContrastIndex).matNPR];
		matNoResp = [sContrastDataAP(intContrastIndex).matPN; sContrastDataAP(intContrastIndex).matNPN];
		
		%plot over contrasts
		handleFig = figure;
		set(handleFig,'Color',[1 1 1]);
		figure(handleFig);
		
		% Plot single-sided amplitude spectrum.
		%pre-compute variables
		dblSamplingFreq = ses.samplingFreq;
		intNFFT = size(matResp,2);
		vecLineX = dblSamplingFreq/2*linspace(0,1,intNFFT/2);
		vecInvert = length(vecLineX):-1:1;
		vecX = [vecLineX vecLineX(vecInvert)];
		
		%get resp data
		matY_R = 2*abs(matResp(:,1:intNFFT/2));
		vecMeanTrace_R = mean(matY_R,1);
		vecSE_R = std(matY_R,[],1)/sqrt(size(matY_R,1));
		
		%get no resp data
		matY_NR = 2*abs(matNoResp(:,1:intNFFT/2));
		vecMeanTrace_NR = mean(matY_NR,1);
		vecSE_NR = std(matY_NR,[],1)/sqrt(size(matY_NR,1));
		
		%diff
		vecPercChange = ((vecMeanTrace_R - vecMeanTrace_NR)./vecMeanTrace_NR)*100;
		vecColorLine = [0 0 1];
		matPower = [matPower;vecPercChange];
			
		%plot
		plot(vecLineX,vecPercChange,'LineWidth',2,'Color',vecColorLine)
		title(sprintf('Contrast %.3f; Single-Sided Amplitude Spectrum of y(t) for %s',sMetaDataAP.vecContrasts(intContrastIndex),strTitle))
		xlabel('Frequency (Hz)')
		ylabel('% higher power preceding response trials')
	end
	
	vecSelectAlpha = 10:15;%7.6091 - 11.8363 Hz
	vecAlphaPowerChange = mean(matPower(2:end,vecSelectAlpha)',1)
	vecSelectDelta = 2:5;%0.8455 - 3.38 Hz
	vecDeltaPowerChange = mean(matPower(2:end,vecSelectDelta)',1)
end

%% act during stim frames with preceding baseline subtraction
if sParams.boolDoStimMinusPreStimFrames
	cellSaveBaseSubtrStimAct = cell(5,2);% dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
	
	%pre-allocate
	intContrasts = length(sMetaDataSB.cellSelectC);
	vecResp = [];
	vecNoResp = [];
	vecMetaRespY = nan(1,intContrasts);
	vecMetaRespE = nan(1,intContrasts);
	vecMetaNoRespY = nan(1,intContrasts);
	vecMetaNoRespE = nan(1,intContrasts);
	intContrastCounter = 0;
	for intContrastIndex=1:length(sMetaDataSB.cellSelectC)
		%get contrast
		dblContrast = sMetaDataSB.vecContrasts(intContrastIndex);
		if dblContrast == 0,continue;end
		intContrastCounter = intContrastCounter + 1;
		
		%figure
		handleFig = figure;
		set(handleFig,'Color',[1 1 1]);
		figure(handleFig);
		
		intPR = length(sContrastDataSB(intContrastIndex).vecPR);
		intNPR = length(sContrastDataSB(intContrastIndex).vecNPR);
		intPN = length(sContrastDataSB(intContrastIndex).vecPN);
		intNPN = length(sContrastDataSB(intContrastIndex).vecNPN);
		
		%perform analyses per contrast
		[h,p,ci] = ttest2(sContrastDataSB(intContrastIndex).vecPR,sContrastDataSB(intContrastIndex).vecPN);
		vecY = [mean(sContrastDataSB(intContrastIndex).vecPR) mean(sContrastDataSB(intContrastIndex).vecPN)];
		vecE = [std(sContrastDataSB(intContrastIndex).vecPR)/sqrt(length(sContrastDataSB(intContrastIndex).vecPR)) std(sContrastDataSB(intContrastIndex).vecPN)/sqrt(length(sContrastDataSB(intContrastIndex).vecPN))];
		errorbar(vecY,vecE,'bx')
		set(gca,'XTick',[1 2],'XTickLabel',{'Resp','No resp'})
		ylabel('dF/F preceding stimulus frames')
		
		strTitle = sprintf('Contrast %.3f Response Pref Pop preceding stimulus [p=%.3f]',dblContrast,p);
		title(strTitle);
		%ylim([0 0.11])
		drawnow;
		
		%put in output
		cellSaveBaseSubtrStimAct{intContrastCounter,1} = sContrastDataSB(intContrastIndex).vecPR;
		cellSaveBaseSubtrStimAct{intContrastCounter,2} = sContrastDataSB(intContrastIndex).vecPN;
		
		%put in meta vector
		vecMetaRespY(intContrastIndex) = mean(sContrastDataSB(intContrastIndex).vecPR);
		vecMetaRespE(intContrastIndex) = std(sContrastDataSB(intContrastIndex).vecPR)/sqrt(length(sContrastDataSB(intContrastIndex).vecPR));
		vecMetaNoRespY(intContrastIndex) = mean(sContrastDataSB(intContrastIndex).vecPN);
		vecMetaNoRespE(intContrastIndex) = std(sContrastDataSB(intContrastIndex).vecPN)/sqrt(length(sContrastDataSB(intContrastIndex).vecPN));
		
		%normalize responses per contrast
		if dblContrast ~= 0 && dblContrast ~= 1
			vecNormalized = imnorm([sContrastDataSB(intContrastIndex).vecPR sContrastDataSB(intContrastIndex).vecPN]);
			vecRespThisContrast = vecNormalized(1:length(sContrastDataSB(intContrastIndex).vecPR));
			vecNoRespThisContrast = vecNormalized((length(sContrastDataSB(intContrastIndex).vecPR) + 1):end);
			
			vecResp = [vecResp vecRespThisContrast];
			vecNoResp = [vecNoResp vecNoRespThisContrast];
		end
	end
	
	%group all responses (except highest & lowest contrasts (0 & 100%))
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);
	
	[h,p,ci] = ttest2(vecResp,vecNoResp);
	vecY = [mean(vecResp) mean(vecNoResp)];
	vecE = [std(vecResp)/sqrt(length(vecResp)) std(vecNoResp)/sqrt(length(vecNoResp))];
	errorbar(vecY,vecE,'bx')
	set(gca,'XTick',[1 2],'XTickLabel',{'Resp','No resp'})
	ylabel('Normalized baseline-subtracted activity during stimulus frames')
	
	strTitle = sprintf('Normalized responses over contrasts for Pref Pop [p=%.3f]',p);
	title(strTitle);
	%ylim([0 0.11])
	
	%plot over contrasts
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);
	
	%pre-compute variables
	vecWindow = [2 6];
	vecWindowSelect = vecWindow(1):vecWindow(end);
	intWL = length(vecWindowSelect);
	vecWindowInv = intWL:-1:1;
	vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
	vecX = [sMetaDataSB.vecContrasts(vecWindowSelect) sMetaDataSB.vecContrasts(vecWindowPlotInv)];
	vecLineX = sMetaDataSB.vecContrasts(vecWindowSelect);
	
	%get data
	for intResp=[0 1]
		if intResp == 1
			vecMeanTrace = vecMetaRespY(vecWindowSelect);
			vecSE = vecMetaRespE(vecWindowSelect);
			vecColorFill = [0.7 0.7 1];
			vecColorLine = [0 0 1];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = vecMetaNoRespY(vecWindowSelect);
			vecSE = vecMetaNoRespE(vecWindowSelect);
		end
		
		
		vecMinTrace = vecMeanTrace-vecSE;
		vecMaxTrace = vecMeanTrace+vecSE;
		vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
		
		%plot
		hold on
		fill(vecX,vecY,vecColorFill,'EdgeColor','none');
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		hold off
		
		
	end
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
	title('Mean +/- st err of population baseline-subtracted dF/F')
	grid on
	xlabel('Contrast')
	ylabel('Mean baseline-subtracted dF/F')
	xlim(sMetaDataSB.vecContrasts(vecWindow))
	%ylim([-0.01 0.06])
	legend({'StErr','No Resp','StErr','Resp'},'Location','Best')
	
	drawnow;
	
	if sParams.boolSavePlots
		drawnow;
		strFig = sprintf('%s_normact_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end


%% save data structures
vecClock = fix(clock);
strDate = num2str(vecClock(1:3));
strRecs = num2str(vecRecordings);
strFile = ['data_aggregate' strSes '_' strrep(strDate,'    ','_') '_' strRecs(~isspace(strRecs))];
%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
	load(['D:\Data\Results\stimdetection\' strFile]);
end
save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');

%save workspace
clear cellSes sLoad sSesAggregate ses
save(['D:\Data\Results\stimdetection\dataRawPre_aggregate' strSes '.mat'],'-v7.3');
cd(strOldDir);

%{
figures to save and take mean over animals
- cellSavePreStim: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveBehavDetect: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
- cellSaveBehavRT: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
- cellSaveBaseSubtrStimAct: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveDuringStim: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveTraceAct: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
- AD: cellSaveNormActDissim: within-group z-scored activation dissimilarity [2 (high/low DCAE) x 2 (detect/no-detect)] x n (animals) with every cell = vector of dissimilarity values
- AD: cellSaveSignalCorrs: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals)] with every cell = vector of pairwise correlation values
- AD: cellSaveNoiseCorrs: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals)] with every cell = vector of pairwise correlation values
- AD: cellSaveDissimCorrITD: inter-trial-distance dependence of assembly consistency [2 (high/low DCAE) x 2 (hits/misses) x n (animals)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	- extra graph of above with means comparison of hits/misses high/low DCAE
- AD: cellSaveDCAE: normalized increase in dF/F for detection trials [5 (c) x n (animals)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials

%}