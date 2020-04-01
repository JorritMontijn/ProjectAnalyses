%%{
%define target dir
clear all;
boolClust=0;
if boolClust
	strFile = mfilename;
	strHome = mfilename('fullpath');
	strHome = strHome(1:(end-length(mfilename)));
	strHome = strHome(1:(end-length('Scripts\')));
	strOutputDir = [strHome 'Data' filesep];
	strDataDir = [strHome 'SessionData' filesep];
	%add recursive subdirectories
	cellPaths = getSubDirs(strHome,inf,{'old','backup'});
	fprintf('path 1: %s; nr of paths=%d\n',cellPaths{1},numel(cellPaths));
	for intPath=1:numel(cellPaths)
		addpath(cellPaths{intPath});
	end
else
	strDataDir = 'D:\Data\Guido\';
	strOutputDir = 'D:\Data\Results\correlationAnalysis\';
end

%starting message
clc;
fprintf('\n		WELCOME TO CONDOR :: THE MIND DOMINATOR :: CORRELATION IS FUTILE\n\n')
fprintf('Process started at [%s] %s. Please be patient.\n\n',getTime,date);

%load session file database
sLoad = load([strDataDir 'SessionDatabase.mat']);
cellRecGratings = sLoad.cellRecGratings;

%get which mouse/mice to process
strPath = mfilename('fullpath');
strPath = strPath(1:(end-length(mfilename)));
sFiles = dir([strPath '*.queueOrder']);
cellSplit = strparse(sFiles(1).name, '_');
vecMice = cellfun(@str2double,cellSplit(2:(end-1)));
boolSavePlots = true;

for intMouse=[1:5 vecMice]
	%mouse name
	strMouse = num2str(intMouse);
	fprintf('Retrieving data for mouse %s [%s]\n',strMouse,getTime);
	
	%get data
	cellRec = cellRecGratings{intMouse};
	[cellSes, cellSesInclude, indInclude] = getIncludeNeurons( strDataDir, cellRec);
	
	%build aggregate
	sSesAggregate = [];
	for intSes=1:length(cellSes)
		%load session file
		ses = cellSes{intSes};
		
		%check for SecsOn/SecsOff
		if ~isfield(ses.structStim,'SecsOn')
			ses.structStim.SecsOff = ses.structStim.FrameOff * ses.samplingFreq;
			ses.structStim.SecsOn = ses.structStim.FrameOn * ses.samplingFreq;
		end
		
		%create aggregate
		sSesAggregate = buildMultiSesAggregate(ses,sSesAggregate);
		clear ses;
	end
	fprintf('Reformatting data [%s]\n',getTime);
	dblFraction = 0.5;
	dblSecWindowSize = 30;
	%[sSesAggregate,indKeepList] = doRecalcdFoF(sSesAggregate,3,[],'neuron',dblFraction,dblSecWindowSize);
	fprintf('Performing preparatory calculations [%s]\n',getTime);
	
	%get info
	sTypes = getStimulusTypes(sSesAggregate);
	cellSelect = getSelectionVectors(sSesAggregate.structStim,sTypes);
	intTrials = length(sSesAggregate.structStim.Orientation);
	intNeurons = numel(sSesAggregate.neuron);
	intStimTypes = length(cellSelect);
	
	%get noise correlations
	structOut = calcStimCorrs(sSesAggregate);
	matAllNoise = structOut.matAllNoise;
	
	%get response matrix
	matResp = getNeuronResponse(sSesAggregate);
	varSave_matResp = matResp;
	vecClasses=sSesAggregate.structStim.Orientation';
	matData = matResp';
	vecOris = unique(sSesAggregate.structStim.Orientation);
	
	%% analyze how much information is present (on average) in a single group of neurons with the specified size
	boolPlot = false;
	
	%get ori indices & trial index
	vecTrialOriIndex = sum((repmat(vecClasses,[1 intStimTypes])==repmat(vecOris,[intTrials 1])).*repmat(1:intStimTypes,[intTrials 1]),2);

	%% calculate multidimensional orthogonal/parallel distance
	vecDimensionalities = 2:intNeurons;
	if vecDimensionalities(end) ~= intNeurons,vecDimensionalities(end+1) = intNeurons;end %#ok<SAGROW>
	vecLookupDimensionalities = vecDimensionalities;
	intDims = length(vecDimensionalities);
	intMaxSampleSize=1000;
	
	matOrthSD = nan(intStimTypes,intDims,intMaxSampleSize);
	matParaSD = nan(intStimTypes,intDims,intMaxSampleSize);
	matOrthSD_Shuf = nan(intStimTypes,intDims,intMaxSampleSize);
	matParaSD_Shuf = nan(intStimTypes,intDims,intMaxSampleSize);
	
	for intDimensionality = vecDimensionalities
		%% start msg
		intDimIndex = find(intDimensionality == vecDimensionalities,1);
		intCombinations = round(factorial(intNeurons)/(factorial(intDimensionality)*factorial(intNeurons-intDimensionality)));
		intSampleSize = min([intMaxSampleSize intCombinations]);
		fprintf('Orthogonal/parallel variability; Dimensionality %d with %d cells... Using mean of %d samples [%s]\n',intDimensionality,intNeurons,intSampleSize,getTime);
		
		%check if we can use nchoosek
		if intCombinations <= intMaxSampleSize
			matNeuronCombs = nchoosek(1:intNeurons,intDimensionality);
		else
			%otherwise, create random combinations
			matNeuronCombs = ones(intSampleSize,intDimensionality);
			for intRandCombo=1:intSampleSize
				vecRand = sort(randperm(intNeurons,intDimensionality));
				
				%check if this combination already exists
				intOverflow = 0;
				while any(all(bsxfun(@eq,matNeuronCombs,vecRand),2)) && intOverflow < 10000
					vecRand = sort(randperm(intNeurons,intDimensionality));
					intOverflow = intOverflow + 1;
				end
				%assign last combination
				matNeuronCombs(intRandCombo,:) = vecRand;
			end
		end
		
		%perform random estimations
		for intComboCounter=1:intSampleSize
			%get neuron group
			vecNeuronGroup = matNeuronCombs(intComboCounter,:);
			matDataPoints = matData(:,vecNeuronGroup);
			
			%calc orth/para
			[vecSD_Orth,vecSD_Para] = calcOrthPara(matDataPoints,vecTrialOriIndex);
			[vecShuffSD_Orth,vecShuffSD_Para] = calcOrthPara(matDataPoints,vecTrialOriIndex,true);
			
			%assign
			matOrthSD(:,intDimIndex,intComboCounter) = vecSD_Orth;
			matParaSD(:,intDimIndex,intComboCounter) = vecSD_Para;
			matOrthSD_Shuf(:,intDimIndex,intComboCounter) = vecShuffSD_Orth;
			matParaSD_Shuf(:,intDimIndex,intComboCounter) = vecShuffSD_Para;
		end
		%relative
		matSDRel_Orth = (nanmean(nanmean(matOrthSD,3),1)./nanmean(nanmean(matOrthSD_Shuf,3),1))*100;
		matSDRel_Para = (nanmean(nanmean(matParaSD,3),1)./nanmean(nanmean(matParaSD_Shuf,3),1))*100;
	end
	
	%%
	close all;
	figure
	subplot(2,2,1)
	plot(vecDimensionalities,nanmean(nanmean(matOrthSD,3),1),'r');
	hold on
	plot(vecDimensionalities,nanmean(nanmean(matParaSD,3),1),'b');
	hold off
	title('Raw data; red=orth; blue=para')
	xlabel('Dimensionality')
	ylabel('Variability (dF/F0)')
	
	subplot(2,2,2)
	plot(vecDimensionalities,nanmean(nanmean(matOrthSD_Shuf,3),1),'r');
	hold on
	plot(vecDimensionalities,nanmean(nanmean(matParaSD_Shuf,3),1),'b');
	hold off
	title('Shuffled data')
	xlabel('Dimensionality')
	ylabel('Variability (dF/F0)')
	
	subplot(2,2,3)
	%relative
	vecSDRel_Orth = (nanmean(nanmean(matOrthSD,3),1)./nanmean(nanmean(matOrthSD_Shuf,3),1))*100;
	vecSDRel_Para = (nanmean(nanmean(matParaSD,3),1)./nanmean(nanmean(matParaSD_Shuf,3),1))*100;
	plot(vecDimensionalities,vecSDRel_Orth,'r')
	hold on
	plot(vecDimensionalities,vecSDRel_Para,'b')
	hold off
	drawnow;
	title('Data normalized to shuffled')
	xlabel('Dimensionality')
	ylabel('Normalized varability (%)')
	
	if boolSavePlots
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('mouse%s_HDV_ori_raw',strMouse);
		export_fig([strOutputDir 'mouse' strMouse filesep strFig '.tif']);
		export_fig([strOutputDir 'mouse' strMouse filesep strFig '.pdf']);
	end
	
	%% save data
	vecClock = fix(clock);
	strDate = num2str(vecClock(1:3));
	strFile = ['LDA_HDV' strMouse '_' strrep(strDate,'    ','_')];
	save([strOutputDir strFile],'matOrthSD','matParaSD','matOrthSD_Shuf','matParaSD_Shuf','vecSDRel_Orth','vecSDRel_Para','vecDimensionalities');
end




