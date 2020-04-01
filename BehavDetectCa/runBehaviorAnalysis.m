%function [ output_args ] = runGetStimResponse
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Currently uses Orientation and Contrast fields

clear all;
close all;
%% set which analyses to run
sParams = struct;
boolDoStimCorrPlot = true;
boolDoStimDecoding = false;
boolDoStimDecodingBehavior = true;
boolDoTraceDecoding = false;
boolDoTraceBehaviorDecoding = false;
boolDoTracePresenceDecoding = false;
boolDoTracePEP = true; %peri-event plot
boolDoTraceBehaviorPEP = true;
boolDoTraceStimBehaviorPEP = true;
boolDoStimDecodingPEP = false;
boolDoCorrPECM = false; %peri-event correlation matrix %TODO
if isfield(sParams,'boolCombineOppositeStimuli'), boolCombineOppositeStimuli = sParams.boolCombineOppositeStimuli; else boolCombineOppositeStimuli = true;end
	
intDataType = 3;%0=s-np;3=s;4=np

%% set which recordings to process
intProcess=1;
if intProcess == 1
	%behavior
	%cellSes{1} = '20130627';
	cellSes{1} = '20140207';
	%cellRec{1} = 'xyt01';
	cellRec{1} = 'xyt01';
elseif intProcess == 2
	%dir tuning
	%np: 60, 39, 50, 70, 80, 66
	%s: 98, 75, 90, 96, 99, 89
	%s-np: 100, 98, 100, 99, 100, 100
	cellSes{1}= '20120718';
	cellRec{1}= 'xyt01';
	
	cellSes{2}= '20120718';
	cellRec{2}= 'xyt03';
	
	cellSes{3}= '20120720';
	cellRec{3}= 'xyt01';
	
	cellSes{4}= '20120720';
	cellRec{4}= 'xyt03';
	
	cellSes{5}= '20121207';
	cellRec{5}= 'xyt01';
	
	cellSes{6}= '20121207';
	cellRec{6}= 'xyt02';
	
elseif intProcess == 3
	%ori tuning
	%np: 24, 29
	%s: 80, 63
	%s-np: 99, 95
	
	cellSes{1}= '20130612';
	cellRec{1}= 'xyt02';
	
	cellSes{2}= '20130625';
	cellRec{2}= 'xyt02';
elseif intProcess == 4
	%plaids
	cellSes{1}= '20130612';
	cellRec{1}= 'xyt01';
	
	cellSes{2}= '20130625';
	cellRec{2}= 'xyt01';
elseif intProcess == 5
	%old dir/plaids
	cellSes{1}= '20130313';
	cellRec{1}= 'xyt02';
	
	cellSes{2}= '20130313';
	cellRec{2}= 'xyt03';
end

%% start processing
for intRec=1:1%numel(cellSes)
	%% load data
	clear('ses','sParamsSD','sOutDB','sOutTD','sTC','sBR')
	strSes = cellSes{intRec};
	strRec = cellRec{intRec};
	
	fprintf('\n   >>> Loading and preparing data file for session %s%s\n\n',strSes,strRec);
	
	strDir = ['D:' filesep 'Data' filesep 'Processed' filesep 'imagingdata' filesep strSes filesep strRec filesep];
	strFile = sprintf('%s%s_ses.mat',strSes,strRec);
	
	load([strDir strFile]);
	
	
	%combine opposite directions 
	if boolCombineOppositeStimuli
		ses.structStim.Orientation = mod(ses.structStim.Orientation,180);
	end
	
	%transform pulses to frames
	cellFields = fields(ses.structStim);
	for intField=1:length(cellFields)
		strField = cellFields{intField};
		strNewField = strrep(strField,'Pulse','Frame');
		ses.structStim.(strNewField) = ses.structStim.(strField);
	end
	
	if isfield(ses.structStim,'cellRespFrames')
		boolBehavior = true;
	else
		boolBehavior = false;
	end
	
	if intDataType == 0
		strDataType='snp';
	elseif intDataType == 3
		strDataType='soma';
		ses = doRectifydFoF(ses,3); %none=s-np;3=s;4=np
	elseif intDataType == 4
		strDataType='np';
		ses = doRectifydFoF(ses,4); %none=s-np;3=s;4=np
	end
	%ses = doRectifydFoF(ses,2);
	
	
	%% compare mouse response to decoding performance
	if boolDoStimDecodingBehavior && boolBehavior
		% get indexing vectors for unique stimulus combinations
		sTypes = getStimulusTypes(ses);
		cellSelect = getSelectionVectors(ses.structStim,sTypes);
		
		%remove empty types
		intTypes = length(cellSelect);
		vecKeep = true(1,intTypes);
		for intType = 1:intTypes
			if sum(cellSelect{intType}) == 0
				vecKeep(intType) = false;
			end
		end
		cellSelect = cellSelect(vecKeep);
		sTypes = getStimulusTypes(ses);
		sTypes.matTypes = sTypes.matTypes(:,vecKeep);
		
		%bootstrap decoding
		sParamsDB.sTypes = sTypes;
		sParamsDB.verbose = 0;
		sParamsDB.cellSelect = cellSelect;
		intNeurons = numel(ses.neuron);
		intStims = numel(ses.structStim.TrialNumber);
		intIters = 100;
		matCorrect = false(intIters,intStims);
		intSampleSize = 10;%round(intNeurons/2);
		for intIter=1:intIters
			%random sample
			vecR = randsample(1:intNeurons,intSampleSize);
			sParamsDB.vecIncludeCells = vecR;
			
			%decode
			sOutDB = doStimDecoding(ses,sParamsDB);
			
			%save output
			vecCorrect = sOutDB.vecDecodedStimType == sOutDB.vecStimType;
			matCorrect(intIter,:) = vecCorrect;
			
			%msg
			if mod(intIter,round(intIters/10)) == 0
				fprintf('Now at bootstrap iteration %d of %d\n',intIter,intIters)
			end
		end
		
		%plot
		vecDecodingPercentage = (sum(matCorrect,1)/size(matCorrect,1))*100;
		figure;plot(vecDecodingPercentage);
		title(sprintf('Bootstrapped decoding performance, iters=%d, sample size=%d/%d',intIters,intSampleSize,intNeurons))
		ylim([0 100])
		xlim([1 intStims])
		xlabel('Stimulus')
		ylabel('Percentage correct')
		set(gcf,'color','w');
		export_fig(sprintf('stim_decoding_bootstrapped%s%s.tif',strSes,strRec))
		
		% plot decoding performance vs stimulus duration split by correct/incorrect detection
		vecStimDur = (ses.structStim.FrameOff - ses.structStim.FrameOn)/ses.samplingFreq;
		vecResp = ses.structStim.vecTrialResponse;
		
		figure;
		plot(vecStimDur(vecResp==1),vecDecodingPercentage(vecResp==1),'bo')
		hold on
		plot(vecStimDur(vecResp==0),vecDecodingPercentage(vecResp==0),'rx')
		hold off
		title(sprintf('Decoding performance vs stim duration; %s%s',strSes,strRec))
		ylim([0 100])
		xlabel('Stimulus duration (secs)')
		ylabel('Decoding percentage correct')
		legend('Correct Response','Incorrect Response','Location','Best')
		set(gcf,'color','w');
		export_fig(sprintf('stim_decoding_duration_detection%s%s.tif',strSes,strRec))
	end
	
	%% calculate correlation trace / prep trace decoding
	if boolDoTraceDecoding || boolDoTracePresenceDecoding || boolDoTracePEP || boolDoStimDecodingPEP || boolDoTraceBehaviorPEP || boolDoTraceStimBehaviorPEP || boolDoTraceBehaviorDecoding
		%define window size
		intWindowLength = round(ses.samplingFreq);
		intFrames = length(ses.neuron(1).dFoF);
		
		%set input parameters
		sParamsTC.vecEpochStart = 1:(intFrames-intWindowLength+1);
		sParamsTC.intEpochDuration = intWindowLength;
		
		%msg
		ptrTime = tic;
		fprintf('\nCalculating sliding window correlations for entire trace [%d Frames]. This might take a while. Please be patient...\n',intFrames)
		
		%get correlation trace
		sTC = processActivityMatrix(ses,sParamsTC);
		matDistances = sTC.matDistances;
		
		%msg
		fprintf('Done! Calculation took %.1f seconds\n',toc(ptrTime))
		
		% get indexing vectors for unique stimulus combinations
		sTypes = getStimulusTypes(ses);
		cellSelect = getSelectionVectors(ses.structStim,sTypes);
		
		%remove empty types
		intTypes = length(cellSelect);
		vecKeep = true(1,intTypes);
		for intType = 1:intTypes
			if sum(cellSelect{intType}) == 0
				vecKeep(intType) = false;
			end
		end
		cellSelect = cellSelect(vecKeep);
		sTypes = getStimulusTypes(ses);
		sTypes.matTypes = sTypes.matTypes(:,vecKeep);
		
		%set params
		sParamsTD.sTypes = sTypes;
		sParamsTD.cellSelect = cellSelect;
		sParamsTD.intWindowLength = intWindowLength;
	end
	
	%% decode stimulus types continuously and compare with behavior and correlation
	if boolDoTraceDecoding
		%do decoding
		sOutTD = doTraceDecoding(ses,sParamsTD);
		
		%get behavioral data
		if boolBehavior
			%get responses
			vecRespFrames = [];
			vecRespTypes = [];
			for intTrial=1:numel(ses.structStim.cellRespFrames)
				vecSelect = ~isnan(ses.structStim.cellRespFrames{intTrial});
				vecRespFrames = [vecRespFrames ses.structStim.cellRespFrames{intTrial}(vecSelect)];
				vecRespTypes = [vecRespTypes ses.structStim.cellResponse{intTrial}(vecSelect)];
			end
			sBR.vecRespFrames = vecRespFrames;
			sBR.vecRespTypes = vecRespTypes;
			
			%get movement
			if isfield(ses.structStim,'cellMovement')
				
			end
		else
			sBR = [];
		end
		
		%plot
		doPlotTraceDecoding(sOutTD,sTC,sBR);
	end
	
	%% decode behavioral responses
	if boolDoTraceBehaviorDecoding
		%set params
		sParamsTBD.sTypes = sTypes;
		sParamsTBD.cellSelect = cellSelect;
		sParamsTBD.intWindowLength = intWindowLength;
		
		%do decoding
		sOutTBD = doTraceDecoding(ses,sParamsTBD);
		
		%plot
		doPlotTraceDecoding(sOutTBD,sTC,sBR);
	end
	
	%% PSTH of stimuli/continuous decoding/mean pop correlation for stimuli
	if boolDoTracePEP
		%get dF/F trace
		vecMeanAct=mean(sTC.matEpoch,1);
		vecSdAct=std(sTC.matEpoch,[],1);
		
		%get correlation trace
		intNeurons = size(sTC.matEpoch,1);
		intFrames = length(vecMeanAct);
		matSelect = tril(true(intNeurons,intNeurons),-1);
		vecCorrMean = nan(1,intFrames);
		vecCorrSD = nan(1,intFrames);
		for intFrame=1:size(sTC.matRawCovar,3)
			matCorr = sTC.matRawCovar(:,:,intFrame);
			vecCorrMean(intFrame+floor(intWindowLength/2)) = mean(matCorr(matSelect));
			vecCorrSD(intFrame+floor(intWindowLength/2)) = std(matCorr(matSelect));
		end
		
		%set plotting parameters
		sEvents = struct;
		sEvents.vecOn = ses.structStim.FrameOn;
		sEvents.vecOff = ses.structStim.FrameOff;
		%sEvents.cellSelect = sParamsTD.cellSelect;
		%sEvents.sTypes = sParamsTD.sTypes;
		sEvents.vecWindow = round([-3*ses.samplingFreq 5*ses.samplingFreq]);
		sEvents.dblFrameRate = ses.samplingFreq;
		
		%plot mean corr
		doPEP(sEvents,vecCorrMean);
		ylabel('Mean population correlation')
		title(sprintf('Stimulus PEP %s%s',strSes,strRec));
		
		%plot std corr
		doPEP(sEvents,vecMeanAct);
		ylabel('Mean population dF/F')
		title(sprintf('Stimulus PEP %s%s',strSes,strRec));
		
		%get pref stim per neuron
		structStimCorrs = calcStimCorrs(ses);
		matSignalCorrs = structStimCorrs.matSignalCorrs;
		matNoiseCorrs = structStimCorrs.matNoiseCorrs;
		matStimResponse = structStimCorrs.matStimResponse;
		matSignalResponse = structStimCorrs.matSignalResponse;
		[vecStimAct,vecPrefStim]=max(matSignalResponse,[],1);
		
		%plot per type
		for intStimType=1:length(sParamsTD.cellSelect)
			vecSelectNeurons = vecPrefStim==intStimType;
			intSelectNeurons = sum(vecSelectNeurons);
			if intSelectNeurons>0
				vecSelectTrials = sParamsTD.cellSelect{intStimType};
				
				%set plotting parameters
				sEvents = struct;
				sEvents.handleFig=figure;
				sEvents.vecOn = ses.structStim.FrameOn(vecSelectTrials);
				sEvents.vecOff = ses.structStim.FrameOff(vecSelectTrials);
				sEvents.vecWindow = round([-3*ses.samplingFreq 5*ses.samplingFreq]);
				sEvents.dblFrameRate = ses.samplingFreq;
				set(gcf,'Color',[1 1 1]);
				
				%get data
				%get dF/F trace
				vecMeanActPref=mean(sTC.matEpoch(vecSelectNeurons,:),1);
				vecSdActPref=std(sTC.matEpoch(vecSelectNeurons,:),[],1);
				vecMeanActNonPref=mean(sTC.matEpoch(~vecSelectNeurons,:),1);
				vecSdActNonPref=std(sTC.matEpoch(~vecSelectNeurons,:),[],1);
				
				%plot mean act pref
				subplot(2,1,1)
				doPEP(sEvents,vecMeanActPref);
				ylabel('Mean pref population dF/F0')
				title(sprintf('Stimulus %s PEP %s%s',num2str(sParamsTD.sTypes.matTypes(end,intStimType)),strSes,strRec));
				
				%plot mean act non-pref
				subplot(2,1,2)
				doPEP(sEvents,vecMeanActNonPref);
				ylabel('Mean non-pref population dF/F0')
				title(sprintf('Stimulus %s PEP %s%s',num2str(sParamsTD.sTypes.matTypes(end,intStimType)),strSes,strRec));
				%{
				%get correlation trace
				intNeurons = size(sTC.matEpoch,1);
				intFrames = size(sTC.matRawCovar,3);
				matSelect = tril(true(intNeurons,intNeurons),-1);
				matSelect(~vecSelectNeurons,:) = false;
				vecCorrMean = zeros(1,intFrames);
				vecCorrSD = zeros(1,intFrames);
				for intFrame=1:intFrames
					matCorr = sTC.matRawCovar(:,:,intFrame);
					vecCorrMean(intFrame) = mean(matCorr(matSelect));
					vecCorrSD(intFrame) = std(matCorr(matSelect));
				end
		
		


				%plot std corr
				doPEP(sEvents,vecMeanActPref);
				%ylabel('Mean population dF/F')
				title(sprintf('Stimulus PEP %s%s',strSes,strRec));
				%}
			end
		end
	end
	
	%% PSTH of continuous stim decoding for stims
	if boolDoStimDecodingPEP
		if ~boolDoTraceDecoding
			%do decoding
			sOutTD = doTraceDecoding(ses,sParamsTD);
		end
		
		%calc decoding trace
		vecStimTrace = getStimTrace(ses,2);
		vecDecodingTrace = sOutTD.vecDecodedType == vecStimTrace;
		
		%set plotting parameters
		sEvents = struct;
		sEvents.vecOn = ses.structStim.FrameOn;
		sEvents.vecOff = ses.structStim.FrameOff;
		%sEvents.cellSelect = sParamsTD.cellSelect;
		%sEvents.sTypes = sParamsTD.sTypes;
		sEvents.vecWindow = round([-3*ses.samplingFreq 5*ses.samplingFreq]);
		sEvents.dblFrameRate = ses.samplingFreq;
		
		%plot mean corr
		doPEP(sEvents,vecDecodingTrace);
		ylabel('Fraction stimulus decoded')
		ylim([0 1])
		title(sprintf('Stimulus PEP %s%s',strSes,strRec));
	end
	
	%% PSTH of continuous decoding/mean pop correlation/responses for response times
	if boolDoTraceBehaviorPEP
		
		%get behavioral data
		if ~exist('sBR','var') || isempty(sBR)
			vecRespFrames = [];
			vecRespTypes = [];
			for intTrial=1:numel(ses.structStim.cellRespFrames)
				vecSelect = ~isnan(ses.structStim.cellRespFrames{intTrial});
				vecRespFrames = [vecRespFrames ses.structStim.cellRespFrames{intTrial}(vecSelect)];
				vecRespTypes = [vecRespTypes ses.structStim.cellResponse{intTrial}(vecSelect)];
			end
			sBR.vecRespFrames = vecRespFrames;
			sBR.vecRespTypes = vecRespTypes;
		else
			vecRespFrames = sBR.vecRespFrames;
			vecRespTypes = sBR.vecRespTypes;
		end
		intRespTypes = length(getUniqueVals(sBR.vecRespTypes));
		
		
		%get good resps
		vecRespFramesCorrect = ses.structStim.vecTrialRespFrames(logical(ses.structStim.vecTrialResponse));
		vecRespTypes(ismember(vecRespFrames,vecRespFramesCorrect)) = 2;
		
		%set plotting parameters
		sEvents = struct;
		sEvents.vecOn = vecRespFrames;
		sEvents.vecTypes = vecRespTypes;
		sEvents.vecWindow = round([-3*ses.samplingFreq 5*ses.samplingFreq]);
		sEvents.dblFrameRate = ses.samplingFreq;
		
		%plot mean corr
		cellHandles = doPEP(sEvents,vecCorrMean);
		cmapResp = colormap(jet(intRespTypes));
		colormap('default');
		figure(cellHandles{1});
		ylabel('Mean population correlation');
		title(sprintf('Incorrect Response PEP %s%s',strSes,strRec));
		figure(cellHandles{2});
		ylabel('Mean population correlation');
		title(sprintf('Correct Response PEP %s%s',strSes,strRec));
		
		%plot std corr
		vecHandles = doPEP(sEvents,vecMeanActPref);
		figure(vecHandles{1});
		ylabel('Mean population dF/F');
		title(sprintf('Incorrect Response PEP %s%s',strSes,strRec));
		figure(vecHandles{2});
		ylabel('Mean population dF/F');
		title(sprintf('Correct Response PEP %s%s',strSes,strRec));
	end
	
	%% PSTH of continuous decoding/mean pop correlation/responses for response times separated per stim type
	if boolDoTraceStimBehaviorPEP
		
		%get pref stim per neuron
		structStimCorrs = calcStimCorrs(ses);
		matSignalCorrs = structStimCorrs.matSignalCorrs;
		matNoiseCorrs = structStimCorrs.matNoiseCorrs;
		matStimResponse = structStimCorrs.matStimResponse;
		matSignalResponse = structStimCorrs.matSignalResponse;
		[vecStimAct,vecPrefStim]=max(matSignalResponse,[],1);
		
		%prep pref-non-pref figs
		intTypes = length(sParamsTD.cellSelect);
		handleFigPNP = figure;
		set(handleFigPNP,'Color',[1 1 1]);
		intNumSubX=ceil(sqrt(intTypes));
		intNumSubY=ceil(sqrt(intTypes));
		
		%plot per type
		for intStimType=1:intTypes
			vecSelectNeurons = vecPrefStim==intStimType;
			intSelectNeurons = sum(vecSelectNeurons);
			if intSelectNeurons>0
				%figure
				handleFig = figure;
				set(handleFig,'Color',[1 1 1]);
				figure(handleFig);
				
				%set plotting parameters for correct trials
				vecCorrResp = logical(ses.structStim.vecTrialResponse);
				vecSelectTrials = sParamsTD.cellSelect{intStimType} & vecCorrResp;
				sEvents = struct;
				sEvents.handleFig=handleFig;
				sEvents.vecOn = ses.structStim.FrameOn(vecSelectTrials);
				sEvents.vecOff = ses.structStim.FrameOff(vecSelectTrials);
				sEvents.vecWindow = round([-3*ses.samplingFreq 5*ses.samplingFreq]);
				sEvents.dblFrameRate = ses.samplingFreq;
				
				%get data
				vecMeanActPref=mean(sTC.matEpoch(vecSelectNeurons,:),1);
				vecSdActPref=std(sTC.matEpoch(vecSelectNeurons,:),[],1);
				vecMeanActNonPref=mean(sTC.matEpoch(~vecSelectNeurons,:),1);
				vecSdActNonPref=std(sTC.matEpoch(~vecSelectNeurons,:),[],1);
				
				%plot mean act pref cor
				subplot(2,2,1)
				doPEP(sEvents,vecMeanActPref);
				ylabel('Mean pref population dF/F0 correct')
				title(sprintf('Stimulus %s PEP %s%s',num2str(sParamsTD.sTypes.matTypes(end,intStimType)),strSes,strRec));
				
				%plot mean act non-pref corr
				subplot(2,2,2)
				doPEP(sEvents,vecMeanActNonPref);
				ylabel('Mean non-pref population dF/F0 correct')
				title(sprintf('Stimulus %s PEP %s%s',num2str(sParamsTD.sTypes.matTypes(end,intStimType)),strSes,strRec));
				
				%plot corr diff
				figure(handleFigPNP);
				subplot(intNumSubX,intNumSubY,intStimType);
				sEvents.handleFig=handleFigPNP;
				sEvents.vecColorFill = [0.75 0.75 1];
				sEvents.vecColorLine =  [0 0 1];
				doPEP(sEvents,vecMeanActPref-vecMeanActNonPref);
				ylabel('Mean (pref) - (non-pref) population dF/F0 correct')
				title(sprintf('Stimulus %s PEP %s%s',num2str(sParamsTD.sTypes.matTypes(end,intStimType)),strSes,strRec));
				
				%set plotting parameters for incorrect trials
				vecCorrResp = logical(ses.structStim.vecTrialResponse);
				vecSelectTrials = sParamsTD.cellSelect{intStimType} & ~vecCorrResp;
				if sum(vecSelectTrials) > 0
					figure(handleFig);
					sEvents = struct;
					sEvents.handleFig=handleFig;
					sEvents.vecOn = ses.structStim.FrameOn(vecSelectTrials);
					sEvents.vecOff = ses.structStim.FrameOff(vecSelectTrials);
					sEvents.vecWindow = round([-3*ses.samplingFreq 5*ses.samplingFreq]);
					sEvents.dblFrameRate = ses.samplingFreq;
					
					%plot mean act pref incorr
					subplot(2,2,3)
					doPEP(sEvents,vecMeanActPref);
					ylabel('Mean pref population dF/F0 incorrect')
					title(sprintf('Stimulus %s PEP %s%s',num2str(sParamsTD.sTypes.matTypes(end,intStimType)),strSes,strRec));
					
					%plot mean act non-pref incorr
					subplot(2,2,4)
					doPEP(sEvents,vecMeanActNonPref);
					ylabel('Mean non-pref population dF/F0 incorrect')
					title(sprintf('Stimulus %s PEP %s%s',num2str(sParamsTD.sTypes.matTypes(end,intStimType)),strSes,strRec));
					
					%plot incorr diff
					figure(handleFigPNP);
					subplot(intNumSubX,intNumSubY,intStimType);
					sEvents.handleFig=handleFigPNP;
					sEvents.vecColorFill = [1 0.75 0.75];
					sEvents.vecColorLine =  [1 0 0];
					doPEP(sEvents,vecMeanActPref-vecMeanActNonPref);
					ylabel('Mean (pref) - (non-pref) population dF/F0 correct')
					title(sprintf('Stimulus %s PEP %s%s',num2str(sParamsTD.sTypes.matTypes(end,intStimType)),strSes,strRec));
				end
			end
		end
	end
	
	%% plot mean correlation matrix for correct decodings vs incorrect decodings
	
	%% plot mean correlation matrix for correct decodings vs incorrect decodings [split by stimulus type]
	
	%% plot mean correlation matrix for correct lever press vs incorrect lever press
	
	%% plot mean correlation matrix for correct lever press vs incorrect lever press [split by stimulus type]
	
	
	drawnow;
end

