%% load data and define groups
%strDataPath
%cellUseForEyeTracking
%strTargetPath
%cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
%cellUseAreas{2} = {'nucleus of the optic tract'};
%cellUseAreas{3} = {'superior colliculus'};
%cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
%cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};
%cellSubjectGroups = {'BL6','DBA'};
clear all;
runHeaderNOT;
intRunMaxArea = 2;

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);
%vecSTA_edges = [-4:0.1:4];
vecSTA_edges = [-0.5:0.05:0.5];

%% get data
if 0%exist(fullpath(strTargetPath,'PupilPSTH.mat'),'file')
	%load
	sLoad = load(fullpath(strTargetPath,'PupilPSTH.mat'));
	cellNames = sLoad.cellNames;
	cellSTA = sLoad.cellSTA;
	cellZetaP = sLoad.cellZetaP;
	matSTAZetaP = sLoad.matSTAZetaP;
	vecSTA_edges = sLoad.vecSTA_edges;
else
	%pre-allocate
	dblMoveEventThreshold = 0.05;
	cellSTA = cellfill(nan(0,numel(vecSTA_edges)-1),[3 intRunMaxArea 0]); %[move-type x V1/NOT x rec]
	matSTAZetaP = nan([3 intRunMaxArea 0]); %[move-type x V1/NOT x rec]
	cellZetaP = cell([3 intRunMaxArea 0]); %[move-type x V1/NOT x rec]
	
	%run
	cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
	cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
	cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
	indRemRecs = contains(cellExperiment,cellRemove);
	indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
	cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
	for intSubType=1
		if intSubType == 1
			strSubjectType = 'BL6';
			dblOffsetT=0;
			dblAverageMouseHeadTiltInSetup = -15; %average tilt due to head-bar placement; procedures changed between BL6 and DBA experiments
			boolInvertX = 1; %other eye was recorded, so temporonasal is nasotemporal
		elseif intSubType == 2
			strSubjectType = 'DBA';
			dblOffsetT=0;
			dblAverageMouseHeadTiltInSetup = -15;
			boolInvertX = 0;
		end
		indUseRecs = contains(cellSubjectType,strSubjectType);
		vecRunRecs = find(indUseRecs & ~(indRemRecs | indRemRecs2));
		
		for intRecIdx=1:numel(vecRunRecs)
			intRec=vecRunRecs(intRecIdx);
			sRec = sExp(intRec);
			strRec = sRec.Name(1:(end-7));
			fprintf('Running %s (%d/%d) [%s]\n',strRec,intRecIdx,numel(vecRunRecs),getTime);
			
			%% get stim absence
			cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
			vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
			vecBlocksNM = find(contains(cellStimType,'naturalmovie','IgnoreCase',true));
			%get timing for DG
			intBlock = vecBlocksDG(1);
			sBlock = sRec.cellBlock{intBlock};
			if isfield(sBlock,'vecPupilStimOn')
				vecPupilLatency = sBlock.vecPupilStimOn-sBlock.vecStimOnTime;
			else
				vecPupilLatency = 0;
			end
			dblPL = median(vecPupilLatency);
			vecAllStimOnsets = [];
			vecAllStimOffsets = [];
			vecAllStimBlock = [];
			dblRecDur = max(cell2vec({sRec.sCluster.SpikeTimes}));
			for intBlockIdx=1:numel(sRec.cellBlock)
				sBlock = sRec.cellBlock{intBlockIdx};
				vecAllStimOnsets = cat(2,vecAllStimOnsets,sBlock.vecStimOnTime);
				vecAllStimOffsets = cat(2,vecAllStimOffsets,sBlock.vecStimOffTime);
				vecAllStimBlock = [];
			end
			vecISI = [vecAllStimOnsets dblRecDur] - [0 vecAllStimOffsets];
			[dblLongestBlankPeriod,intStim] = max(vecISI);
			dblBlankStart = vecAllStimOffsets(intStim-1);
			
			%get eye movement events
			sPupil = sRec.sPupil;
			if ~isfield(sPupil,'vecBlinks')
				sPupil.vecBlinks = false(size(sPupil.vecTime));
			end
			dblPupilStart = dblBlankStart+dblPL+10;
			dblPupilStop = dblBlankStart+dblLongestBlankPeriod+dblPL-10;
			indUsePupil = sPupil.vecTime>dblPupilStart & sPupil.vecTime<dblPupilStop & ~sPupil.vecBlinks;
			vecPupilT = sPupil.vecTime(indUsePupil);
			vecPupilX = sPupil.vecCenterX(indUsePupil);
			vecPupilY = sPupil.vecCenterY(indUsePupil);
			vecdT = diff(vecPupilT);
			vecdX = diff(vecPupilX);
			vecdY = diff(vecPupilY);
			hFig=figure;
			maxfig(hFig);
			for intMoveType = 1:3
				if intMoveType==1
					vecdM = vecdX./vecdT;
					strMoveType = 'horz';
				elseif intMoveType==2
					vecdM = vecdY./vecdT;
					strMoveType = 'vert';
				elseif intMoveType==3
					vecdM = sqrt(vecdX.^2 + vecdY.^2)./vecdT;
					strMoveType = 'tot';
				else
					error('not possible');
				end
				
				%5% cut-off
				[vecSortdM,vecReorder] = sort(vecdM);
				intCutOff = floor(numel(vecSortdM)*dblMoveEventThreshold);
				vecUpperValIdx = vecReorder((end-intCutOff-1):end);
				indUpperVals = false(size(vecdM));
				indUpperVals(vecUpperValIdx) = true;
				vecUpperVals = vecdM(vecUpperValIdx);
				
				%get epochs
				vecStartMoveIdx = 1+find(diff(indUpperVals)==1);
				vecStopMoveIdx = find(diff(indUpperVals)==-1);
				if vecStartMoveIdx(1) > vecStopMoveIdx(1)
					vecStartMoveIdx = cat(2,1,vecStartMoveIdx);
				end
				if vecStartMoveIdx(end) > vecStopMoveIdx(end)
					vecStopMoveIdx = cat(2,vecStopMoveIdx,vecStartMoveIdx(end));
				end
				vecMoveDurSecs = nan(size(vecStopMoveIdx));
				for intEpoch=1:numel(vecMoveDurSecs)
					vecMoveDurSecs(intEpoch) = sum(vecdT(vecStartMoveIdx(intEpoch):vecStopMoveIdx(intEpoch)));
				end
				dblMedianFrameDur = median(vecdT);
				vecCutOffs = dblMedianFrameDur*(0.5+(1:(max(vecStopMoveIdx-vecStartMoveIdx)+1)));
				vecStartMoveT = vecPupilT(vecStartMoveIdx)-dblPL;
				
				%remove overlapping starts
				dblMinDist = range(vecSTA_edges)/2;
				vecCheckEvents = find(diff(vecStartMoveT)<dblMinDist)+1;
				indRemEvent = false(size(vecStartMoveT));
				for intEvIdx=1:numel(vecCheckEvents)
					intEv = vecCheckEvents(intEvIdx);
					vecDiffT = vecStartMoveT(intEv) - vecStartMoveT(~indRemEvent);
					if any(vecDiffT > 0 & vecDiffT < dblMinDist)
						indRemEvent(intEv) = true;
					end
				end
				vecStartMoveT(indRemEvent) = [];
				
				%% spike-triggered average per area
				indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
				cellAreasPerCluster = {sRec.sCluster.Area};
				for intArea=1:intRunMaxArea
					% select cells
					strAreaGroup =  cellAreaGroupsAbbr{intArea};
					vecSelectCells = find(indUseCells(:) & flat(contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true)));
					if isempty(vecSelectCells)
						continue;
					end
					cellSpikes = {sRec.sCluster(vecSelectCells).SpikeTimes};
					
					%select only responsive cells
					vecZetaP = ones(1,numel(cellSpikes));
					intResampNum=100;
					for intCell=1:numel(cellSpikes)
						vecZetaP(intCell) = zetatest(cellSpikes{intCell},vecStartMoveT+vecSTA_edges(1),range(vecSTA_edges),intResampNum,0);
					end
					vecZetaP(vecZetaP>0.05)=[];
					cellSpikes(vecZetaP>0.05)=[];
					intCellNum = numel(cellSpikes);
					if intCellNum==0,continue;end
					vecSpikeT = sort(cell2vec(cellSpikes));
					%reduce spikes to max 100k
					intMaxSpikes = inf;
					if numel(vecSpikeT) > intMaxSpikes
						vecUseSpikes = unique(round(linspace(1,numel(vecSpikeT),intMaxSpikes)));
					else
						vecUseSpikes = 1:numel(vecSpikeT);
					end
					
					%% plot grand average over all cells
					figure(hFig);
					hAx = subplot(2,3,intMoveType+(intArea-1)*3);
					[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecSpikeT(vecUseSpikes),vecSTA_edges,vecStartMoveT,hAx);
					intResampNum=250;
					intPlot=0;
					[dblZetaP,sZeta] = zetatest(vecSpikeT(vecUseSpikes),vecStartMoveT+vecSTA_edges(1),range(vecSTA_edges),intResampNum,intPlot);
					%dblZetaP=nan;
					ylabel(sprintf('Total firing rate in %s (Hz)',strAreaGroup));
					xlabel(sprintf('Time after %s eye-movement (s)',strMoveType));
					title(sprintf('%s (%s) - %s %s; Z-p=%.3e [N=%d cells]',strRec,strSubjectType,strAreaGroup,strMoveType,dblZetaP,intCellNum),'interpreter','none');
					
					%% save data
					cellSTA{intMoveType,intArea,intRecIdx} = matPET; %[move-type x V1/NOT x rec]
					matSTAZetaP(intMoveType,intArea,intRecIdx) = dblZetaP; %[move-type x V1/NOT x rec]
					cellZetaP{intMoveType,intArea,intRecIdx} = vecZetaP; %[move-type x V1/NOT x rec]
				end
			end
			
			%% save figure
			drawnow;
			strFigName = ['PupilMoveSpikePSTH_' strRec];
			export_fig(fullpath(strTargetPath,[strFigName '.tif']));
			export_fig(fullpath(strTargetPath,[strFigName '.pdf']));
		end
	end
	
	%% save data
	cellNames = {sExp(vecRunRecs).Name};
	save(fullpath(strTargetPath,'PupilPSTH.mat'),'cellSTA','matSTAZetaP','cellZetaP','cellNames','vecSTA_edges');
end
%% plot 1: cell-level
%transform p2z
cellZetaZ = cellfun(@(x) -norminv(x/2),cellZetaP,'UniformOutput',false);
indRemRec = squeeze(any(any(cellfun(@isempty,cellZetaZ),2),1));
cellZetaZ_red = cellZetaZ(:,:,~indRemRec);
cellCompZ = cell(3,2);
vecBinsE=0:0.5:6;
vecBinsC = vecBinsE(2:end)-median(diff(vecBinsE));
figure;maxfig;
for intMoveType=1:3
	subplot(2,3,intMoveType);
	hold on;
	if intMoveType==1
		strMoveType = 'horz';
	elseif intMoveType==2
		strMoveType = 'vert';
	elseif intMoveType==3
		strMoveType = 'tot';
	else
		error('not possible');
	end
	for intArea=1:intRunMaxArea
		if intArea==1
			varCol = 'b';
		elseif intArea==2
			varCol = 'r';
		else
			error('not possible');
		end
		
		cellCompZ{intMoveType,intArea} = cell2vec(cellZetaZ_red(intMoveType,intArea,:));
		vecBinsV = histcounts(cellCompZ{intMoveType,intArea},vecBinsE);
		plot(vecBinsC,vecBinsV./sum(vecBinsV),varCol);
	end
	xlabel('Z-score');
	ylabel('Normalized fraction of cells');
	title(sprintf('%s; spiking modulation by eye movement',strMoveType));
	fixfig; grid off;
	
	%plot significant cells
	subplot(2,3,intMoveType+3);
	hold on;
	dbl2Sd = -norminv(0.05/2);
	dblAlphaSd = (1-normcdf(1))*2;
	k1 = sum(cellCompZ{intMoveType,1}>dbl2Sd);
	n1 = numel(cellCompZ{intMoveType,1});
	k2 = sum(cellCompZ{intMoveType,2}>dbl2Sd);
	n2 = numel(cellCompZ{intMoveType,2});
	[phatV1,pciV1] = binofit(k1,n1,dblAlphaSd);
	[phatNot,pciNot] = binofit(k2,n2,dblAlphaSd);
	[pBino,z] = bino2test(k1,n1,k2,n2);
	
	errorbar(0.2,phatV1,phatV1-pciV1(1),phatV1-pciV1(2),'xb');
	errorbar(0.8,phatNot,phatNot-pciNot(1),phatNot-pciNot(2),'xr');
	xlim([0 1])
	set(gca,'xtick',[0.2 0.8],'xticklabel',{'V1','NOT'});
	ylabel('Fraction of significant cells (mean +/- sd)');
	title(sprintf('%s; 2-sample binomial test, p=%.3f',strMoveType,pBino));
	fixfig; grid off;
	
end

%save fig
drawnow;
strFigName = ['PupilMoveCellModulation'];
export_fig(fullpath(strTargetPath,[strFigName '.tif']));
export_fig(fullpath(strTargetPath,[strFigName '.pdf']));

%% plot 2: rec-level
vecSTAC = vecSTA_edges(2:end)-mean(diff(vecSTA_edges))/2;
figure;maxfig;
for intMoveType=1:3
	if intMoveType==1
		strMoveType = 'horz';
	elseif intMoveType==2
		strMoveType = 'vert';
	elseif intMoveType==3
		strMoveType = 'tot';
	else
		error('not possible');
	end
	matMeanCtx = [];
	matMeanNot = [];
	for intRec=1:size(cellSTA,3)
		if ~isempty(cellSTA{intMoveType,1,intRec}),matMeanCtx = cat(1,matMeanCtx,zscore(mean(cellSTA{intMoveType,1,intRec})));end
		if ~isempty(cellSTA{intMoveType,2,intRec}),matMeanNot= cat(1,matMeanNot,zscore(mean(cellSTA{intMoveType,2,intRec})));end
	end
	
	
	%ctx
	subplot(2,3,intMoveType)
	hold on;
	errorfill(vecSTAC,mean(matMeanCtx,1),std(matMeanCtx,[],1)./sqrt(size(matMeanCtx,1)),[0 0 1]);
	hold off
	%p-values
	intBins = size(matMeanCtx,2);
	vecCtxP = nan(1,intBins);
	for intBin=1:intBins
		[h,p]=ttest(matMeanCtx(:,intBin));
		vecCtxP(intBin)=p;
	end
	xlabel('Time after pupil movement (s)');
	ylabel('Spiking rate (z-score)');
	
	[hCtx,c,vecCtxP_corr] = fdr_bh(vecCtxP);
	strSig = '';
	for i=find(hCtx)
		strSig = [strSig sprintf('%.3fs=%.4f; ',vecSTAC(i),vecCtxP_corr(i))];
	end
	title([sprintf('Ctx, %s, BH-FDR p<0.05: ',strMoveType) strSig]);
	fixfig;
	
	%not
	subplot(2,3,3+intMoveType)
	hold on;
	errorfill(vecSTAC,mean(matMeanNot,1),std(matMeanNot,[],1)./sqrt(size(matMeanNot,1)),[1 0 0]);
	hold off
	xlabel('Time after pupil movement (s)');
	ylabel('Spiking rate (z-score)');
	
	
	vecNotP = nan(1,intBins);
	for intBin=1:intBins
		[h,p]=ttest(matMeanNot(:,intBin));
		vecNotP(intBin)=p;
	end
	vecNotP(intBin)=p;
	
	[hNot,c,vecNotP_corr] = fdr_bh(vecNotP);
	strSig = '';
	for i=find(hNot)
		strSig = [strSig sprintf('%.3fs=%.3f; ',vecSTAC(i),vecNotP_corr(i))];
	end
	title([sprintf('NOT, %s, BH-FDR p<0.05: ',strMoveType) strSig]);
	fixfig;
end
%% plot 3
matAvg=cell2mat(cellfun(@(x) reshape(zscore(mean(x)),[1 1 1 numel(vecSTA_edges)-1]),cellSTA(:,:,~indRemRec),'UniformOutput',false));
intUseRecs = sum(~indRemRec);

%smooth
matAvgSmooth = matAvg;
dblSmoothWidth = 0.4;
vecSmooth = normpdf(-10:10,0,dblSmoothWidth) ./ sum(normpdf(-10:10,0,dblSmoothWidth));
for intMoveType=1:3
	for intArea=1:intRunMaxArea
		for intRec=1:intUseRecs
			matAvgSmooth(intMoveType,intArea,intRec,:) = imfilt(squeeze(matAvg(intMoveType,intArea,intRec,:))',vecSmooth);
		end
	end
end
matMean = squeeze(mean(matAvgSmooth,3));
matSd = squeeze(std(matAvgSmooth,[],3));
figure;maxfig;
for intMoveType=1:3
	subplot(2,3,intMoveType);
	hold on;
	if intMoveType==1
		strMoveType = 'horz';
	elseif intMoveType==2
		strMoveType = 'vert';
	elseif intMoveType==3
		strMoveType = 'tot';
	else
		error('not possible');
	end
	for intArea=1:intRunMaxArea
		if intArea==1
			varCol = [0 0 1];
		elseif intArea==2
			varCol = [1 0 0];
		else
			error('not possible');
		end
		
		%plot
		errorfill(vecSTAC,squeeze(matMean(intMoveType,intArea,:)),squeeze(matSd(intMoveType,intArea,:))./sqrt(intUseRecs),varCol);
	end
	legend({cellAreaGroupsAbbr{1},'',cellAreaGroupsAbbr{2},''})
end