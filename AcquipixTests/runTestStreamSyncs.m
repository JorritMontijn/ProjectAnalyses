

%% load meta files
strRoot = 'P:\Montijn\DataNeuropixels';
sMetaFilesNI = dir([strRoot filesep '**' filesep '*nidq.meta']);
sMetaFilesImAp = dir([strRoot filesep '**' filesep '*ap.meta']);
sSyncFiles = dir([strRoot filesep '**' filesep 'syncSY.mat']);
cellPathsImAp = {sMetaFilesImAp.folder};
cellPathsSync = {sSyncFiles.folder};

vecDiffsT0_Corr = nan(size(sMetaFilesNI));
vecDiffsT0 = nan(size(sMetaFilesNI));
for intRecNI=1:numel(sMetaFilesNI)
	%% get ImAp file
	fprintf('%d/%d: %s [%s]\n',intRecNI,numel(sMetaFilesNI),sMetaFilesNI(intRecNI).name,getTime);
	strTargetFolderNI = sMetaFilesNI(intRecNI).folder;
	intRecIm = find(contains(cellPathsImAp,strTargetFolderNI));
	if isempty(intRecIm)
		fprintf('%s has no ImAp Entry\n',strTargetFolderNI);
		continue;
	end
	
	%% get NI clock
	strFileMetaNI = fullpath(sMetaFilesNI(intRecNI).folder,sMetaFilesNI(intRecNI).name);
	sMetaNI = DP_ReadMeta(strFileMetaNI);
	dblSampRateReportedNI = DP_SampRate(sMetaNI);
	intFirstSampleNI = str2double(sMetaNI.firstSample);
	dblT0_NI_Reported = intFirstSampleNI/dblSampRateReportedNI;
	
	%% get ImAp clock
	strFileMetaImAp = fullpath(sMetaFilesImAp(intRecIm).folder,sMetaFilesImAp(intRecIm).name);
	sMetaImAp = DP_ReadMeta(strFileMetaImAp);
	dblSampRateReportedImAp = DP_SampRate(sMetaImAp);
	intFirstSampleImAp = str2double(sMetaImAp.firstSample);
	dblT0_ImAp_Reported = intFirstSampleImAp/dblSampRateReportedImAp;
	
	
	%% calc diff
	vecDiffsT0(intRecNI) = dblT0_ImAp_Reported - dblT0_NI_Reported;
	
	%% get sync file
	intRecSync = find(contains(cellPathsSync,strTargetFolderNI));
	if isempty(intRecSync)
		fprintf('%s has no syncSY entry\n',strTargetFolderNI);
		continue;
	end
	
	%% get corrected clocks
	%get file
	strPathSync = cellPathsSync{intRecSync};
	sSyncImAp = load(fullpath(strPathSync,'syncSY.mat'));
	boolVecSyncPulsesImAp = sSyncImAp.syncSY;
	
	%pre-allocate
	vecT0_new = nan(1,2);
	vecCorrFactor_new = nan(1,2);
	cellStreamStr = {'NI','ImAp'};
	strT0 = '    ';
	for intSource=1:2
		%% get source data
		strStream = cellStreamStr{intSource};
		if intSource == 1
			%NI
			
			%reported rates
			intFirstSample = str2double(sMetaNI.firstSample);
			dblRateFromMetaData = str2double(sMetaNI.niSampRate);
			dblRecLength = str2double(sMetaNI.fileTimeSecs);
			
			%actual pulses
			matDataNI = -DP_ReadBin(-inf, inf, sMetaNI, strrep(sMetaFilesNI(intRecNI).name,'.meta','.bin'), strTargetFolderNI); %1=PD,2=sync pulse
			intSyncPulseCh = PP_FindSyncPulseCh(matDataNI);
			[boolVecSyncPulses,dblCritValSP] = DP_GetUpDown(matDataNI(intSyncPulseCh,:));
			
		elseif intSource == 2
			%ImAp
			
			%reported rates
			intFirstSample = str2double(sMetaImAp.firstSample);
			dblRateFromMetaData = str2double(sMetaImAp.imSampRate);
			dblRecLength = str2double(sMetaImAp.fileTimeSecs);
			
			%actual pulses
			strPathSync = cellPathsSync{intRecSync};
			sSyncImAp = load(fullpath(strPathSync,'syncSY.mat'));
			boolVecSyncPulses = sSyncImAp.syncSY;
		else
			%unknown
			error
		end
		
		%% calc real rate
		%get ImAp sync pulses
		dblSampRate = PP_GetPulseIntervalFromBinVec(boolVecSyncPulses);
		
		%compare real and pre-calibrated rate
		dblRateError=(1-(dblSampRate/dblRateFromMetaData));
		dblRateErrorPercentage  = dblRateError*100;
		
		%max deviation
		dblMaxFault = dblRecLength*dblRateError;
		fprintf('%s stream; %.4f%% error gives max fault of %.0f ms; reported rate is %.6f Hz, real rate is %.6f Hz\n',...
			strStream,dblRateErrorPercentage,dblMaxFault*1000,dblRateFromMetaData,dblSampRate);
		
		%corrected
		dblT0_new = intFirstSample/dblSampRate; %the true onset
		dblCorrectionfactor = dblRateFromMetaData/dblSampRate;
		strT0 = [strT0 sprintf('T0 %s = %.3f s; ',strStream,dblT0_new)];
		
		%% save data
		vecT0_new(intSource) = dblT0_new;
		vecCorrFactor_new(intSource) = dblCorrectionfactor;
	end
		
	%% calc diff
	vecDiffsT0_Corr(intRecNI) = diff(vecT0_new);
	strT0 = [strT0 sprintf('delta-T0 = %.1f ms',1000*vecDiffsT0_Corr(intRecNI))];
	disp(strT0);
end

%% make plot
histogram(vecDiffsT0_Corr*1000);
xlim(max(abs(get(gca,'xlim')))*[-1 1]);
xlabel('NI - ImAp \DeltaT0 (ms)');
ylabel('# of recordings');