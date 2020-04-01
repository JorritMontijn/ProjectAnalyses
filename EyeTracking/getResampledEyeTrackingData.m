function sEyeTracking = getResampledEyeTrackingData(sEyeTracking,ses,boolEndPresent)
	
	%frames to realtime seconds
	dblCamHz = 25.05;
	dblMicHz = ses.samplingFreq;
	dblMicDur = ses.time.dur;
	intMicFrames = round(dblMicHz*dblMicDur);
	if boolEndPresent
		intReqCamFrames = length(sEyeTracking.vecPosX);
		dblCamHz = intReqCamFrames/dblMicDur;
		fprintf('Actual cam Hz %sxyt%02d=%f\n',ses.session,ses.recording,dblCamHz);
	else
		intReqCamFrames = round(dblMicDur*dblCamHz);
	end
	
	%get fields
	cellFields=fieldnames(sEyeTracking);
	for intField = 1:length(cellFields)
		strField = cellFields{intField};
		if strcmp(strField(1:3),'vec') || strcmp(strField(1:3),'ind')
			%get data
			vecCamData = double(sEyeTracking.(strField));
			intCamFrames = length(vecCamData);
			if intCamFrames < intReqCamFrames
				vecCamDummy = nan(1,intReqCamFrames);
				vecCamDummy(1:intCamFrames) = vecCamData;
				vecCamData = vecCamDummy;
			else
				vecCamData = vecCamData(1:intReqCamFrames);
			end
			
			
			%resample data
			vecResampledCamData = resample(vecCamData,intMicFrames,intReqCamFrames);
			
			%rebinarize
			if strcmp(strField(1:3),'ind'),vecResampledCamData = vecResampledCamData>0.5;end
				
			%put in struct
			sEyeTracking.(strField) = vecResampledCamData;
			
			%z-score
			%resample data
			vecCamDataZ = nanzscore(vecCamData);
			vecResampledCamDataZ = resample(vecCamDataZ,intMicFrames,intReqCamFrames);
			
			%rebinarize
			if strcmp(strField(1:3),'ind'),vecResampledCamDataZ = vecResampledCamDataZ>0.5;end
			sEyeTracking.([strField '_Z']) = vecResampledCamDataZ;
		end
	end
end