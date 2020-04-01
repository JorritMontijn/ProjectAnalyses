function sOut = calcPsychometricCurve(sIn)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	%sBehavior can be either ses.structStim or sParTP structure
	
	
	
	%get data location
	if isfield(sIn,'Contrast')
		sBehavior = sIn;
	elseif isfield(sIn.structStim,'Contrast')
		sBehavior = sIn.structStim;
	end
	
	%retrieve and reformat data
	cellFields = cell(1);
	cellFields{1} = 'Contrast';
	sTypes = getStimulusTypes(sBehavior,cellFields);
	cellSelect = getSelectionVectors(sBehavior,sTypes);
	intNumContrasts = length(cellSelect);
	vecCorrect=zeros(1,intNumContrasts);
	vecRTs=nan(1,intNumContrasts);
	matRTs=nan(sum(cellSelect{1}),intNumContrasts);
	vecAllRTs = sBehavior.vecTrialRespSecs - sBehavior.SecsOn;
	
	%loop through contrasts
	for intContrast=1:intNumContrasts
		vecTempResponses = sBehavior.vecTrialResponse(cellSelect{intContrast});
		vecTempRTs = vecAllRTs(cellSelect{intContrast});
		
		%remove responses <100ms
		vecTempResponses(vecTempRTs<0.1) = 0;
		
		%get behavioral mean
		dblPercCorrect = sum(vecTempResponses==1)/numel(vecTempResponses);
		dblMeanRT = mean(vecTempRTs(vecTempResponses~=0));
		vecData = vecTempRTs(vecTempResponses~=0);
		
		%put in output
		cellResponses{intContrast} = vecTempResponses;
		vecCorrect(intContrast) = dblPercCorrect;
		vecRTs(intContrast) = dblMeanRT;
		matRTs(1:length(vecData),intContrast) = vecData;
	end
	
	%perform chi-square tests
	matP = nan(intNumContrasts,intNumContrasts);
	for intC1=1:intNumContrasts
		vecRespC1 = cellResponses{intC1};
		for intC2=intC1:intNumContrasts
			vecRespC2 = cellResponses{intC2};
			matP(intC2,intC1) = chiSquare(vecRespC1,vecRespC2);
		end
	end
	
	%put in output
	sOut.cellResponses = cellResponses;
	sOut.vecCorrect = vecCorrect;
	sOut.vecRTs = vecRTs;
	sOut.matRTs = matRTs;
	sOut.matP = matP;
end

