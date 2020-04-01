%mean vs variance; blue=baseline, red=stimulus
close all
clear all
%get block data
vecMice = [1 2 3 5 6 7];
intAnimalNum = length(vecMice);
intMouseCounter = 0;

%pre-allocate required parameters
intContrastNum = 6;
sParameters.matMeanPopResp = nan(2,intContrastNum,intAnimalNum);
sParameters.matStDevPopResp = nan(2,intContrastNum,intAnimalNum);
sParameters.matMeanPopBandwidth = nan(2,intContrastNum,intAnimalNum);
sParameters.matStDevPopBandwidth = nan(2,intContrastNum,intAnimalNum);

sParameters.matNeuronMeanStDevResp = nan(2,intContrastNum,intAnimalNum);
sParameters.matNeuronStDevStDevResp = nan(2,intContrastNum,intAnimalNum);
sParameters.matNeuronMeanFano = nan(2,intContrastNum,intAnimalNum);
sParameters.matNeuronStDevFano = nan(2,intContrastNum,intAnimalNum);
	
%loop
for intMouse = vecMice
	%get data
	if intMouse == 1
		strSes = '20140207';
	elseif intMouse == 2
		strSes = '20140314';
	elseif intMouse == 3
		strSes = '20140425';
	elseif intMouse == 4
		%return; %exclude?; bad behavior, weird signals
		strSes = '20140430';
	elseif intMouse == 5
		strSes = '20140507';
	elseif intMouse == 6
		strSes = '20140530';
	elseif intMouse == 7
		strSes = '20140604';
	end
	load(['D:\Data\Results\stimdetection\dataRawPre_aggregate' strSes '.mat']);

	%increment counter
	intMouseCounter = intMouseCounter + 1;
	fprintf('Processing mouse %d/%d [%d]: %s\n',intMouseCounter,length(vecMice),intMouse,strSes);
	
	%transform session file to ori mod 180
	ses = cellMultiSes{1};
	vecOrientationPerTrial = mod(ses.structStim.Orientation,180);
	
	%get neuronal responses per trial
	[matTrialResponse,cellSelectContrasts] = getTrialResponseData(ses,ses.structStim);
	
	%split for different contrasts
	for intC=1:length(cellSelectContrasts)
		indContrastTrials = cellSelectContrasts{intC};
		
		%split for different response trials
		for intRespType=1:2
			if intRespType == 1
				indSelectResp = ses.structStim.vecTrialResponse;
			elseif intRespType == 2
				indSelectResp = ses.structStim.vecTrialResponse~=1;
			end
			
			%get trials
			indTrials = indContrastTrials & indSelectResp;
			matResp = matTrialResponse(:,indTrials);
			vecOris = vecOrientationPerTrial(indTrials);
			
			%calculate selection vectors
			vecOriList = unique(vecOrientationPerTrial);
			intOris = length(vecOriList);
			cellSelect = cell(1,intOris);
			for intOri=1:intOris
				cellSelect{intOri} = vecOris == vecOriList(intOri);
			end
			
			%get seed parameters
			vecParameters = getSeedParameters(matResp,cellSelect,vecOriList);
			
			%assign to output
			sParameters.matMeanPopResp(intRespType,intC,intMouseCounter) = vecParameters(1);
			sParameters.matStDevPopResp(intRespType,intC,intMouseCounter) = vecParameters(2);
			sParameters.matMeanPopBandwidth(intRespType,intC,intMouseCounter) = vecParameters(3);
			sParameters.matStDevPopBandwidth(intRespType,intC,intMouseCounter) = vecParameters(4);
			
			sParameters.matNeuronMeanStDevResp(intRespType,intC,intMouseCounter) = vecParameters(5);
			sParameters.matNeuronStDevStDevResp(intRespType,intC,intMouseCounter) = vecParameters(6);
			sParameters.matNeuronMeanFano(intRespType,intC,intMouseCounter) = vecParameters(7);
			sParameters.matNeuronStDevFano(intRespType,intC,intMouseCounter) = vecParameters(8);
		end
	end
end