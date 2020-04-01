function [meanCalciumNeuronMat, vecAngle] = Plaids_responses_subfunction(filenamePlaids, filenameTuning, percNeuronsInclude)
	

%% Calculates neurons with 50% highest OSI's
load(filenameTuning);

OSImat = [];
prefAngles = zeros(1, length(ses.neuron));

%calculates OSI's for all neurons
h = waitbar(0, 'Please wait while calculating Orientation Selectivity Indices');
for i = 1 : length(ses.neuron)
	[structActivity,matAct,vecAngle,sOut] = calcTuning(ses,i);
	close figure 1
	[pHOTori,pHOTdir,vecTuningParams,structTuningProperties] = testTuning(matAct,vecAngle,false);
	OSI = structTuningProperties.OSI;
	OSImat = [OSImat OSI];
	
	prefAngle = structTuningProperties.dblPrefAngle;
	prefAngles(i) = prefAngle;
	
	waitbar(i/length(ses.neuron));
end
close(h);

%sorts neurons based on OSI and puts neurons with 50% hightest OSI's in
%matrix
[OSIs ind] = sort(OSImat, 2, 'descend');
upperHalf = floor(length(ses.neuron)*percNeuronsInclude/100); %only 100/n % of highest OSI neurons included
vecIncludeNeurons = ind(1:upperHalf);
vecIncludeNeurons = sort(vecIncludeNeurons);

% if ~isempty(strfind(filenameTuning, '20130307xyt04')); f=find(vecIncludeNeurons==61); vecIncludeNeurons(f)=[]; end %remove outlier neuron #61 
% if ~isempty(strfind(filenameTuning, '20130307xyt01')); f=find(vecIncludeNeurons==61); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130307xyt01')); f=find(vecIncludeNeurons==106); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130307xyt01')); f=find(vecIncludeNeurons==52); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130313xyt02')); f=find(vecIncludeNeurons==102); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130313xyt02')); f=find(vecIncludeNeurons==101); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130313xyt02')); f=find(vecIncludeNeurons==38); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130313xyt02')); f=find(vecIncludeNeurons==111); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130315xyt01')); f=find(vecIncludeNeurons==6); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130315xyt01')); f=find(vecIncludeNeurons==42); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130315xyt01')); f=find(vecIncludeNeurons==44); vecIncludeNeurons(f)=[]; end %remove outlier
% if ~isempty(strfind(filenameTuning, '20130315xyt01')); f=find(vecIncludeNeurons==61); vecIncludeNeurons(f)=[]; end %remove outlier

vecPrefAngles = prefAngles(vecIncludeNeurons);

neuronsMat = [vecIncludeNeurons; vecPrefAngles]';
	
	
%%  Calculates the frame numbers
load(filenamePlaids);

	
	%{
		T1(3s)		T2(3)		T3(3s)		T4(3s)		Tot dur (+mask(1s)+ITI(5s))
	1)	Low on		High on		high off	low off			18s
	2)	High on		high off								12s
	3)	Low on								low off			18s
	4)	High on		low on		low off		high off		18s
	5)	low on		low off									12s
	6)	High on								high off		18s
	7)	H+L on		high off	low off						15s
	8)	H+L on		low off		high off					15s
	9)	Low on		low on		low off		low off			18s
	10)	High on		High on		high off	high off		18s
	%}

trialMat = stimDetection(ses);

% Create cell array for frame numbers
framesIncl = cell(10,1);
for i = 1:10
	
	framesIncl{i,1} = cell(5,1);
	for j = 1:5
		framesIncl{i,1}{j,1} = cell(2,1);
	end
	
end

% Fill in the frame numbers 
for i = 1 : length(trialMat)
	
	if trialMat(3,i) == 0
		k = 1;
	elseif trialMat(3,i) == 90
		k = 2;
	end
		
	for j = 1 : 5  
		frameStart = trialMat(5,i) + ses.samplingFreq*3*(j-1);
		framesIncl{trialMat(2,i),1}{j,1}{k,1} = [framesIncl{trialMat(2,i),1}{j,1}{k,1}; frameStart];
	end
	
end


%% Calculate dF/F values for plaids

% Situations to include
% [situation (see table above) / interval / primary grating orientation (1=0 % deg, 2= 90 deg)]
situationsMat = [10 2 1 ; ...
	10 2 2
	7 1 1
	8 2 1
	7 1 2
	9 2 1
	9 2 2
	7 2 1
	8 2 2
	7 2 2 ];

meanCalciumNeuronMat = nan(length(neuronsMat), length(situationsMat)+3);
meanCalciumNeuronMat(:,1:2) = neuronsMat;

% Calculate mean dFoF for all neurons for every situation
for i = 1 : length(neuronsMat)
	neu = neuronsMat(i);
	
	for j = 1 : length(situationsMat)
		
		frameCell = framesIncl{situationsMat(j,1), 1}{situationsMat(j,2), 1}{situationsMat(j,3), 1};
		
		dFoFVec = [];
		for k = 1 : length(frameCell)
			dFoF = mean(ses.neuron(1, neu).dFoF(round(frameCell(k)) : round(frameCell(k)+3*ses.samplingFreq)));
			dFoFVec = [dFoFVec dFoF]; %all situations are used several times, so this vector incluces dF/F values for all trials
		end
		meanCalciumNeuronMat(i, j+2) = mean(dFoFVec); %mean all the dF/F values for one situation for one neuron 
	end
end


%% Include 100% contrast gratings in columns 14 and 15 of 'meanCalciumNeuronMat'

% for 0 degree gratings
for j = 1 : length(meanCalciumNeuronMat)
	zeroVec = [];
	for i = 1 : length(ses.structStim.FrameOn)
				
		if ses.structStim.Orientation(i) == 0
			dFoFVec = ses.neuron(1, meanCalciumNeuronMat(j,1)).dFoF(ses.structStim.FrameOn(i):ses.structStim.FrameOff(i));
			zeroVec = [zeroVec dFoFVec']; %all dF/F values in 1 vector
			
		end
	end
 	meanCalciumNeuronMat(j,14) = mean(zeroVec); %calculates the mean of all dF/F values for the gratings
end

% for 90 degree gratings
for j = 1 : length(meanCalciumNeuronMat)
	zeroVec = [];
	for i = 1 : length(ses.structStim.FrameOn)
				
		if ses.structStim.Orientation(i) == 90
			dFoFVec = ses.neuron(1, meanCalciumNeuronMat(j,1)).dFoF(ses.structStim.FrameOn(i):ses.structStim.FrameOff(i));
			zeroVec = [zeroVec dFoFVec'];
			
		end
	end
 	meanCalciumNeuronMat(j,15) = mean(zeroVec);
end

end
