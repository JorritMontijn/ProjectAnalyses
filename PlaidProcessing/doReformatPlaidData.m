function sTrialType = doReformatPlaidData(ses)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%stimuli:
	%{
		Dur1 (3s)	Dur2 (3s) 	Dur3 (3s)					[+3 ITI]
	1)	Low on		High on		high off	low off			12s
	2)	High on		high off								6s
	3)	Low on								low off			12s
	4)	High on		low on		low off		high off		12s
	5)	low on		low off									6s
	6)	High on								high off		12s
	7)	H+L on		high off	low off						9s
	8)	H+L on		low off		high off					9s
	9)	Low on		low on		low off		low off			12s
	10)	High on		High on		high off	high off		12s
	11)	L+L on		L+L off									6s
	12)	H+H on		H+H off									6s
	%}
	
	%{
%step 1
load data

%step 2

%reformat data; output:
sTrialType(intTrialType).sPrimOri(intOriType).sEpoch(intEpoch).matFrames(intRepetition,intFrame)

intTrialType = [1:10]; trial types as above
intOriType = [1 2]; 1=0; 2=90
intEpoch = [1:3]; 1-3: short, mid, long
intFrame = frame number (1 - max in epoch)
intRepetition = repetition number (1-4)

matFrames: matrix containing corresponding original frame numbers
	%}
	
	%get frames
	%%
	%get durations
	vecDursAll = ses.structStim.FrameOff - ses.structStim.FrameOn;
	vecDursAll(vecDursAll==78) = 77;
	
	%get different possible durations; lowest first; 1=mask; 2=short (1 interval); 3=mid (2 ints); 4=long (3 ints)
	vecDurTypes = findvals(vecDursAll);
	if length(vecDurTypes) ~= 3
		error([mfilename ':DurationError'],'Number of possible stimulus durations is not equal to 4: please check data and/or scripts')
	end
	intEpochDur = vecDurTypes(1);
	intMaskDur = 0;
	
	%get trials
	intTrials = max(ses.structStim.TrialNumber);
	
	%get repetitions
	intReps = (intTrials/12)/2;
	
	%pre-allocate
	for intTrialType=1:12
		for intOriType = 1:2
			for intEpoch = 1:5
				sTrialType(intTrialType).sPrimOri(intOriType).sEpoch(intEpoch).matFrames = [];
			end
		end
	end
	
	%run loop
	for intTrial=1:intTrials
		%get stimuli belonging to this trial
		indStims = ses.structStim.TrialNumber == intTrial;
		
		%number of stims in this trial
		intStims = sum(indStims);
		
		%get data from this trial
		vecOn = ses.structStim.FrameOn(indStims);
		vecOff = ses.structStim.FrameOff(indStims);
		vecContrast = ses.structStim.Contrast(indStims);
		vecOri = ses.structStim.Orientation(indStims);
		vecDurs = vecDursAll(indStims);
		intTrialStart = min(vecOn);
		
		%determine trial type
		[intTrialType,intTrialDur] = getPlaidTrialType(vecDurs,vecDurTypes,vecOn,vecOff,vecContrast,vecOri);
		
		%get primary orientation
		if vecOri(1) == 0
			intOriType = 1; %0 degrees
		else
			intOriType = 2; %90 degrees
		end
		
		%go through epochs
		for intEpoch = 1:3
			if intEpoch <= intTrialDur
				
				%get frame numbers
				intStart = intTrialStart+(intEpochDur*(intEpoch-1));
				intStop = intStart+intEpochDur-1;
				
				%put in matrix
				sTrialType(intTrialType).sPrimOri(intOriType).sEpoch(intEpoch).matFrames(end+1,:) = intStart:intStop;
			end
		end

		%add mask epoch
		intStart = intStop + 1;
		intStop = intStart + intMaskDur - 1;
		sTrialType(intTrialType).sPrimOri(intOriType).sEpoch(5).matFrames(end+1,:) = intStart:intStop;
	end
end

