function [intTrialType,intTrialDur] = getPlaidTrialType(vecDurs,vecDurTypes,vecOn,vecOff,vecContrast,vecOri)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	
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
	
	%get nr of stims
	intStims = length(vecOn);
	
	%set vars
	intLow = 25;
	intHigh = 50;
	
	%get type
	intTrialType = [];
	intTrialDur = [];
	for intTrial = 1:12
		if intTrial == 1
			%1)	Low on		High on		high off	low off
			if intStims > 1 && vecContrast(1) == intLow && vecContrast(2) == intHigh && vecOn(1) < vecOn(2)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 4;
			end
		elseif intTrial == 2
			%2)	High on		high off
			if intStims == 1 && vecContrast(1) == intHigh && vecDurs(1) == vecDurTypes(1)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 2;
			end
		elseif intTrial == 3
			%3)	Low on								low off
			if intStims == 1 && vecContrast(1) == intLow && vecDurs(1) == vecDurTypes(3)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 4;
			end
		elseif intTrial == 4
			%4)	High on		low on		low off		high off
			if intStims > 1 && vecContrast(1) == intHigh && vecContrast(2) == intLow && vecOn(1) < vecOn(2)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 4;
			end
		elseif intTrial == 5
			%5)	low on		low off
			if intStims == 1 && vecContrast(1) == intLow && vecDurs(1) == vecDurTypes(1)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 2;
			end
		elseif intTrial == 6
			%6)	High on								high off
			if intStims == 1 && vecContrast(1) == intHigh && vecDurs(1) == vecDurTypes(3)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 4;
			end
		elseif intTrial == 7
			%7)	H+L on		high off	low off
			if intStims > 1 && vecOn(1) == vecOn(2) && vecContrast(1) == intHigh && vecDurs(2) == vecDurTypes(2)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 3;
			end
		elseif intTrial == 8
			%8)	H+L on		low off		high off
			if intStims > 1 && vecOn(1) == vecOn(2) && vecContrast(1) == intLow && vecDurs(2) == vecDurTypes(2)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 3;
			end
		elseif intTrial == 9
			%9)	Low on		low on		low off		low off
			if intStims > 1 && vecContrast(1) == intLow && vecContrast(2) == intLow && vecOn(1) < vecOn(2)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 4;
			end
		elseif intTrial == 10
			%10)	High on		High on		high off	high off
			if intStims > 1 && vecContrast(1) == intHigh && vecContrast(2) == intHigh && vecOn(1) < vecOn(2)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 4;
			end
		elseif intTrial == 11
			%11)	L+L on		L+L off	
			if intStims > 1 && vecContrast(1) == intLow && vecContrast(2) == intLow && vecOn(1) == vecOn(2)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 1;
			end
		elseif intTrial == 12
			%12)	H+H on		H+H off	
			if intStims > 1 && vecContrast(1) == intHigh && vecContrast(2) == intHigh && vecOn(1) == vecOn(2)
				intTrialType(end+1) = intTrial;
				intTrialDur(end+1) = 1;
			end
		end
	end

end

