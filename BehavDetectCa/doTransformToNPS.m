for intMouse=1:8
	close all
	clearvars -except intMouse
	boolUseNeuropilSubtraction = false;
	
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
			vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
		elseif intMouse == 2
			strSes = '20140314';
			vecBlock = [1 1 1 1 1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 3
			strSes = '20140425';
			vecBlock = [1 1 1 1 1 1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 4
			strSes = '20140507';
			vecBlock = [1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 5
			strSes = '20140530';
			vecBlock = [1 1 1 1 1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 6
			strSes = '20140604';
			vecBlock = [1 1 1 1]; %define whether neurons are in same population or not
		elseif intMouse == 7
			strSes = '20140711';
			vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
		elseif intMouse == 8
			strSes = '20140715';
			vecBlock = [1 1 1 1]; %define whether neurons are in same population or not
		end
		
		load(['D:\Data\Results\stimdetection\dataPreProAggregate' strSes '.mat']);
	end
	
	%recalc with nps
	fprintf('Calculating neuropil-subtracted dF/F0 for session %s [%s]\n',strSes,getTime);
	for intPopulation=1:numel(cellMultiSes)
		cellMultiSes{intPopulation} = doRecalcdFoF(cellMultiSes{intPopulation},5);
	end
	
	%save
	strSes = ['NPS' strSes];
	save(['D:\Data\Results\stimdetection\dataPreProAggregate' strSes '.mat'],'cellMultiSes','cellAggregate','vecBlock','-v7.3');
end