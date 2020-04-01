function  doPlotTuningCurve(vecIn,structStim,boolDoubleCurve)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%get responses per orientation
	sTypes = getStimulusTypes(structStim,{'Orientation'});
	cellSelect = getSelectionVectors(structStim,sTypes);
	if ~exist('boolDoubleCurve','var')
		boolDoubleCurve = true;
	end
	
	%
	intTypes = sTypes.vecNumTypes(1);
	cellResps = cell(1,intTypes);
	for intOri=1:intTypes
		cellResps{intOri} = vecIn(cellSelect{intOri});
	end
	vecMean = cellfun(@mean,cellResps);
	vecSEM = cellfun(@std,cellResps)./sqrt(cellfun(@numel,cellResps));
	vecOris = sTypes.matTypes;
	
	
	%check if we have to mirror
	if max(vecOris<pi) && boolDoubleCurve
		vecOris = [vecOris vecOris+pi];
		vecMean = [vecMean vecMean];
		vecSEM = [vecSEM vecSEM];
	end
	
	%close loop
	vecOris = [vecOris vecOris(1)];
	vecMean = [vecMean vecMean(1)];
	vecSEM = [vecSEM vecSEM(1)];
	
	
	%plot
	%{
	subplot(2,2,1)
	polar(vecOris,vecMean,'b-');
	hold on
	polar(vecOris,vecMean-vecSEM,'b--');
	polar(vecOris,vecMean+vecSEM,'b--');
	hold off
	
	subplot(2,2,2)
	%}
	plot(vecOris(1:end-1),vecMean(1:end-1),'b-');
	hold on
	plot(vecOris(1:end-1),vecMean(1:end-1)-vecSEM(1:end-1),'b--');
	plot(vecOris(1:end-1),vecMean(1:end-1)+vecSEM(1:end-1),'b--');
	hold off
	
%end

