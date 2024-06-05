function dblError = fGenQuantDistros(dblMuFactor,dblSdFactor)
	%fGenQuantDistros Summary of this function goes here
	%   dblError = fGenQuantDistros(dblMuFactor,dblSdFactor)
	
	global cellAct;
	
	if numel(dblMuFactor)==2
		dblSdFactor=dblMuFactor(2);
		dblMuFactor=dblMuFactor(1);
	end
	
	vecMu1 = cellfun(@mean,cellAct(1,:));
	vecSd1 = cellfun(@std,cellAct(1,:));
	vecMu2 = cellfun(@mean,cellAct(2,:));
	vecSd2 = cellfun(@std,cellAct(2,:));
	vecMu = (vecMu2 - vecMu1)/2;
	vecSd = (vecSd1 + vecSd2)/2;
	vecCV = vecSd./vecMu;
	dblAvgCV = mean(vecCV);
	
	vecDeltaMu = vecMu-mean(vecMu);
	vecGenMu = dblMuFactor*vecDeltaMu+mean(vecMu);
	vecGenSd = dblSdFactor*dblAvgCV*(dblMuFactor*vecDeltaMu+mean(vecMu));
	
	vecError = nan(1,size(cellAct,2));
	for i=1:size(cellAct,2)
		vecAct = sort(cat(1,-cellAct{1,i},cellAct{2,i}))';
		
		vecP = linspace(1/numel(vecAct),1-1/numel(vecAct),numel(vecAct));
		x = norminv(vecP,vecGenMu(i),vecGenSd(i));
		
		vecError(i) = sum((vecAct - x).^2);
	end
	dblError = sum(vecError);
end

