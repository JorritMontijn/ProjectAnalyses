%function doPlotStimDetectDecoding
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

matMeanCorrectDecoded = nan(4,6); %rand/miss rand/hit HCAR/miss HCAR/hit
matSdCorrectDecoded = nan(4,6);

sDecoding = cellSaveDecoding;

%calc mean
matMeanCorrectDecoded(4,:) = mean(sDecoding(1).matMixBootstrappedDecodingOutput(:,:,1),1); %HCAR hit
matMeanCorrectDecoded(3,:) = mean(sDecoding(1).matMixBootstrappedDecodingOutput(:,:,2),1); %HCAR miss
matMeanCorrectDecoded(2,:) = mean(sDecoding(2).matMixBootstrappedDecodingOutput(:,:,1),1); %rand hit
matMeanCorrectDecoded(1,:) = mean(sDecoding(2).matMixBootstrappedDecodingOutput(:,:,2),1); %rand miss

%calc sd
matSdCorrectDecoded(4,:) = std(sDecoding(1).matMixBootstrappedDecodingOutput(:,:,1),[],1); %HCAR hit
matSdCorrectDecoded(3,:) = std(sDecoding(1).matMixBootstrappedDecodingOutput(:,:,2),[],1); %HCAR miss
matSdCorrectDecoded(2,:) = std(sDecoding(2).matMixBootstrappedDecodingOutput(:,:,1),[],1); %rand hit
matSdCorrectDecoded(1,:) = std(sDecoding(2).matMixBootstrappedDecodingOutput(:,:,2),[],1); %rand miss

%set vars
vecLineX = 1:6;
vecWindowInv = 6:-1:1;
vecX = [vecLineX vecWindowInv];

%plot
h=figure;
for intPlotType=1:4
	
	vecMeanTrace = matMeanCorrectDecoded(intPlotType,:);
	vecSD = matSdCorrectDecoded(intPlotType,:);
	
	if intPlotType == 1
		vecColorFill = [1.0 0.7 0.7];
		vecColorLine = [1 0 0];
	elseif intPlotType == 2
		vecColorFill = [0.7 1.0 0.7];
		vecColorLine = [0 1 0];
	elseif intPlotType == 3
		vecColorFill = [1.0 0.7 1.0];
		vecColorLine = [1 0 1];
	elseif intPlotType == 4
		vecColorFill = [0.7 1.0 1.0];
		vecColorLine = [0 1 1];
	end
	
	vecMinTrace = vecMeanTrace-vecSD;
	vecMaxTrace = vecMeanTrace+vecSD;
	vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
	
	%plot
	hold on
	%fill(vecX,vecY,vecColorFill,'EdgeColor','none');
	plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
	hold off
	%end
end
%legend(gca,'','Miss Random','','Hit Random','','Miss HCAR','','Hit HCAR')
legend(gca,'Miss Random','Hit Random','Miss HCAR','Hit HCAR')