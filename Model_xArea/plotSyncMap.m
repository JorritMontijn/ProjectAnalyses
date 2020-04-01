function plotSyncMap(vecEdgesX,vecEdgesY,matValMeans,dblWindow_ms,boolRemVals)
	
	%check input
	if ~exist('boolRemVals','var') || isempty(boolRemVals)
		boolRemVals = true;
	end
	
	%plot synchronization matrix
	vecPlotX = vecEdgesX(2:end)-median(diff(vecEdgesX))/2;
	vecPlotY = vecEdgesY(2:end)-median(diff(vecEdgesY))/2;
	if boolRemVals
		indKeepX = ~all(isnan(matValMeans),1);
		indKeepY = ~all(isnan(matValMeans),2);
		matValMeans = matValMeans(indKeepY,indKeepX);
		vecPlotX = vecPlotX(indKeepX);
		vecPlotY = vecPlotY(indKeepY);
		vecEdgesX = vecEdgesX(indKeepX);
		vecEdgesY = vecEdgesY(indKeepY);
	end
	
	subplot(2,2,2)
	matC = colormap(parula(numel(vecPlotX)));
	hold all;
	for intL=1:size(matC,1)
		plot(vecPlotY,matValMeans(:,intL),'Color',matC(intL,:));
	end
	xlim([vecEdgesY(1) vecEdgesY(end)]);
	hold off;
	xlabel('Synchronization rate (Hz)');
	ylabel('V2 spiking rate (Hz)');
	title('Lines: various V1 spiking rate bins');
	fixfig;
	
	subplot(2,2,3)
	matC = colormap(parula(numel(vecPlotY)));
	hold all;
	for intL=1:size(matC,1)
		plot(vecPlotX,matValMeans(intL,:),'Color',matC(intL,:));
	end
	xlim([vecEdgesX(1) vecEdgesX(end)]);
	hold off
	ylabel('V2 spiking rate (Hz)');
	xlabel('V1 spiking rate (Hz)');
	title('Lines: various synchronization rate bins');
	fixfig;
	
	subplot(2,2,1)
	imagesc(vecPlotX,vecPlotY,matValMeans);
	axis xy;
	nancolorbar(matValMeans,[],parula(256));
	title(sprintf('Average V2 spiking rate (Hz); window: %dms',dblWindow_ms));
	xlabel('Average V1 spiking rate (Hz)');
	ylabel('Average sync. rate per neuron (Hz)');
	fixfig;
	grid off;
	
	%maximize
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
end