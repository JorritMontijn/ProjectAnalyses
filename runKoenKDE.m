%% options
intSizeN = 64; %size of grid: 32 for each neuron, 64 for standardized axes
boolNormRunDens = true; %normalize for running density? (true)
boolRemRun0 = false; %instead of normalizing for running density, you can also remove speed=0
boolPlotPerNeuron = false; %make plot for every neuron?
boolRectifyAxes = true;
%log-transform
matAllData = log(signals+1);
matRunSpeed = log(runspeed+1);

%define range
if boolRectifyAxes
	vecMinXY = [0 -0.2];
	vecMaxXY = [3.3 1.1];
else
	vecMinXY = [min(matRunSpeed(:)) min(matAllData(:))];
	vecMaxXY = [max(matRunSpeed(:)) max(matAllData(:))];
end

matAllDens = zeros(intSizeN, intSizeN, size(matAllData,2), 2);

for i = 1:2
	%% prep data
	matSubData = matAllData(:,:,i);
	vecRunSpeed = matRunSpeed(:,i);
	
	%remove zero run speed?
	if boolRemRun0
		indRem = vecRunSpeed==0;
	else
		indRem = false(size(vecRunSpeed));
	end
	matSubData = matSubData(~indRem,:);
	vecRunSpeed = vecRunSpeed(~indRem);
	
	
	%% run
	for intNeuron=1:size(matSubData,2)
		%% single neuron
		vecNeuronAct = matSubData(:,intNeuron);
		
		%make matrix
		matData = cat(2,vecRunSpeed,vecNeuronAct);
		
		%get 2D diffusion KDE
		[bandwidth,matDens,X,Y]=kde2d(matData,intSizeN,vecMinXY,vecMaxXY);
		
		%normalize density for # of running speed samples
		if boolNormRunDens
			matDens = matDens ./ (sum(matDens,1));
		end
		vecX = X(1,1:end);
		vecY = Y(1:end,1);
		vecLimX = [vecMinXY(1) vecMaxXY(1)];
		vecLimY = [vecMinXY(2) vecMaxXY(2)];
		
		%save to output
		matAllDens(:,:,intNeuron,i) = matDens;
		
		% plot
		if boolPlotPerNeuron
			%define axes in real units
			vecTickX = round(linspace(1,intSizeN,4));
			vecTickX_val = roundi(exp(vecX(vecTickX))-1,1);
			
			vecTickY = round(linspace(1,intSizeN,4));
			vecTickY_val = round(100*(exp(vecY(vecTickY))-1));
			
			
			subplot(2,3,1);maxfig;
			histx(matData(:,1));
			xlim(vecLimX);
			xlabel('Running speed');
			drawnow;
			set(gca,'xtick',roundi(log(vecTickX_val+1),1),'xticklabel',vecTickX_val);
			fixfig;
			
			subplot(2,3,2);
			histx(matData(:,2));
			xlim(vecLimY);
			xlabel('dF/F0 (%)');
			set(gca,'xtick',roundi(log(vecTickY_val/100+1),1),'xticklabel',vecTickY_val);
			fixfig;
			
			
			
			subplot(2,3,3);
			imagesc(log(matDens+1));
			colorbar;
			set(gca,'xtick',vecTickX,'xticklabel',vecTickX_val,'ytick',vecTickY,'yticklabel',vecTickY_val);
			ylabel('dF/F0 (%)');
			xlabel('Running speed');
			axis xy
			title(sprintf('Neuron %d',intNeuron));
			fixfig;
			pause
		end
	end
end

%% plot overall fig
figure;maxfig;
vecTickX = round(linspace(1,intSizeN,4));
vecTickX_val = roundi(exp(vecX(vecTickX))-1,1);

vecTickY = round(linspace(1,intSizeN,4));
vecTickY_val = round(100*(exp(vecY(vecTickY))-1));

matDensDiff = diff(matAllDens,[],4);
matDensDiffPyr = mean(matDensDiff(:,:,~reds),3);
matDensDiffChand = mean(matDensDiff(:,:,reds),3);
subplot(2,3,1)
imagesc(log(matDensDiffPyr+1));
colorbar;
set(gca,'xtick',vecTickX,'xticklabel',vecTickX_val,'ytick',vecTickY,'yticklabel',vecTickY_val);
ylabel('dF/F0 (%)');
xlabel('Running speed');
axis xy
title(sprintf('Pyr, post-pre'));
fixfig;

subplot(2,3,2)
imagesc(log(matDensDiffChand+1));
colorbar;
set(gca,'xtick',vecTickX,'xticklabel',vecTickX_val,'ytick',vecTickY,'yticklabel',vecTickY_val);
ylabel('dF/F0 (%)');
xlabel('Running speed');
axis xy
title(sprintf('Chand, post-pre'));
fixfig;
