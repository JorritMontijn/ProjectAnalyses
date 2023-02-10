%% options
intSizeN = 64; %size of grid: 32 for each neuron, 64 for standardized axes
boolNormRunDens = true; %normalize for running density? (true)
boolRemRun0 = false; %instead of normalizing for running density, you can also remove speed=0
boolStandardizeAxes = true; %if true, you can sum all KDE maps over neurons, otherwise each neuron gets its own axis
boolPlotPerNeuron = true; %make plot for every neuron?

%standardize axes?
matAllData = signals;

%remove zero run speed?
if boolRemRun0
	indRem = runspeed==0;
else
	indRem = false(size(runspeed));
end

%log-transform
matAllData = log(matAllData(~indRem,:)+1);
vecRunSpeed = log(runspeed(~indRem)+1);

%standardize axes?
if boolStandardizeAxes
%define range
	vecMinXY = [min(vecRunSpeed(:)) min(matAllData(:))];
	vecMaxXY = [max(vecRunSpeed(:)) max(matAllData(:))];
end

%% run
for intNeuron=1:size(matAllData,2)
	vecNeuronAct = matAllData(:,intNeuron);
	
	%make matrix
	matData = cat(2,vecRunSpeed,vecNeuronAct);
	
	%define range
	if ~boolStandardizeAxes
		vecMinXY = min(matData,[],1);
		vecMaxXY = max(matData,[],1);
	end
	
	%get 2D diffusion KDE
	[bandwidth,matDens,X,Y]=kde2d(matData,intSizeN,vecMinXY,vecMaxXY);
	
	%normalize density for # of running speed samples
	if boolNormRunDens
		matDens = matDens ./ sum(matDens,1);
	end
	vecX = X(1,1:end);
	vecY = Y(1:end,1);
	vecLimX = [min(vecX) max(vecX)];
	vecLimY = [min(vecY) max(vecY)];
	
	%% plot
	if boolPlotPerNeuron
		subplot(2,3,1)
		histx(matData(:,1));
		xlim(vecLimX);
		xlabel('Running speed');
		vecTickX = get(gca,'xtick');
		vecRealX = exp(vecTickX)-1;
		set(gca,'xticklabel',roundi(vecRealX,1));
		fixfig;
		
		subplot(2,3,2);
		histx(matData(:,2));
		xlim(vecLimY);
		xlabel('dF/F0 (%)');
		vecTickY = get(gca,'xtick');
		vecRealY = exp(vecTickY)-1;
		set(gca,'xticklabel',round(vecRealY*100));
		fixfig;
		
		
		%define axes in real units
		vecTickX = round(linspace(1,intSizeN,4));
		vecTickX_val = roundi(exp(vecX(vecTickX))-1,1);
		
		vecTickY = round(linspace(1,intSizeN,4));
		vecTickY_val = round(100*(exp(vecY(vecTickY))-1));
		
		subplot(2,3,3);
		imagesc(log(matDens+1));
		set(gca,'xtick',vecTickX,'xticklabel',vecTickX_val,'ytick',vecTickY,'yticklabel',vecTickY_val);
		ylabel('dF/F0 (%)');
		xlabel('Running speed');
		axis xy
		title(sprintf('Neuron %d',intNeuron));
		fixfig;
		pause
	end
end