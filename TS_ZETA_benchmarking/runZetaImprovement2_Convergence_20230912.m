%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
cellUniqueAreas = {...
    'HeteroPoissonPeak',...Area 1
    '',...Area 2
    '',...Area 3
    '',...Area 4
    '',...Area 5
    '',...Area 6
    'lateral geniculate',...Area 7
    'Primary visual',...Area 8
    };
intArea = 1;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'E:\DataPreProcessed\';
end

strDataTargetPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1];%[1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecResamples = 50:50:500;%[10:10:100];
intUseGenN = 1;
intResampRepeats = 10000;
vecTrialNum = [];
%set var

for intRandType=vecRandTypes
    %reset vars
    clearvars -except strDataSourcePath intResampRepeats intUseGenN vecTrialNum vecRestrictRange cellRepStr intRandType vecRandTypes intRunStim vecRunStim cellRunStim intArea vecRunAreas cellUniqueAreas boolSave vecResamples strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRunTypes
    strArea = cellUniqueAreas{intArea};
   
    if intRandType == 1
        strRunType = strArea;
        fprintf('Prepping normal... [%s]\n',getTime);
    elseif intRandType ==2
        strRunType = [strArea '-Rand'];
        fprintf('Prepping random... [%s]\n',getTime);
    end


    %% pre-allocate output variables
    intResamps = numel(vecResamples);
    intNeurons = intUseGenN;
    cellNeuron = cell(1,intNeurons);
    vecNumSpikes = zeros(1,intNeurons);
    matZetaLinear = zeros(intResamps,intNeurons,intResampRepeats);
    matZetaUniform = zeros(intResamps,intNeurons,intResampRepeats);
    
    %% analyze
    hTic = tic;
    for intNeuron=1:intNeurons%31
        %% message
        if toc(hTic) > 5
            fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
            hTic=tic;
        end
        clear vecTrialStarts;

        %% generate
        strRecIdx = 'x';
        strMouse = 'Artificial';
        strBlock = '1';
        strDate = getDate();
        intSU = intNeuron;
        intClust = intNeuron;

        %set parameters
        dblBaseRate = exprnd(5);
        dblPrefRate = dblBaseRate+exprnd(20);
        dblKappa = rand(1)*5+5;
        vecTrialAngles=repmat([0:45:359],[1 20]);
        vecTrialDur=linspace(0.2,2,numel(vecTrialAngles));
        vecStimOnTime = cumsum(vecTrialDur)+1;
        vecStimOffTime = vecStimOnTime + 1;
        dblJitter = 0.05;
        boolDoublePeaked = true;
        intAddSpikes = 0;%numel(vecTrialDur)/10;

        vecTrialStarts(:,1) = vecStimOnTime;
        vecTrialStarts(:,2) = vecStimOffTime;
        [vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,vecTrialStarts,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,[],intAddSpikes);

        for intResampIdx = 1:numel(vecResamples)
            intResampNum = vecResamples(intResampIdx);
            %% message
            fprintf('Processing %s, N%d/%d, resampling %d (%d/%d) [%s]\n',strRunType,intNeuron,intNeurons,intResampNum,intResampIdx,numel(vecResamples),getTime);
            hTic=tic;


            %% get visual responsiveness
            %get trial dur
            dblUseMaxDur = round(median(diff(vecTrialStarts(:,1)))*2)/2;
            %set derivative params
            if contains(strRunType,'Rand')
                dblDur = dblUseMaxDur;
                vecJitter = (2*dblDur*rand([numel(vecStimOnTime) 1])-dblDur);
                matEventTimes = bsxfun(@plus,vecTrialStarts,vecJitter);
            else
                matEventTimes = vecTrialStarts;
            end
            vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
            intTrials = size(matEventTimes,1);
            intSpikeNum = numel(vecSpikeTimes);
            if intSpikeNum>50000 || intSpikeNum<3,continue;end

            intLatencyPeaks = 0;
            intPlot = 0;
            boolDirectQuantile = false;
            dblJitterSize = [];
            boolStitch = true;

            %run repeats
            for intResampRepeat=1:intResampRepeats
                intJitterDistro = 1;
                dblZetaP_Linear = zetatest(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolDirectQuantile,dblJitterSize,boolStitch,intJitterDistro);
                intJitterDistro = 2;
                dblZetaP_Uniform = zetatest(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolDirectQuantile,dblJitterSize,boolStitch,intJitterDistro);

                matZetaLinear(intResampIdx,intNeuron,intResampRepeat) = dblZetaP_Linear;
                matZetaUniform(intResampIdx,intNeuron,intResampRepeat) = dblZetaP_Uniform;
            end

            %%
            % assign data
            cellNeuron{intNeuron} = [strArea strDate 'N' num2str(intSU)];
            vecNumSpikes(intNeuron) = intSpikeNum;
        end
    end
    if boolSave
        save([strDataTargetPath 'ZetaConvergence' strRunType 'Repeats' num2str(intResampRepeat) '.mat' ],...
            'vecResamples','cellNeuron','vecNumSpikes','matZetaLinear','matZetaUniform');
    end
end