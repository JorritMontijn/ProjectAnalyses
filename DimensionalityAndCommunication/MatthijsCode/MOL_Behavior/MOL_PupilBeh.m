%% This script visualize the pupil change time-locked to stimulus appearance with respect to type of mouse respnse

[Data] = MOL_SelectData();

sParams.secBeforeStim   = 2;        % (in sec) determine how much time before stimuli you want to visualize
sParams.secAfterStim    = 3;    	% (in sec) determine how much time after stimuli you want to visualize

sParams.responseType    = 'hit';    % (hit, miss, correct_rej, false_alarm) choose which condition will be visualized
sParams.sort            = 'after';  % ('before, 'after') sort the heatmap using mean of activity before or after stimuli

%% GET ALL EPOCHS FROM ALL SESSION
epochMatrix = MOL_PupilBeh_epoch(Data, sParams);

%% ORGANIZE EPOCHS 
secBef      = sParams.secBeforeStim;
secAft      = sParams.secAfterStim;

probeTrial  = strcmp(Data.trialData.trialType, 'P');
correct     = Data.trialData.correctResponse;

if strcmp(sParams.responseType, 'hit')
    response =  find(correct&~probeTrial);
elseif strcmp(sParams.responseType, 'miss')
    response        = find(~correct&~probeTrial);
elseif strcmp(sParams.responseType, 'correct_rej')
    response = find(correct&probeTrial);
elseif strcmp(sParams.responseType, 'false_alarm')
    response = find(~correct&probeTrial);
end

pupilFs     = Data.pupilData.fs(1);
pupilSectionL = pupilFs*(secBef+secAft); 

Time        = linspace(-secBef,secAft,pupilSectionL);

%% TAKE EPOCHS WITH RESPECT TO PARTICULAR RESPONSE
pupMatResp   = epochMatrix(response,:);

%% SORT THE MATRIX 

addCol = pupilSectionL+1;
matCondL = length(pupMatResp);

if strcmp(sParams.sort, 'before')
    sortVec = 1:secBef*pupilFs;
elseif strcmp(sParams.sort, 'after')
    sortVec = secBef+pupilFs:pupilSectionL;
end

pupMatResp(:,addCol) = mean(pupMatResp(:,sortVec),2);
pupMatConSort = sortrows(pupMatResp,addCol);
pupMatConSort = pupMatConSort(:,1:pupilSectionL);

%% HEATMAP
figure; set(gcf,'units','normalized','Position',[0.1 0.1 0.5 0.7],'color','w')
cscale = [-2 2.5];
imagesc(pupMatConSort, cscale); hold on;
%suptitle(sprintf())
temp = ones(1,matCondL);
plot(temp*(pupilFs*secBef),1:1:matCondL ,'k','LineWidth',3)
xticks(-secBef*pupilFs:1:secAft*pupilFs-1)
set(gca, 'XTick', 0:pupilFs:pupilSectionL,'XTickLabels', num2cell(-secBef:1:secAft),'FontSize', 10)
set(gca, 'YTick', 0:50:length(pupMatConSort), 'FontSize', 10)
colorbar

%% PLOT MEAN PUPIL SIZE
pupMeanCon  = mean(pupMatResp);

figure;
plot(Time,pupMeanCon(:,1:pupilSectionL))
