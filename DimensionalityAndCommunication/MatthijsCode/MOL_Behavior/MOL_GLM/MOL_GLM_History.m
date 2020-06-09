function MOL_GLM_History(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};


%%

% ToDO
%- Conflict trials
%- Normalized intensity levels
%- Ridge regression L2 (how to choose lambda)

%% General settings:
% showIndFig          = 0;
% showResFig          = 1;
% sessioncounter      = 0;

%% Loop over mice and then sessions:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    
    %Select all the data from this mouse:
    mouseid         = mouseids{mou};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    
    [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each mouse individually:
    
    %% Initialize
    
    %Modality (trial type):
    X_Modality  = NaN(size(temptrialData.session_ID));
    X_Modality(strcmp(temptrialData.trialType,'Y')) = 1;
    X_Modality(strcmp(temptrialData.trialType,'X')) = 2;
    X_Modality(strcmp(temptrialData.trialType,'P')) = 3;
    X_Modality(strcmp(temptrialData.trialType,'C')) = 4;
    if any(isnan(X_Modality))
        error('unknown trial types')
    end
    
    %Visual Intensity:
    X_visualOriChange       = zeros(size(temptrialData.session_ID));
    X_visualOriChange(strcmp(temptrialData.trialType,'X')) = abs(temptrialData.visualOriChange(strcmp(temptrialData.trialType,'X')));
    X_visualOriChange(strcmp(temptrialData.trialType,'C')) = abs(temptrialData.visualOriChange(strcmp(temptrialData.trialType,'C')));
    %Log Visual Intensity:
    X_logvisualOriChange    = X_visualOriChange;
    X_logvisualOriChange(X_logvisualOriChange==0) = 0.1;
    X_logvisualOriChange    = log(X_logvisualOriChange);
    
    %Auditory Intensity:
    X_audioFreqChange       = zeros(size(temptrialData.session_ID));
    X_audioFreqChange(strcmp(temptrialData.trialType,'Y')) = abs(temptrialData.audioFreqChange(strcmp(temptrialData.trialType,'Y')));
    X_audioFreqChange(strcmp(temptrialData.trialType,'C')) = abs(temptrialData.audioFreqChange(strcmp(temptrialData.trialType,'C')));
    %Log Auditory Intensity:
    X_logaudioFreqChange    = X_audioFreqChange;
    X_logaudioFreqChange(X_logaudioFreqChange==0) = 0.1;
    X_logaudioFreqChange    = log(X_logaudioFreqChange);
    
    %% Rewards: (binary)
    X_Correct       = temptrialData.correctResponse;
    if tempsessionData.VisualLeftCorrectSide
        X_visualCorrect = temptrialData.leftCorrect & strcmp(temptrialData.responseSide,'L');
        X_audioCorrect  = temptrialData.rightCorrect & strcmp(temptrialData.responseSide,'R');
    else
        X_visualCorrect = temptrialData.rightCorrect & strcmp(temptrialData.responseSide,'R');
        X_audioCorrect  = temptrialData.leftCorrect & strcmp(temptrialData.responseSide,'L');
    end
    
    %% Photostim:
%     X_photoStim         = temptrialData.hasphotostim==1;
    
    %% Response vector:
    %Auditory response  = 1
    %Visual response    = 2
    %No response        = 3
    
    Y_Response = NaN(size(temptrialData.session_ID));
    if tempsessionData.VisualLeftCorrectSide
        Y_Response(strcmp(temptrialData.responseSide,'L')) = 2;
        Y_Response(strcmp(temptrialData.responseSide,'R')) = 1;
    else
        Y_Response(strcmp(temptrialData.responseSide,'L')) = 1;
        Y_Response(strcmp(temptrialData.responseSide,'R')) = 2;
    end
    Y_Response(temptrialData.noResponse==1) = 3;
    
    %1) Null model: random variable
    iModel              = 1;
    X_string{iModel}    = {'Null'};
    X                   = [rand(size(X_visualOriChange))];
    [B{iModel,mou},CV_Perf{iModel,mou}] = MOL_fitGLM_session(X,Y_Response);
    
    %2) Purely sensory model:
    iModel              = 2;
    X_string{iModel}    = {'VisualChange' 'AudioChange'};
    X                   = [X_visualOriChange X_audioFreqChange];
    [B{iModel,mou},CV_Perf{iModel,mou}] = MOL_fitGLM_session(X,Y_Response);
    
    %3) Log Sensory model:
    iModel              = 3;
    X_string{iModel}    = {'LogVisualChange' 'LogAudioChange'};
    X                   = [X_logvisualOriChange X_logaudioFreqChange];
    [B{iModel,mou},CV_Perf{iModel,mou}] = MOL_fitGLM_session(X,Y_Response);
    
    %4) Sensory + reward-history model:
    iModel              = 4;
    X_string{iModel}    = {'LogVisualChange' 'LogAudioChange' 'VisualCorrect-1' 'AudioCorrect-1'};
    X                   = [X_logvisualOriChange X_logaudioFreqChange [X_visualCorrect(2:end); 0] [X_audioCorrect(2:end); 0]];
    [B{iModel,mou},CV_Perf{iModel,mou}] = MOL_fitGLM_session(X,Y_Response);
    
    %5) Sensory +  sensory history model:
    iModel              = 5;
    X_string{iModel}    = {'LogVisualChange' 'LogAudioChange' 'VisualChange-1' 'AudioChange-1'};
    X                   = [X_visualOriChange X_audioFreqChange [X_visualOriChange(2:end); 0] [X_audioFreqChange(2:end); 0]];
    [B{iModel,mou},CV_Perf{iModel,mou}] = MOL_fitGLM_session(X,Y_Response);
    
    %6) Log Sensory + Opto model:
%     iModel              = 6;
%     X_string{iModel}    = {'LogVisualChange' 'LogAudioChange' 'Opto'};
%     X                   = [X_logvisualOriChange X_logaudioFreqChange X_photoStim];
%     [B{iModel,mou},CV_Perf{iModel,mou}] = MOL_fitGLM_session(X,Y_Response);
    
    %7) Reward History model:
    %         iModel              = 7;
    %         X_string{iModel}    = {'VisualCorrect-1' 'AudioCorrect-1'};
    %         X                   = [[X_visualCorrect(2:end); 0] [X_audioCorrect(2:end); 0]];
    %         [B{iModel,mou},CV_Perf{iModel,mou}] = MOL_fitGLM_session(X,Y_Response);
%     Export data to .csv file
    X                   = [X_visualOriChange X_audioFreqChange X_visualCorrect X_audioCorrect];
    rootdir             = 'E:\Documents\PhD\TempExportPietro\';
    csvwrite(sprintf('%s%s.csv',rootdir,mouseids{mou}),[Y_Response X],0) %writes matrix M to file filename as comma-separated values.

end

nModels = size(B,1);

%% Make figure of cross-validated performance on test data:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]); hold all;

Colors  = parula(length(mouseids));
xloc    = 1; %position of performance per mouse/mean

for iMod=1:nModels
    for iMou = 1:length(mouseids)
        h = plot(xloc,squeeze(cell2mat(CV_Perf(iMod,iMou,:))),'.','Color',Colors(iMou,:),'MarkerSize',40); hold all;
        xloc = xloc+1;
    end
    
    mPerf        = nanmean([CV_Perf{iMod,:,:}]);
    stdPerf      = nanstd([CV_Perf{iMod,:,:}]);
    plot(xloc,mPerf,'k.','MarkerSize',30)
    errorbar(xloc,mPerf,stdPerf,stdPerf,'k','LineWidth',2)
    xloc = xloc+2;
    
end

%Make up:
set(gca, 'XTick', 1:(length(mouseids)+2)*nModels,'XTickLabels', repmat([mouseids 'Mean' ' '],1,nModels),'XTickLabelRotation',45)
ylim([0.33 0.8])
xlim([0.5 (length(mouseids)+2)*nModels+0.5])
ylabel('Cross-validated % Correct')

%% Make figure of coefficients:
[~,nMice,~]  = size(B);
Modality_string = {'Auditory' 'Visual'};

for iMod=1:nModels
    nBeta       = length(X_string{iMod});
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.2 0.25*nBeta 0.65]); hold all;
    for iModality = 1:2
        for iB = 1:nBeta
            subplot(1,2*nBeta,iModality+(iB-1)*2)
            for iMou = 1:nMice
                allmodelcoef            = [B{iMod,iMou,:}];
                nSes = size(allmodelcoef,2)/2;
                plot(allmodelcoef(iB+1,iModality:2:nSes*2),repmat(iMou,1,nSes)/10,'.','Color',Colors(iMou,:),'MarkerSize',45); hold all;
            end
            set(gca, 'YTick', [1:nMice]/10,'YTickLabels', mouseids,'YTickLabelRotation',15,'Fontsize',10)
            ylim([0 (nMice+1)/10])
            title(sprintf('%s: weights on %s',X_string{iMod}{iB},Modality_string{iModality}),'Fontsize',15)
            plot([0,0],[0 (nMice+1)/10],'k:')
            xlim([-1 1])
        end
    end
end


end