function [visconditions,auconditions,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMatOpto(varargin)
%% Get input arguments:
tempsessionData     = varargin{1};
temptrialData       = varargin{2};

%% Get the trialtypes per session:
trialtypes = unique(temptrialData.trialType);

%% Define probe, auditory and visual trials and associated settings:
seshasprobe = 0; seshasau = 0; seshasvis = 0; seshasconflict = 0;
for tt = 1:length(trialtypes)
    switch trialtypes{tt}
        case 'P'
            seshasprobe         = 1;
            probetrialtype      = 'P';
        case 'Q'
            seshasprobe         = 1;
            probetrialtype      = 'Q';
        case 'R'
            seshasprobe         = 1;
            probetrialtype      = 'R';
        case 'V'
            if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                seshasvis           = 1;
                vistrialtype        = 'V';
                vissessionfield     = 'vecVisStimIntensities';
                vistrialfield       = 'visualInt';
            elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                seshasvis           = 1;
                vistrialtype        = 'V';
                vissessionfield     = 'vecOriChange';
                vistrialfield       = 'visualOriChange';
            end
        case 'X'
            seshasvis           = 1;
            vistrialtype        = 'X';
            vissessionfield     = 'vecOriChange';
            vistrialfield       = 'visualOriChange';
        case 'A'
            if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                seshasau            = 1;
                autrialtype         = 'A';
                ausessionfield      = 'vecAuStimIntensities';
                autrialfield        = 'audioInt';
            elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                seshasau            = 1;
                autrialtype         = 'A';
                if strcmp(tempsessionData.auChangeUnit,'Hz') || isempty(tempsessionData.auChangeUnit{1})
                    ausessionfield      = 'vecFreqChange';
                    autrialfield        = 'audioFreqChange';
                else strcmp(tempsessionData.auChangeUnit,'Oct')
                    ausessionfield      = 'vecOctChange';
                    autrialfield        = 'audioOctChange';
                end
            end
        case 'Y'
            seshasau            = 1;
            autrialtype         = 'Y';
            if strcmp(tempsessionData.auChangeUnit(1),'Hz') || isempty(tempsessionData.auChangeUnit{1})
                ausessionfield      = 'vecFreqChange';
                autrialfield        = 'audioFreqChange';
            else strcmp(tempsessionData.auChangeUnit(1),'Oct')
                ausessionfield      = 'vecOctChange';
                autrialfield        = 'audioOctChange';
            end
        case 'C'
            seshasconflict      = 1;
            conflicttrialtype   = 'C';
    end
end

%% Visual conditions
if seshasvis
    if iscell(tempsessionData.(vissessionfield))
        visconditions   = tempsessionData.(vissessionfield){1};
    else
        visconditions   = tempsessionData.(vissessionfield);
    end
end

%% Auditory conditions:
if seshasau
    if iscell(tempsessionData.(ausessionfield))
        auconditions    = tempsessionData.(ausessionfield){1};
    else
        auconditions    = tempsessionData.(ausessionfield);
    end
end

%% Opto conditions:
if isfield(tempsessionData,'UseOpto')
    seshasopto = tempsessionData.UseOpto;
else         seshasopto = 0;
end

selectionopto = [];
if seshasopto && isfield(temptrialData,'hasphotostim') && any(temptrialData.hasphotostim==1)
    selectionopto(:,1) = ~(temptrialData.hasphotostim==1);
    selectionopto(:,2) = (temptrialData.hasphotostim==1);
else selectionopto = ones(length(temptrialData.trialStart),1);
end

%Init matrices to store output:
FullRespMat     = NaN(length(auconditions)+1,length(visconditions)+1,size(selectionopto,2),3);
FullnTrialsMat  = NaN(length(auconditions)+1,length(visconditions)+1,size(selectionopto,2));

%Init matrices to store output (without opto):
% FullRespMat     = NaN(length(auconditions)+1,length(visconditions)+1,3);
% FullnTrialsMat  = NaN(length(auconditions)+1,length(visconditions)+1);

%% Selection of probe trials:
if seshasprobe
    selectionprobe      = strcmp(temptrialData.trialType,probetrialtype);
else selectionprobe     = strcmp(temptrialData.trialType,'Z'); %Nonsensical index
end

%% Get responses:
for iOpto = 1:size(selectionopto,2)
    if seshasprobe
        FullnTrialsMat(1,1,iOpto)         = sum(selectionprobe & selectionopto(:,iOpto));
        for iResp = 1:3
            FullRespMat(1,1,iOpto,iResp)          = sum(temptrialData.vecResponse==iResp & selectionprobe & selectionopto(:,iOpto));
        end
    end
    
    if seshasau % Auditory only trials:
        for iAu = 1:length(auconditions)
            idx = strcmp(temptrialData.trialType,autrialtype) & ismember(abs(temptrialData.(autrialfield)),auconditions(iAu)) & selectionopto(:,iOpto);
            FullnTrialsMat(iAu+1,1,iOpto)    = sum(idx);
            for iResp = 1:3
                FullRespMat(iAu+1,1,iOpto,iResp)          = sum(temptrialData.vecResponse==iResp & idx);
            end
        end
    end
    
    if seshasvis % Visual only trials:
        for iVis = 1:length(visconditions)
            idx = strcmp(temptrialData.trialType,vistrialtype) & ismember(abs(temptrialData.(vistrialfield)),visconditions(iVis)) & selectionopto(:,iOpto);
            FullnTrialsMat(1,iVis+1,iOpto)    = sum(idx);
            for iResp = 1:3
                FullRespMat(1,iVis+1,iOpto,iResp)          = sum(temptrialData.vecResponse==iResp & idx);
            end
        end
    end
    
    if seshasconflict
        for iAu = 1:length(auconditions)
            for iVis = 1:length(visconditions)
                %Get index of trials as crosssection of both au and vis conditions:
                idx = ismember(abs(temptrialData.(autrialfield)),auconditions(iAu))...
                    & ismember(abs(temptrialData.(vistrialfield)),visconditions(iVis))...
                    & strcmp(temptrialData.trialType,conflicttrialtype)...
                    & selectionopto(:,iOpto);
                FullnTrialsMat(iAu+1,iVis+1,iOpto)    = sum(idx);
                for iResp = 1:3
                    FullRespMat(iAu+1,iVis+1,iOpto,iResp)          = sum(temptrialData.vecResponse==iResp & idx);
                end
            end
        end
    end
end

%% Checks and balances:
if FullnTrialsMat ~= sum(FullRespMat,4)
    error('trials dont match')
end

if nansum(FullRespMat(:)) ~= length(temptrialData.session_ID)
    error('total amount of trials dont match')
end

FullRespMat             = FullRespMat ./ repmat(FullnTrialsMat,1,1,1,3); %Make fractional (as percentage correct etc.)

end