function MOL_SessionRast(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};

set(0,'defaultAxesFontSize',20)

states          = {'iti'            'stim'         'respwin'      'timeout'};
statecolors     = {[0.95 0.95 0.95], [0.9 0.7 0.3], [0.2 0.9 0.4], [0.8 0.1 0.1]};

trialtypes      = {'A'          'V'          'P'            'C'         'Y'              'X'};
trialcolors     = {[1 0.2 0.2], [0.2 0.2 1], [0.5 0.5 0.5], [0.9 0 0.9] ,[0.95 0.5 0.5] ,[0.5 0.5 0.95]};

RewardColor     = [0.1 0.1 0.1];
AudioLickColor  = [1 0 0];
VisualLickColor = [0 0 1];

SortBy           = 'stimChange';
% SortBy          = 'stimStart'; %Sort the trials to get a clean overview
% SortBy          = 'audioInt'; %Sort the trials to get a clean overview

AlignOn         = 'stimChange'; %On which timestamp to align as t=0
% AlignOn         = 'stimStart'; %On which timestamp to align as t=0
% SplitBy         = {'trialType'};


for sesid = unique(sessionData.session_ID)'
    
    [temptrialData] = MOL_getTempPerSes(sesid,trialData);
        
    %% Sort the trials
    if isfield(trialData,SortBy)
        Delay = temptrialData.(SortBy) - temptrialData.trialStart;
        [~,trialorder] = sort(Delay);
    else
        trialorder = temptrialData.trialNum;
    end
    
    figure;
    set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
    
    for i = 1:length(trialorder)
        trial = trialorder(i);
%         rectangle('Position',[0 i temptrialData.trialEnd(trial)-temptrialData.trialStart(trial) 0.25],'FaceColor','k','EdgeColor','k'); hold all;
        
        for state = 1:length(states)
            if isfield(temptrialData,strcat(states{state},'Start')) && isfield(temptrialData,strcat(states{state},'End'))
                
                ts_start    = temptrialData.(strcat(states{state},'Start'))(trial) - temptrialData.(AlignOn)(trial);
                ts_end      = temptrialData.(strcat(states{state},'End'))(trial)  - temptrialData.(AlignOn)(trial);
                
                if ~isnan(ts_start) && ~isnan(ts_end)
                    if strcmp(states{state},'stim')
                        idx     = strcmp(trialtypes,temptrialData.trialType{trial});
                        rectangle('Position',[ts_start/1e6 i+0.25 (ts_end-ts_start)/1e6 0.5],'FaceColor',trialcolors{idx},'EdgeColor',trialcolors{idx}); hold all;
                    else
                        rectangle('Position',[ts_start/1e6 i+0.25 (ts_end-ts_start)/1e6 0.5],'FaceColor',statecolors{state},'EdgeColor',statecolors{state}); hold all;
                    end
                end
            end
        end
        
        if isfield(temptrialData,'stimChange') && ~isnan(temptrialData.stimChange(trial))
            rectangle('Position',[temptrialData.stimChange(trial)/1e6 - temptrialData.(AlignOn)(trial)/1e6 i 0.1 1],'FaceColor',[0.9 0.1 0.9],'EdgeColor',[0.1 0.1 0.1]); hold all;
        end
        
        if isfield(temptrialData,'rewardTime') && ~isnan(temptrialData.rewardTime(trial))
            rectangle('Position',[temptrialData.rewardTime(trial)/1e6-temptrialData.(AlignOn)(trial)/1e6 i 0.2 1],'FaceColor',RewardColor,'EdgeColor',RewardColor); hold all;
        end
        
        if isfield(temptrialData,'passiveRewardTime')
            if ~isnan(temptrialData.passiveRewardTime(trial))
                rectangle('Position',[temptrialData.passiveRewardTime(trial)/1e6-temptrialData.(AlignOn)(trial)/1e6 i 0.2 1],'FaceColor',RewardColor,'EdgeColor',RewardColor); hold all;
            end
        end
        
        for lick = 1:length(temptrialData.lickTime{trial})
            if sessionData.VisualLeftCorrectSide  
                RightLickColor = AudioLickColor; LeftLickColor = VisualLickColor;
            else RightLickColor = VisualLickColor; LeftLickColor = AudioLickColor;
            end
            licktime = temptrialData.lickTime{trial}(lick) - temptrialData.(AlignOn)(trial);
            if ~isnan(licktime) %Sometimes last trial of a session did not get TTL of stimChange
                if temptrialData.lickSide{trial}(lick) == 'R'
                    rectangle('Position',[licktime/1e6 i 0.005 1],'FaceColor',RightLickColor,'EdgeColor',RightLickColor,'LineWidth', 0.8); hold all;
                elseif temptrialData.lickSide{trial}(lick) == 'L'
                    rectangle('Position',[licktime/1e6 i 0.005 1],'FaceColor',LeftLickColor,'EdgeColor',LeftLickColor,'LineWidth', 0.8); hold all;
                end
            end
        end
    end
    
    p = patch([0 0.1 0.1 0],[1 1 i+1 i+1],[0.6 0.6 0.6]);
    set(p,'FaceAlpha',0.8,'EdgeColor','none');

    xmax = max(temptrialData.trialEnd - temptrialData.trialStart);
    if xmax>15; xmax = 15; end
    xlim([-8 2])
    ylim([1 i+1])
    ylabel('Trials')
    xlabel('Time relative to stimulus onset (s)')
    set(gca,'Xtick',-8:2:floor(xmax))
%     set(gca,'Xtick',2:2:floor(xmax),'Ytick',trialorder'+0.5,'yticklabel',cellstr(num2str(trialorder'))')
end


end
