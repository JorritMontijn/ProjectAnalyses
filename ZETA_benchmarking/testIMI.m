%{
2.4 Fitting the IMI Model

We fit the m-IMI and TRRP models to spike trains derived from the LIF model
via maximum likelihood. 

2.4.1 Fitting the m-IMI Model

We follow Kass and Ventura (2001) to fit the m-IMI model to data. For
fitting m-IMI model in equation 2.4, we use the additive form: 
log ?(t, s?(t)) = log ?1(t) + log g1(t ? s?(t)).
(2.16)

We first represent a spike train as a binary sequence of 0s and 1s by
discretizing time into small intervals of length ?, letting 1 indicate that
a spike occurred within the corresponding time interval. We represent log
?1(t) and log g 1(t ? s *(t)) with cubic splines. Given suitable knots for
both terms, cubic splines may be described by linear combinations of
B-spline basis functions (de Boor, 2001),     
log?1(i?)=?k=1M?kAk(i?),
(2.17)

logg1(i??s?(i?))=?k=1L?kBk(i??s?(i?)),
(2.18)

where M and L are the numbers of basis functions that are determined by the
order of splines and the number of knots. Note that the shapes of B-spline
basis functions {Ak} and {Bk} also depend on the location of knots. Fitting
of the model is accomplished easily via maximum likelihood: for fixed
knots, the model is binary generalized linear model (McCullagh & Nelder,
1989) with     
log?(i?,s?(i?))=?k=1M?kAk(i?)+?k=1L?kBk(i??s?(i?)),
(2.19)

where {Ak(i?)} and {Bk(i?)} play the role of explanatory variables. The
coefficients of the spline basis elements, {?k} and {?k}, are determined
via maximum likelihood. This can be performed by using a standard software
such as R and Matlab Statistics Toolbox.   

In the following simulations we chose the knots by preliminary examination
of data. We conducted the fitting procedure for several candidates of knots
and then chose the one that minimizes the KL divergence between the
estimate and the LIF model.   
    
%}





%% select all neurons in LP and drifting grating stimuli
[sAggStim,sAggNeuron]=loadDataNpx('lateral posterior nucleus','driftinggrating');

%% get data for neuron #16 (original neuron #30)
intNeuronToAnalyze = 16;
sNeuron = sAggNeuron(intNeuronToAnalyze);
vecSpikeTimes = sNeuron.SpikeTimes;
vecStimOnTime = sAggStim(1).cellStim{1}.structEP.vecStimOnTime;
vecStimOffTime = sAggStim(1).cellStim{1}.structEP.vecStimOffTime;
dblTrialDur = min(vecStimOffTime-vecStimOnTime);
vecOrientation = sAggStim(1).cellStim{1}.structEP.Orientation;

%jitter
%vecStimOnTime = vecStimOnTime + 2*(rand(size(vecStimOnTime))-0.5)*dblTrialDur*2;
vecStimOnTime = sort(vecStimOnTime);

%% fit
tic
[dblMIMI_P,vecLatencies,sMIMI,sRate] = getMIMI(vecSpikeTimes,vecStimOnTime,dblTrialDur,4,2);
toc
return
%% plot
subplot(2,3,3)
	plot(sMIMI.vecX,sMIMI.vecY);
	%errorbar(vecX,vecY,vecY_sem);
	hold on
	plot(sMIMI.vecX,sMIMI.vecFitY,'k--');
	hold off
xlabel('Time after stimulus onset (s)');
ylabel('Mean firing rate (Hz)');
legend('Data','MIMI model fit');

%plot mean rates
subplot(2,3,4)
scatter(sort(sMIMI.vecY),sort(sMIMI.vecFitY))
vecLim = [min([get(gca,'xlim') get(gca,'ylim')]) max([get(gca,'xlim') get(gca,'ylim')])];
xlim(vecLim);ylim(vecLim);
hold on
plot(vecLim,vecLim,'--','color',[0.5 0.5 0.5])

%% transform to probabilities from poisson process
subplot(2,3,5)
vecResid = sMIMI.residual;
vecResidNoMod = sMIMI.vecY-mean(sMIMI.vecY);
histx(sMIMI.vecFitY)
title(sprintf('mean fit y=%.3f,sd fit y=%.3f',mean(sMIMI.vecFitY),std(sMIMI.vecFitY)))

%% are fitted values normal?
%k-s test
[h,pKS,ksstat,cv] = kstest(zscore(sMIMI.vecFitY));


subplot(2,3,6)
vecQ = normcdf(sMIMI.vecFitY,mean(sMIMI.vecFitY));

hold on
scatter(linspace(0,1,numel(vecQ)),sort(vecQ))

vecLim = [0 1];
xlim(vecLim);ylim(vecLim);
plot(vecLim,vecLim,'--','color',[0.5 0.5 0.5]);
hold off
title(sprintf('K-S test,p=%.3f',pKS))
fixfig;