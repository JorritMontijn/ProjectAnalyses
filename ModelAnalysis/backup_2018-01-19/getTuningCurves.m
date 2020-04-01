function sOut = getTuningCurves(matResp,vecStimTypeList)
    %Get tuning curves for neurons
    %

   %neuron x trial
   intNeurons = size(matResp,1);
   
   %get stimulus response by repetition; from [N x T] to [N x S x R]
   matRespNSR = getStimulusResponses(matResp,vecStimTypeList);
   
   %fit von Mises
   
   %get bandwidth, OSI, etc
   
   
end

