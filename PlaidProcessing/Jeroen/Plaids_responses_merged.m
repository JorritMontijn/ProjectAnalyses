% This script makes a figure like in the Carandini paper
% dependencies: - plaids_responses_subfunction.m
%               - stimDetection.m

clear all

%% INPUT
filenamesPlaids = {'D:\Data\Processed\imagingdata\20130313\xyt01\20130313xyt01_ses.mat'};


%filenames orientation tuning files corresponding with plaid files
filenamesTuning = {'D:\Data\Processed\imagingdata\20130313\xyt02\20130313xyt02_ses.mat'};

%percentage of neurons with highest OSI to include in analysis
percNeuronsInclude = 35;


%% Calculate mean dF/F values for all neurons and all stimulus 'situations'
meanCalciumNeuronMatComb=[];

for s = 1 : length(filenamesPlaids)
	[meanCalciumNeuronMat, vecAngle] = Plaids_responses_subfunction(filenamesPlaids{s}, filenamesTuning{s}, percNeuronsInclude);
	
		
	%normalise values
	[r c] = size(meanCalciumNeuronMat);
	for i = 3 : c
		maxVal = max(meanCalciumNeuronMat(:,i));
		meanCalciumNeuronMat(:,i) = meanCalciumNeuronMat(:,i)/maxVal;
	end
	
	%merge matrices from different sessions
	meanCalciumNeuronMatComb = [meanCalciumNeuronMatComb; meanCalciumNeuronMat];
	
end

	

%replace negative dF/F values with 0 and outliers with 'NaN'
for i = 2*length(meanCalciumNeuronMatComb)+1 : numel(meanCalciumNeuronMatComb)
	if meanCalciumNeuronMatComb(i)<=0
		meanCalciumNeuronMatComb(i) = 0;
	elseif meanCalciumNeuronMatComb(i) >=2
		meanCalciumNeuronMatComb(i) = NaN;
	end
end


%merge 0-179 with 180-359 degrees
for i = 1 : length(meanCalciumNeuronMatComb)
	if meanCalciumNeuronMatComb(i,2) == 180
		meanCalciumNeuronMatComb(i,2) = 0;
	elseif meanCalciumNeuronMatComb(i,2) == 225
		meanCalciumNeuronMatComb(i,2) = 45;
	elseif meanCalciumNeuronMatComb(i,2) == 270
		meanCalciumNeuronMatComb(i,2) = 90;
	elseif meanCalciumNeuronMatComb(i,2) == 315
		meanCalciumNeuronMatComb(i,2) = 135;
	end
end


%Sort all neurons based on their preferred orientation and calculate mean
%dF/F for all neurons with the same pref orientation
meanMatAngles = cell(13,1);
for i = [3:12 14 15]
	neurs = length(meanCalciumNeuronMatComb);
	meanMat = vecAngle';
	meanMat = meanMat(1:4);
	
	
	for j = 1 : length(meanMat)
		angle = meanMat(j,1);
		newvec = [];
		
		
		
		for k = 1 : neurs
			if meanCalciumNeuronMatComb(k,2) == angle
				
				if i == 3 || i == 8	%combine situations 3 and 4, and situations 8 and 9 because they're the same
					newvec = [newvec meanCalciumNeuronMatComb(k, i) meanCalciumNeuronMatComb(k, i+1)];
				else
					newvec = [newvec meanCalciumNeuronMatComb(k, i)];
				end
				
			end
			
			%calculate mean dF/F value and the sem
			meanMat(j,2) = nanmean(newvec); %mean of mean dF/F values for all neurons with a certain preferred orientation
			meanMat(j,3) = nanstd(newvec)/sqrt(length(newvec)); %sem

			
		end
	end
	meanMatAngles{i-2, 1} = meanMat; %put in cell array
	
end


%% Make plot

figure

subplot(5,5,13)
errorbar(meanMatAngles{1,1}(:,1), meanMatAngles{1,1}(:,2), meanMatAngles{1,1}(:,3))
axis([-10 145 min(meanMatAngles{1,1}(:,2))-max(meanMatAngles{1,1}(:,3)) max(meanMatAngles{1,1}(:,2))+max(meanMatAngles{1,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])


subplot(5,5,14)
errorbar(meanMatAngles{5,1}(:,1), meanMatAngles{5,1}(:,2), meanMatAngles{5,1}(:,3))
axis([-10 145 min(meanMatAngles{5,1}(:,2))-max(meanMatAngles{5,1}(:,3)) max(meanMatAngles{5,1}(:,2))+max(meanMatAngles{5,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])

subplot(5,5,15)
errorbar(meanMatAngles{9,1}(:,1), meanMatAngles{9,1}(:,2), meanMatAngles{9,1}(:,3))
axis([-10 145 min(meanMatAngles{9,1}(:,2))-max(meanMatAngles{9,1}(:,3)) max(meanMatAngles{9,1}(:,2))+max(meanMatAngles{9,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])

subplot(5,5,18)
errorbar(meanMatAngles{3,1}(:,1), meanMatAngles{3,1}(:,2), meanMatAngles{3,1}(:,3))
axis([-10 145 min(meanMatAngles{3,1}(:,2))-max(meanMatAngles{3,1}(:,3)) max(meanMatAngles{3,1}(:,2))+max(meanMatAngles{3,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])

subplot(5,5,19)
errorbar(meanMatAngles{6,1}(:,1), meanMatAngles{6,1}(:,2), meanMatAngles{6,1}(:,3))
axis([-10 145 min(meanMatAngles{6,1}(:,2))-max(meanMatAngles{6,1}(:,3)) max(meanMatAngles{6,1}(:,2))+max(meanMatAngles{6,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])

subplot(5,5,20)
errorbar(meanMatAngles{10,1}(:,1), meanMatAngles{10,1}(:,2), meanMatAngles{10,1}(:,3))
axis([-10 145 min(meanMatAngles{10,1}(:,2))-max(meanMatAngles{10,1}(:,3)) max(meanMatAngles{10,1}(:,2))+max(meanMatAngles{10,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])

subplot(5,5,23)
errorbar(meanMatAngles{4,1}(:,1), meanMatAngles{4,1}(:,2), meanMatAngles{4,1}(:,3))
axis([-10 145 min(meanMatAngles{4,1}(:,2))-max(meanMatAngles{4,1}(:,3)) max(meanMatAngles{4,1}(:,2))+max(meanMatAngles{4,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])

subplot(5,5,24)
errorbar(meanMatAngles{8,1}(:,1), meanMatAngles{8,1}(:,2), meanMatAngles{8,1}(:,3))
axis([-10 145 min(meanMatAngles{8,1}(:,2))-max(meanMatAngles{8,1}(:,3)) max(meanMatAngles{8,1}(:,2))+max(meanMatAngles{8,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])

subplot(5,5,22)
errorbar(meanMatAngles{12,1}(:,1), meanMatAngles{12,1}(:,2), meanMatAngles{12,1}(:,3))
axis([-10 145 min(meanMatAngles{12,1}(:,2))-max(meanMatAngles{12,1}(:,3)) max(meanMatAngles{12,1}(:,2))+max(meanMatAngles{12,1}(:,3))])
xlabel('Preferred orientation')
ylabel('Normalised response')
set(gca, 'XTick', [0 45 90 135])

subplot(5,5,10)
errorbar(meanMatAngles{13,1}(:,1), meanMatAngles{13,1}(:,2), meanMatAngles{13,1}(:,3))
axis([-10 145 min(meanMatAngles{13,1}(:,2))-max(meanMatAngles{13,1}(:,3)) max(meanMatAngles{13,1}(:,2))+max(meanMatAngles{13,1}(:,3))])
set(gca, 'XTick', [0 45 90 135])



%% Images

[X1] = imread('90_50.png');
subplot(5,5,11)
subimage(X1)
axis off

[X2] = imread('90_25.png');
subplot(5,5,16)
subimage(X2)
axis off

[X3] = imread('0_50.png');
subplot(5,5,3)
subimage(X3)
axis off

[X4] = imread('0_25.png');
subplot(5,5,4)
subimage(X4)
axis off

[X5] = imread('0_0.png');
subplot(5,5,5)
subimage(X5)
axis off

[X6] = imread('90_0.png');
subplot(5,5,21)
subimage(X6)
axis off

[X7] = imread('0_100.png');
subplot(5,5,2)
subimage(X7)
axis off

[X8] = imread('90_100.png');
subplot(5,5,6)
subimage(X8)
axis off


set(gcf,'color','w')