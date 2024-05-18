%uses files saved by: decodeDrift_direction.m, decodeDriftControl_direction.m, decodeDrift3T_direction.m
close all
clear all
zeroMean = 0;
deriv=0; 

% % For reviewers (3A for direction decoding):
% disjoint=0;%for expNum=2  
% illusion = 0;%0=direction, 1=illusion
% expNum=1;
% crossValidation = 'run';%'trial','run'
% locThresh = 0.2; %only include voxels with coherence>locThresh, for each of the two localizers.


% % %PANEL A left (loc 1):
% disjoint=0;%for expNum=2  
% illusion = 1;%0=direction, 1=illusion
% expNum=1;
% crossValidation = 'run';%'trial','run'
% locThresh = 0.0; %only include voxels with coherence>locThresh, for each of the two localizers.

% %PANEL A right (loc 1):
% disjoint=0;%for expNum=2   
% illusion = 1;%0=direction, 1=illusion
% expNum=2;
% crossValidation = 'run';%'trial','run'
% locThresh = 0.2; %only include voxels with coherence>locThresh, for each of the two localizers.

% %PANEL B:
% disjoint=0;%for expNum=2
% illusion = 0;%0=direction, 1=illusion
% expNum=2;
% crossValidation = 'run';%'trial','run'
% locThresh = 0.2; %only include voxels with coherence>locThresh, for each of the two localizers.

%PANEL C:
disjoint=0;%for expNum=2  
illusion = 0;%0=direction, 1=illusion
expNum=3;
crossValidation = 'run';%'trial','run'
locThresh = 0.2; %only include voxels with coherence>locThresh, for each of the two localizers.

signifThresh1=95;
signifThresh2=99;

meanColor = [0 0 0];
subColor = [0.6 0.6 0.6];

plotColors = {[0, 0.4470, 0.7410], [0.6350, 0.0780, 0.1840],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.8500, 0.3250, 0.0980]};	

nanColor = [0 0 0];
disjointStr = '';
if disjoint
    disjointStr = 'Disjoint';
end
saveFolder = '~/noah/';
switch expNum
    case 1
        expName = '3conds';
    case 2
        expName = 'attn';
    case 3
        expName = '4conds';
    case 4
        expName = '3tdata';
end
thrsh = num2str(locThresh,'%.2f');
locStr = thrsh([1 3:4]);
directionOrIllusionStr = 'direction';
if illusion
    directionOrIllusionStr = 'illusion';
end

%load the appropriate file
if expNum<3
load([saveFolder 'decodeDrift' disjointStr  '_' directionOrIllusionStr '_' expName '_' crossValidation '_thresh' locStr '.mat'],'crossValidation',...
    'r2thresh','locThresh','ROIs','subNames','dataFolder',...
    'subRoiNans','subNumTrials','subNumScans',...
    'subAcc','subRandAcc','subRand95Acc','nperms','locSharedVox','totalVox','locNumVox');
else
    load([saveFolder 'decodeDrift_direction_' expName '_' crossValidation '.mat'],'crossValidation',...
        'r2thresh','locThresh','ROIs','subNames','dataFolder',...
        'subRoiNans','subNumTrials','subNumScans',...
        'subAcc','subRandAcc','subRand95Acc','nperms','numVox','totalVox');
    locStr = '';
end

numSubs = length(subNames);

if expNum==2
    goodSubs = [2:10];
end
if expNum==1
    goodSubs = [1:3 5:12];%keeping only one session of subject 47
end
if expNum==1 | expNum==2
    numSubs = length(goodSubs);
    subAcc = subAcc(goodSubs,:,:,:);
    subRandAcc = subRandAcc(goodSubs,:,:,:,:);
    subRand95Acc = subRand95Acc(goodSubs,:,:,:);
    locSharedVox = locSharedVox(goodSubs,:);
    totalVox = totalVox(goodSubs,:);
    locNumVox = locNumVox(goodSubs,:,:);
    subNumTrials = subNumTrials(2:end);
    subNumScans = subNumScans(goodSubs);
    subRoiNans = subRoiNans(goodSubs,:);
    
    subNames = subNames(goodSubs);
end




if illusion %average over left vs no illusion and right vs no illusion
    subAcc = squeeze(mean(subAcc,4));
    subRandAcc = squeeze(mean(subRandAcc,4));
    subRand95Acc= squeeze(mean(subRand95Acc,4));
    
%     subAcc = squeeze(subAcc(:,:,:,1));
%     subRandAcc = squeeze(subRandAcc(:,:,:,1,:));
%     subRand95Acc= squeeze(subRand95Acc(:,:,:,1));
end

numRois = length(ROIs);

facealpha = 0.3;
ilocTitle = {'stim','eye'};
randDist = squeeze(mean(subRandAcc,1));
meanAcc = squeeze(mean(subAcc,1));
stdAcc = squeeze(std(subAcc,1));
numLocs=2;
if expNum==4
   meanAcc = meanAcc'; 
   numLocs=1;
   clear randDist
   randDist(:,1,:) = squeeze(mean(subRandAcc,1));
end

for iroi=1:numRois
    for iloc=1:numLocs
%         decodingAcc(isub,iroi,iloc) = res{iloc}.accuracy;
%         randAccDist(isub,iroi,iloc,:) = res{iloc}.randAccuracy;
        
%         temp = subRandAcc(:,iroi,iloc,:);
%         randDist(iroi,iloc,:) = temp(:);
        
%         meanAcc(iroi,iloc) = mean(decodingAcc(:,iroi,iloc));
        rand99(iroi,iloc) = prctile(randDist(iroi,iloc,:),signifThresh2);
        rand95(iroi,iloc) = prctile(randDist(iroi,iloc,:),signifThresh1);
        rand50(iroi,iloc) = prctile(randDist(iroi,iloc,:),50);
        rand05(iroi,iloc) = prctile(randDist(iroi,iloc,:),100-signifThresh1);
        
        %get p-value for each ROI
        pval(iroi,iloc) = 1 - sum(randDist(iroi,iloc,:)<meanAcc(iroi,iloc))/nperms;
    end
end

maxRandAcc = max(randDist,[],3);
minRandAcc = min(randDist,[],3);
roiInd = 1:numRois;

% figure(1)
% rows=2;
% cols=1;
% badRoi= zeros(2,numRois);
% 
% for iloc=1:numLocs
%     for iroi=1:numRois
%        if sum(isnan( subAcc(:,iroi,iloc)))>0
%           %mark this ROI
%           badRoi(iloc,iroi) = 1;
%           randDist(iroi,iloc,:) = randDist(1,iloc,:);
%           facecolorVec(iloc,iroi,:) = nanColor;
%        else
%            facecolorVec(iloc,iroi,:) = plotColors{iloc};
%        end
%     end
% end
% 
% 
% figure(2)
% % roiInd = 1:numRois;
% for iloc=1:numLocs
%     shift = (iloc-1.5)*0.2;
%    temp = plot([ roiInd+shift; roiInd+shift], [minRandAcc(:,iloc)'; maxRandAcc(:,iloc)'], 'linewidth',3,'color',squeeze(facecolorVec(iloc,1,:)).^0.5);
%    p(iloc) = temp(1);%for legend, keep plot properties of first ROI
%    hold on
%    plot([ roiInd+shift-0.3; roiInd+shift+0.3], [rand95(:,iloc)'; rand95(:,iloc)'], 'linewidth',3,'color',squeeze(facecolorVec(iloc,1,:)).^0.5);
%    plot([ roiInd+shift-0.3; roiInd+shift+0.3], [rand05(:,iloc)'; rand05(:,iloc)'], 'linewidth',3,'color',squeeze(facecolorVec(iloc,1,:)).^0.5);
%    p(3) = plot(roiInd+shift, meanAcc(:,iloc), '.','markersize',15,'color',[0 0 0]);
% end
% xlim([0 numRois+1]);
% ylim([0.3 0.7]);
% xlabel('ROI');
% xticks([1:numRois])
% xticklabels(ROIs);
% ylabel('decoding accuracy');
% if expNum==1 | expNum==2
%     L=legend([p(1) p(2) p(3)],'stim','eye','data');
% elseif expNum==3
%     L=legend([p(1) p(2) p(3)],'local+global','local','data');
% end
% title([expName, ' ', directionOrIllusionStr, ' ', crossValidation ' cv, ' num2str(size(subRandAcc,1)) ' subjects, co>' num2str(locThresh)]);
% set(gcf,'position',[100 1500 1200 400]);



% roiInd = 1:numRois/2;
% roiInd = [1:12 25:2:length(ROIs)];
roiInd = [25:28];
hemiROIs = ROIs(roiInd);
for i=1:length(hemiROIs)
   hemiROIs{i} =  hemiROIs{i}(2:end);
end

hemiROIs = {'EVC','LO','hMT+','V3A/B'};

roiPlotInd = 1:length(roiInd);
% hemiROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01','lV02'};
% hemiROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01','lV02','lEVC','lLO','lMT+','lV3AB'};
for iloc=1:2
    figure(2+iloc); clf;
    shift = 0;%(iloc-1.5)*0.2;
   plot([ roiPlotInd+shift; roiPlotInd+shift], [minRandAcc(roiInd,iloc)'; maxRandAcc(roiInd,iloc)'], 'linewidth',3,'color',subColor);
   hold on
   plot([ roiPlotInd+shift-0.3; roiPlotInd+shift+0.3], [rand95(roiInd,iloc)'; rand95(roiInd,iloc)'], ':','linewidth',2,'color',subColor);
   plot([ roiPlotInd+shift-0.3; roiPlotInd+shift+0.3], [rand50(roiInd,iloc)'; rand50(roiInd,iloc)'], '-','linewidth',3,'color',subColor);
%    roiX = repmat(roiInd+shift,numSubs,1);

   %plot single subject accuracy
%    plot(roiPlotInd,subAcc(:,roiInd,iloc),'o','markersize',5,'color',plotColors{7});

    %plot decoding accuracy for ROIs with p<0.05
   signifRoiInd = meanAcc(roiInd,iloc)>rand95(roiInd,iloc);%1=significant ROI, 0=insignificant ROI
   plot(roiPlotInd(signifRoiInd)+shift, meanAcc(roiInd(signifRoiInd),iloc), '.','markersize',40,'color',plotColors{7});
   %plot decoding accuracy for ROIs with p<0.01
   signifRoiInd2 = meanAcc(roiInd,iloc)>rand99(roiInd,iloc);%1=significant ROI, 0=insignificant ROI
   plot(roiPlotInd(signifRoiInd2)+shift, meanAcc(roiInd(signifRoiInd2),iloc), '.','markersize',40,'color',plotColors{2});
   %plot decoding accuracy for insignificant ROIs
   plot(roiPlotInd(~signifRoiInd)+shift, meanAcc(roiInd(~signifRoiInd),iloc), '.','markersize',40,'color',[0.4 0.4 0.4]);
   
   xlim([0.5 length(roiPlotInd)+0.5]);
   ylim([0.35 0.7]);
%    xlabel('Left hemisphere ROI');
   xticks([roiPlotInd])
   xticklabels(hemiROIs);
   yticks(0.4:0.1:0.6);
%    ylabel({'Decoding accuracy (% correct)'});
   % L=legend([p(1) p(2) p(3)],'stim','eye','data');
   title([expName, ' ', directionOrIllusionStr(1:3), ' ', crossValidation ', ' num2str(size(subRandAcc,1)) 'subs, co>' num2str(locThresh) ', loc ', num2str(iloc)]);
   fontsize=16;
   drawPublishAxis('xLabelOffset', -6/64,'yLabelOffset', -10/64, 'xAxisMargin', 0/64, 'yAxisMargin', 0/64,'labelFontSize',fontsize);
   set(gcf,'position',[5 5 10 10]);
   legend off
   
   print('-painters','-dpdf',fullfile(saveFolder,'figures/',[expName, '_', directionOrIllusionStr, '_loc', disjointStr, '_', num2str(iloc) zeroMeanStr '.pdf']));
end

%% output results to copy into table
[meanAcc(roiInd,1)'; stdAcc(roiInd,1)']
pval(roiInd,1)'

if expNum==3
    [meanAcc(roiInd,2)'; stdAcc(roiInd,2)']
    pval(roiInd,2)'
end

%% NUMBER OF VOXELS
if expNum<3
    roiTotalVox = totalVox(:,roiInd);%all voxels
    roiSharedVox = locSharedVox(:,roiInd);%
    for iloc=1:2
        roiNumVox{iloc} = locNumVox(:,roiInd,iloc);%chosen for each localizer (after R2).
    end
else
    roiTotalVox = totalVox(:,roiInd);%all voxels
    roiNumVox = numVox(:,roiInd);%chosen by R2.
end
