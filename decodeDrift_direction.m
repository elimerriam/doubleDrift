% decodeDrift_direction.m
%
% associated with the following publication: Steinberg, NJ, Roth, ZN, Movshon, JA, and Merriam, EP (2024).
% Brain representations of motion and position in the double-drift illusion
% 
% DOI:https://doi.org/10.7554/eLife.76803
%
%   usage: decodeDrift_direction()
%   by: zvi roth
%   date: 3/05/2024
%   purpose: Decode illusion direction in Expt 1 and 2.

close all
clear all

saveFolder = '~/data/'; %location of saved timeseries

nperms = 1000; %number of permutations for bootstrapping

zeroMean = 1; %0 = no normalization. 1 = subtract mean.

%experiment number
expNum=1; %Expt 1

%leave-one-out cross-validation: leave out one run, or one trial?
crossValidation = 'run';%'trial','run'

tic
r2thresh = 50;%r2 percentile. Only use voxels with R^2 above this percentile.
locThresh = 0.20; %threshold for localizar coherence

thrsh = num2str(locThresh,'%.2f');
locStr = thrsh([1 3:4]);


switch expNum
    case 1
        expName = '3conds';
    case 2
        expName = 'attn';
end

load([saveFolder 'driftTS_' expName '.mat'],'ROIs','subNames','dataFolder',...
    'stimLoc','eyeLoc','driftBetas','subRoiNans','subNumTrials','subNumScans');

numSubs = length(subNames);
numRois = length(ROIs);

%create new ROIs by combining Benson atlas ROIs
extraROIs = {'lEVC','lLO','lMT+','lV3AB','rEVC','rLO','rMT+','rV3AB'};
extraROIind = {[1 2 3],  [5 6], [7 8], [9 10], [13 14 15], [17 18], [19 20],  [21 22]};
for isub=1:numSubs
    for xroi=1:length(extraROIs)
        combinedROIs = [];
        for roiPart=1:length(extraROIind{xroi})
            combinedROIs = [combinedROIs ' + ' ROIs{extraROIind{xroi}(roiPart)}];
        end
        ['new ' extraROIs{xroi} ' = ' combinedROIs]
        ROIs{numRois+xroi} = extraROIs{xroi};
        subRoiNans{isub,numRois+xroi} = [];
        tempStruct.scm = [];%[nvox 96]
        tempStruct.hdrlen = driftBetas{1,1}.hdrlen;%1
        tempStruct.nhdr = driftBetas{1,1}.nhdr;%96
        tempStruct.dim = [];%[2 2 nvox time]
        tempStruct.volumes = [];%[1 time]
        tempStruct.ehdr = [];%[nvox 96]
        tempStruct.ehdrste = [];%[nvox 96]
        tempStruct.r2 = [];%[nvox 1]
        tempStruct.covar = [];%[96 96]
        %localizer coherence values
        stimLoc{isub,numRois+xroi}.co = [];
        eyeLoc{isub,numRois+xroi}.co = [];
        for roiPart=1:length(extraROIind{xroi})
            curRoi = extraROIind{xroi}(roiPart);
            subRoiNans{isub,numRois+xroi} = [subRoiNans{isub,numRois+xroi}; subRoiNans{isub,curRoi}];
            tempStruct.scm = [tempStruct.scm; driftBetas{isub,curRoi}.scm];
            tempStruct.dim = cat(3,tempStruct.dim,driftBetas{isub,curRoi}.dim);
            tempStruct.volumes = [tempStruct.volumes driftBetas{isub,curRoi}.volumes];
            tempStruct.ehdr = [ tempStruct.ehdr; driftBetas{isub,curRoi}.ehdr];
            tempStruct.ehdrste = [ tempStruct.ehdrste; driftBetas{isub,curRoi}.ehdrste];
            tempStruct.r2 = [tempStruct.r2; driftBetas{isub,curRoi}.r2];
            tempStruct.covar = driftBetas{isub,curRoi}.covar;%just keeping the last ROI
            stimLoc{isub,numRois+xroi}.co = [stimLoc{isub,numRois+xroi}.co; stimLoc{isub,curRoi}.co];
            eyeLoc{isub,numRois+xroi}.co = [eyeLoc{isub,numRois+xroi}.co; eyeLoc{isub,curRoi}.co];
        end
        driftBetas{isub,numRois+xroi} = tempStruct;
    end
end

numRois = length(ROIs);%now including combined ROIs
clear locSharedVox totalVox locNumVox
for isub=1:numSubs
    nTrials = subNumTrials(isub);
    for iroi=1:numRois
        
        roiname = ROIs{iroi};
        roiNans = subRoiNans{isub,iroi};
        trialReg = driftBetas{isub,iroi};
        if ~isempty(trialReg)
            temp = trialReg.r2>prctile(trialReg.r2,50) & stimLoc{isub,iroi}.co(~roiNans)>locThresh & eyeLoc{isub,iroi}.co(~roiNans)>locThresh;
            locSharedVox(isub,iroi) = sum(temp);
            totalVox(isub,iroi) = length(temp);
            for iloc=1:2
                switch iloc
                    case 1
                        locCorrData = stimLoc{isub,iroi};
                    case 2
                        locCorrData = eyeLoc{isub,iroi};
                end
                %     keyboard
                locCorrData.co = locCorrData.co(~roiNans);%data from concatenation excluded these voxels
                goodVox = trialReg.r2>prctile(trialReg.r2,50) & locCorrData.co>locThresh;
                
                locNumVox(isub,iroi,iloc) = sum(goodVox);
                ['sub' num2str(isub) ' ' roiname ' loc ' num2str(iloc) ': using ' num2str(locNumVox(isub,iroi,iloc)) ' out of ' num2str(totalVox(isub,iroi)) ' voxels']
                if sum(goodVox)>1
                    
                    betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials*2)+1:nTrials*3])'; %this is for 1 vs 3
                    betas = zscore(betas,0,1);%zscore each trial
                  
                    if zeroMean %subtract mean from each voxel
                        betas = betas-mean(betas,2);
                    end
                    % create a grouping variable
                    group = cat(1, repmat('l', nTrials, 1), repmat('r', nTrials, 1));
                    
                    if strcmp(crossValidation, 'trial')
                        CVO = cvpartition(group, 'KFold', nTrials*2);
                        err = zeros(CVO.NumTestSets,1);
                        for i = 1:CVO.NumTestSets
                            trIdx = CVO.training(i);
                            teIdx = CVO.test(i);
                            ytest = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
                            err(i) = sum(~strcmp(ytest,group(teIdx)));
                        end
                        accuracy = 1 - sum(err)/sum(CVO.TestSize);
                        disp(sprintf('(plotROI): Classification accuracy: %f', accuracy));
                        
                    elseif strcmp(crossValidation, 'run')
                        
                        clear accuracy; clear ytest;
                        nRuns = length(group)/2/4;
                        runVec = repmat(1:nRuns, 4, 1);
                        runVec = runVec(:);
                        runVec = cat(1, runVec, runVec);
                        for iRun = 1:nRuns
                            trIdx = runVec ~= iRun;
                            teIdx = runVec == iRun;
                            ytest(teIdx) = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
                        end
                        accuracy = sum(ytest==group') / length(group);
                        disp(sprintf('(plotROI): Classification accuracy: %f', accuracy));
                    else
                        keyboard
                    end
                    
                    nScans = subNumScans(isub);
                    
                    % trialsPerCond will be the trials of BOTH CONDITIONS
                    trialsPerCond = length(group) / nScans;
                    groupByScan = cat(1, repmat('l', trialsPerCond/2, 1), repmat('r', trialsPerCond/2, 1));
                    
                    % We have to create 'length(nScans)' amount of randomized labels
                    clear randLabels
                    for iBoot = 1:nperms
                        for iScan = 1:nScans
                            tempList = groupByScan(randperm(size(groupByScan,1)));
                            idx1 = (iScan*4)-3;
                            idx2 = iScan*4;
                            randLabels(idx1:idx2,iBoot) = tempList(1:4);
                            randLabels(idx1+length(group)/2 : idx2+length(group)/2,iBoot) = tempList(5:8);
                        end
                    end
                    disp(sprintf('(plotROI): Classification method = %swise', crossValidation));
                    bootAccuracy = zeros(nperms, 1);
                    
                    if strcmp(crossValidation, 'trial')
                        for iBoot=1:nperms
                            gNull = randLabels(:,iBoot);
                            for i = 1:CVO.NumTestSets
                                trIdx = CVO.training(i);
                                teIdx = CVO.test(i);
                                ytest = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
                                err(i) = sum(~strcmp(ytest,gNull(teIdx)));
                            end
                            bootAccuracy(iBoot) = 1 - sum(err)/sum(CVO.TestSize);
                        end
                    elseif strcmp(crossValidation, 'run')
                        for iBoot=1:nperms
                            gNull = randLabels(:,iBoot);
                            clear ytest;
                            for iRun = 1:nRuns
                                trIdx = runVec ~= iRun;
                                teIdx = runVec == iRun;
                                ytest(teIdx) = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
                            end
                            bootAccuracy(iBoot) = sum(ytest==group') / length(group);
                        end
                    else
                        keyboard
                    end
                    
                    disp(sprintf('(plotROI): Bootstrapped threshold is: %f', prctile(bootAccuracy, 95)));
                    
                    subAcc(isub,iroi,iloc) = accuracy;
                    subRandAcc(isub,iroi,iloc,:) = bootAccuracy;
                    subRand95Acc(isub,iroi,iloc) = prctile(bootAccuracy, 95);
                else
                    'ZERO GOOD VOXELS!!'
                    subAcc(isub,iroi,iloc) = NaN;
                    subRandAcc(isub,iroi,iloc,:) = NaN(nperms,1);
                    subRand95Acc(isub,iroi,iloc) = NaN;%prctile(bootAccuracy, 95);
                end
            end
        else
            'ZERO GOOD VOXELS!!'
            subAcc(isub,iroi,1:2) = NaN;
            subRandAcc(isub,iroi,1:2,:) = NaN(2,nperms);
            subRand95Acc(isub,iroi,1:2) = NaN;%prctile(bootAccuracy, 95);
        end
        
    end
end

zeroMeanStr = '';
if zeroMean
    zeroMeanStr = '_zeroMean';
end

varinfo=whos('subRandAcc');
saveopt='';
if varinfo.bytes >= 2^31
    saveopt='-v7.3';
end

save([saveFolder 'decodeDrift_direction_' expName '_' crossValidation '_thresh' locStr zeroMeanStr '.mat'],'crossValidation',...
    'r2thresh','locThresh','ROIs','subNames','dataFolder',...
    'subRoiNans','subNumTrials','subNumScans',...
    'subAcc','subRandAcc','subRand95Acc','nperms','locSharedVox','totalVox','locNumVox');
toc