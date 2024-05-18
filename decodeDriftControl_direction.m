% decodeDriftControl_direction.m
%
% associated with the following publication: Steinberg, NJ, Roth, ZN, Movshon, JA, and Merriam, EP (2024).
% Brain representations of motion and position in the double-drift illusion
% 
% DOI:https://doi.org/10.7554/eLife.76803
%
%   usage: decodeDriftControl_direction()
%   by: zvi roth
%   date: 3/05/2024
%   purpose: Decode illusion direction in Expt 3.
% 

close all
clear all
saveFolder = '~/data/';%location of saved timeseries


nperms = 1000; %number of permutations for decoding significance

zeroMean = 1; %0 = no normalization. 1 = subtract mean.

expNum=3;
crossValidation = 'run';%'trial','run'
tic
r2thresh = 50;%r2 pbercentile

nconds=4; %number of conditions
expName = '4conds';

load([saveFolder 'driftTS_' expName '.mat'],'ROIs','subNames','dataFolder',...
    'driftBetas','subRoiNans','subNumTrials','subNumScans');

numSubs = length(subNames);
numRois = length(ROIs);


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
        end
        driftBetas{isub,numRois+xroi} = tempStruct;
    end
end

numRois = length(ROIs);%now including combined ROIs

for isub=1:numSubs
    nTrials = subNumTrials(isub);
    for iroi=1:numRois
        roiname = ROIs{iroi};
        roiNans = subRoiNans{isub,iroi};
        trialReg = driftBetas{isub,iroi};
        if ~isempty(trialReg)
            goodVox = trialReg.r2>prctile(trialReg.r2,50);
            totalVox(isub,iroi) = length(goodVox);
            numVox(isub,iroi) = sum(goodVox);
            ['sub' num2str(isub) ' ' roiname ': using ' num2str(numVox(isub,iroi)) ' voxels out of ' num2str(totalVox(isub,iroi)) ' voxels']
            if sum(goodVox)>0

                for idecode=1:2 % 1==illusion,  2 == just local motion
                    if idecode==1
                        betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
                    else
                        betas = trialReg.ehdr(goodVox,[nTrials*2+1:3*nTrials (3*nTrials)+1:nTrials*4])'; %this is for 3 vs 4
                    end
                    betas = zscore(betas);
                    
                    if zeroMean %subtract each voxel mean
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
                        nRuns = length(group)/2/3;
                        runVec = repmat(1:nRuns, 3, 1);
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
                            idx1 = (iScan*trialsPerCond/2)-2;
                            idx2 = iScan*trialsPerCond/2;
                            randLabels(idx1:idx2,iBoot) = tempList(1:trialsPerCond/2);
                            randLabels(idx1+length(group)/2 : idx2+length(group)/2,iBoot) = tempList((trialsPerCond/2)+1:end);
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
                    
                    subAcc(isub,iroi,idecode) = accuracy;
                    subRandAcc(isub,iroi,idecode,:) = bootAccuracy;
                    subRand95Acc(isub,iroi,idecode) = prctile(bootAccuracy, 95);
                end
            else
                'ZERO GOOD VOXELS!!'
                subAcc(isub,iroi,1:2) = NaN;
                subRandAcc(isub,iroi,1:2,:) = NaN(nperms,1);
                subRand95Acc(isub,iroi,1:2) = NaN;
            end
        else
            'ZERO GOOD VOXELS!!'
            subAcc(isub,iroi,1:2) = NaN;
            subRandAcc(isub,iroi,1:2,:) = NaN(2,nperms);
            subRand95Acc(isub,iroi,1:2) = NaN;
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

%%
save([saveFolder 'decodeDrift_direction_' expName '_' crossValidation zeroMeanStr '.mat'],'crossValidation',...
    'r2thresh','ROIs','subNames','dataFolder',...
    'subRoiNans','subNumTrials','subNumScans',...
    'subAcc','subRandAcc','subRand95Acc','nperms','numVox','totalVox');
toc