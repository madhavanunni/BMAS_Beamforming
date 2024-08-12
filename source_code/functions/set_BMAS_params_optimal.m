%-------------------------------------------------------------------------%
%-- Script for configuring the parameters for BMAS beamfomring
%-- "Beam multiply and sum: A tradeoff between delay and sum and filtered 
%-- delay multiply and sum beamforming for ultrafast ultrasound imaging".
%-- Authors: Madhavanunni A N and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 05 - Aug - 2024
%-------------------------------------------------------------------------%

vec_elem = reshape(repmat(sysPara.elem_pos,[size(bf_points,1) 1]),[size(sysPara.elem_pos,1) size(bf_points,1) 3]);
vec_bf = permute(reshape(repmat(bf_points,[size(sysPara.elem_pos,1) 1]),[size(bf_points,1) size(sysPara.elem_pos,1)  3]),[2 1 3]);
vectors.elem_to_bf = vec_elem-vec_bf;
vectors.rx_distance = sqrt(sum((vectors.elem_to_bf).^2,3));
[vectors.tx_distance, txElements] = min(vectors.rx_distance,[],1);

bfParams.txRxAngles = getTxRxAngles(sysPara,bf_points,txElements,vectors);
bfParams.txRxAngles = squeeze(bfParams.txRxAngles)*180/pi;
bfParams.txRxAngles(isnan(bfParams.txRxAngles)) = 90;
%%
bfParams.activeApertureOrg = (bfParams.txRxAngles<=maxRxSteerAng);
bfParams.activeApertLen = single(sum(bfParams.activeApertureOrg,1));
bfParams.activeApertLen(bfParams.activeApertLen<=4) = 4;
bfParams.activeApertLen(bfParams.activeApertLen>64 & bfParams.activeApertLen<=96) = 96;
bfParams.activeApertLen(bfParams.activeApertLen>96 & bfParams.activeApertLen<=128) = 128;

for pixelInd = 1:xl*zl
    startInd(1,pixelInd) = single(find(bfParams.activeApertureOrg(:,pixelInd),1,'first'));
    endInd(1,pixelInd) = single(find(bfParams.activeApertureOrg(:,pixelInd),1,'last'));
    startPos = sysPara.elem_pos(startInd(1,pixelInd),1);
    endPos = sysPara.elem_pos(endInd(1,pixelInd),1);
    if startPos==endPos
        activeApertCentre = startPos;
        activeApertStart = activeApertCentre;
        activeApertEnd = activeApertCentre;
    else
        activeApertCentre = startPos+((endPos-startPos)/2);
        activeApertStart = activeApertCentre-(sysPara.elemPitch*bfParams.activeApertLen(:,pixelInd)/2);
        activeApertEnd = activeApertCentre+(sysPara.elemPitch*bfParams.activeApertLen(:,pixelInd)/2);
    end
    bfParams.activeAperture(:,pixelInd) = (sysPara.elem_pos(:,1)>=activeApertStart & sysPara.elem_pos(:,1)<=activeApertEnd);
    
end

bfParams.rxSubApertLen(bfParams.activeApertLen<=8) = 8;%activeApertLens(1);
bfParams.rxSubApertLen(bfParams.activeApertLen>8 & bfParams.activeApertLen<=32) = 8;%8
bfParams.rxSubApertLen(bfParams.activeApertLen>32 & bfParams.activeApertLen<=64) = 16;
bfParams.rxSubApertLen(bfParams.activeApertLen>64 & bfParams.activeApertLen<=96) = 24;
bfParams.rxSubApertLen(bfParams.activeApertLen>96 & bfParams.activeApertLen<=128) = 32;%32


bfParams.rxSubApertStride(bfParams.rxSubApertLen<=4) = 1;%activeApertLens(1);
bfParams.rxSubApertStride(bfParams.rxSubApertLen>4 & bfParams.rxSubApertLen<=8) = 2;
bfParams.rxSubApertStride(bfParams.rxSubApertLen>8 & bfParams.rxSubApertLen<=16) = 4;
bfParams.rxSubApertStride(bfParams.rxSubApertLen>16 & bfParams.rxSubApertLen<=24) = 6;%6;
bfParams.rxSubApertStride(bfParams.rxSubApertLen>24 & bfParams.rxSubApertLen<=32) = 8;%8;%8
bfParams.rxSubApertStride(bfParams.rxSubApertStride==0)=1;
%% %
bfParams.rxCenterPos = single(nan(128,xl*zl));
for pixelInd = 1:xl*zl
    startInd(1,pixelInd) = single(find(bfParams.activeAperture(:,pixelInd),1,'first'));
    endInd(1,pixelInd) = single(find(bfParams.activeAperture(:,pixelInd),1,'last'));
    startPos = sysPara.elem_pos(startInd(1,pixelInd),1);
    endPos = sysPara.elem_pos(endInd(1,pixelInd),1);
    strideLen = sysPara.elemPitch*(bfParams.rxSubApertStride(pixelInd));
    if startPos==endPos
        thisRxCenterPos = startPos;
    else
        thisRxCenterPos = startPos:strideLen:endPos;
    end
    bfParams.rxCenterPos(1:size(thisRxCenterPos,2),pixelInd) = thisRxCenterPos;
end
bfParams.rxCenterPos(all(isnan(bfParams.rxCenterPos),2),:)=[];
bfParams.subApertMask = ~isnan(bfParams.rxCenterPos);
bfParams.nSubAperts = single(sum(bfParams.subApertMask,1));

%%--Fill in dummy values for empty rx centers for computational simplicity
for ii = 1:numel(bfParams.rxCenterPos)
    if isnan(bfParams.rxCenterPos(ii))
        bfParams.rxCenterPos(ii)=bfParams.rxCenterPos(ii-1);
    end
end
%%%
bfParams.subApertStart = bfParams.rxCenterPos-(repmat(sysPara.elemPitch*bfParams.rxSubApertLen,size(bfParams.rxCenterPos,1),1)/2);
bfParams.subApertEnd = bfParams.rxCenterPos+(repmat(sysPara.elemPitch*bfParams.rxSubApertLen,size(bfParams.rxCenterPos,1),1)/2);

%%%
clear elemPos4apodMask apodMask
for nrx = 1:sysPara.nElements 
    elemPos4apodMask = sysPara.elem_pos(nrx,1);
    apodMask(:,:,nrx) = (elemPos4apodMask>=bfParams.subApertStart) & (elemPos4apodMask<=bfParams.subApertEnd);
    
    clc;
    disp(strcat(num2str(nrx),'/',num2str(sysPara.nElements)))
    
end
bfParams.apodMask = apodMask;
clear apodMask
%% %
complexityMetrics.nAdders = (bfParams.nSubAperts.*(bfParams.rxSubApertLen-1))+(bfParams.nSubAperts-1);
complexityMetrics.maxAdderCount = max(complexityMetrics.nAdders);
complexityMetrics.meanAdderCount = mean(complexityMetrics.nAdders);
complexityMetrics.totalAdders = sum(complexityMetrics.nAdders);
complexityMetrics.nMultipliers = bfParams.nSubAperts-1;
complexityMetrics.maxMultiplierCount = max(complexityMetrics.nMultipliers);
complexityMetrics.meanMultiplierCount = mean(complexityMetrics.nMultipliers);
complexityMetrics.totalMultipliers = sum(complexityMetrics.nMultipliers);