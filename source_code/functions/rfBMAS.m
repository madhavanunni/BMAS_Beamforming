%-------------------------------------------------------------------------%
%-- Script Implements the conventional Beam Multiply And Sum (BMAS) beamforming 
%-- utilizing DAS algebra to reduce computation time
%-- "Beam multiply and sum: A tradeoff between delay and sum and filtered 
%-- delay multiply and sum beamforming for ultrafast ultrasound imaging".
%-- Authors: Madhavanunni A N and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 05 - Aug - 2024
%-------------------------------------------------------------------------%
function beamformedDataBMAS=rfBMAS(rawData,timeVector,focal_delay,x_grid,z_grid,probe_geometry,txAngle,sysPara,bfParams,vectors,beamApod_enable)

N_elements=size(rawData,2);
x = x_grid(:);
z = z_grid(:);
pixels=size(z_grid,1)*size(z_grid,2);
%% Distance and Time of Flight calculation
% Total time of flight calculation
tx_d = repmat((z*cos(txAngle)+x*sin(txAngle)),[1,size(sysPara.elem_pos,1)]).';
vectors.tof_total = ((vectors.rx_distance+tx_d)/sysPara.c)+focal_delay;
clear vec_elem vec_bf tx_d
%% Delay compensation  (Time to Space Mapping)
for nrx=1:sysPara.nElements
    tof_samples(nrx,:) = interp1(timeVector,rawData(:,nrx),vectors.tof_total(nrx,:),'spline',0);
end

%% Beam Apodization
if beamApod_enable==1
    rx_f_number = 1.75;%1.2;1.7
    rx_aperture = vectors.rx_distance/rx_f_number;
    elem_idx = 1:sysPara.nElements;
    elemIdx = single(ones(pixels,1)*elem_idx).';
    for nrx=1:sysPara.nElements
        rx_aperture_distance = abs(sysPara.elem_pos(nrx)-sysPara.elem_pos(:,1)).*ones(N_elements,pixels);
        clc;
        disp(strcat('Beam Apodization: ', num2str(nrx)))
        beamApod = getElemApod(rx_aperture_distance,rx_aperture,elemIdx,'hanning');%'tukey25'
        beamApod(isnan(beamApod))=0;
        beamApod_samples(nrx,:) = squeeze(sum((tof_samples.*beamApod),1)).';
    end
else
    beamApod_samples = tof_samples;
end
clear tof_samples
%% receive apodization
%%-- dynamically expanding receive aperture with a selected apodization

rx_f_number = 1.75;%1.75;
rx_aperture = (vectors.rx_distance/rx_f_number).';%


num_rx_apert = size(bfParams.rxCenterPos,1);
for nrx = 1:num_rx_apert
    rxCenterPos = repmat(bfParams.rxCenterPos(nrx,:).',[1,N_elements]);
    rx_aperture_distance = abs(rxCenterPos-repmat(probe_geometry.',pixels,1));
    receive_apodization = apodization(rx_aperture_distance,rx_aperture,'hanning');%'tukey25'
    receive_apodization(isnan(receive_apodization))=0;
%%--Apply aperture mask to cutoff the weights out od the subaperture spread
%%-- (not mandatory though). This will reduce the required computations
    apertureMask = squeeze(bfParams.apodMask(nrx,:,:));
    receive_apodization = receive_apodization.*apertureMask;

    beamformedData(:,nrx) = sum(receive_apodization.*beamApod_samples.',2);
    clc
    disp(strcat('DAS_aperture:',num2str(nrx),' / ',num2str(num_rx_apert)))
end

clear startElem endElem
%%
apertGain = (bfParams.subApertMask./(repmat(bfParams.nSubAperts,num_rx_apert,1))).';
beamformedData = beamformedData.*apertGain;

%% Upsampling by 2
xl = size(z_grid,2);
zl = size(z_grid,1);
thisBFdata = double(reshape(beamformedData,zl,xl,num_rx_apert));
upSampledData = zeros(num_rx_apert,2*zl,xl);
Fs=2*sysPara.fs;
for nrx=1:num_rx_apert
    upSampledData(nrx,:,:)=resample(squeeze(thisBFdata(:,:,nrx)),2,1);
end
z_axis = linspace(min(z_grid(:,1)),max(z_grid(:,1)),2*zl);
x_axis = x_grid(1,:);

bfSamples = reshape(upSampledData,num_rx_apert,[]);
beamformedDataDAS = zeros(1,2*pixels);
beamformedDataFMAS = zeros(1,2*pixels);
for rxIdx = num_rx_apert:-1:2
    xi = bfSamples(rxIdx,:);
    xi_sqrt=sign(xi).*sqrt(abs(xi));
    beamformedDataDAS = beamformedDataDAS+xi_sqrt;
    
    xj = bfSamples(rxIdx-1,:);
    xj_sqrt=sign(xj).*sqrt(abs(xj));
    beamformedDataFMAS = beamformedDataFMAS+(xj_sqrt.*beamformedDataDAS);
    clc;
    disp(strcat('step: ',num2str(rxIdx)));
end

%% custom filter
% f0new=2*sysPara.fc;
% f0Norm=2*f0new/Fs;
% % h=firpm(199,[0 f0Norm-0.2 f0Norm-0.1 f0Norm+0.1 f0Norm+0.2 1],[0 0 1 1 0 0]);
% h=firpm(199,[0 f0Norm-0.3 f0Norm-0.1 f0Norm+0.1 f0Norm+0.3 1],[0 0 1 1 0 0]);
% % % bfDataDAS = reshape(beamformedDataDAS,xl,zl*2);
% bfDataFMAS = reshape(beamformedDataFMAS,zl*2,xl);
% for xIdx = 1:size(x_grid,2)
% beamformedDataBMAS(:,xIdx)=filtfilt(h,1,(bfDataFMAS(:,xIdx)));
% end
% beamformedDataBMAS = single(beamformedDataBMAS);
% % % beamformedDataBMAS = bfDataFMAS;%reshape(bfDataFMAS,size(z_grid));
% % % figure();imagesc(20*log10(abs(hilbert(beamformedDataBMAS))))

%% BPF 5-20 MHz
beamformedDataFMAS = reshape(beamformedDataFMAS,2*pixels,1);
fStart = sysPara.fc+0.5e6;%12e6;%7e6;
fTransition = 1e6;
fStop = 3*sysPara.fc-1e6;%18e6;%17e6;%19e6;
f1 = [fStart fStart+fTransition fStop fStop+fTransition];
beamformedDataFMAS = band_pass(beamformedDataFMAS,Fs,f1);
beamformedDataBMAS=reshape(beamformedDataFMAS,numel(z_axis),numel(x_axis));
clc
