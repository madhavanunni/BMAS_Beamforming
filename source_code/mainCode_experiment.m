%-------------------------------------------------------------------------%
%-- Script for DAS and Beam Multiply and Sum (BMAS) beamfomring
%-- "Beam multiply and sum: A tradeoff between delay and sum and filtered 
%-- delay multiply and sum beamforming for ultrafast ultrasound imaging".
%-- Authors: Madhavanunni A N and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified: July 2024
%-------------------------------------------------------------------------%
%-- Dependencies:
%-- PICMUS datasets (carotid_cross_expe_dataset_rf)
%-------------------------------------------------------------------------%
%-- Acknowledgements for the dataset, and apodization function:
%-- 1. H. Liebgott, A. Rodriguez-Molares, F. Cervenansky, J. A. Jensen, 
%-- and O. Bernard, “Plane-Wave Imaging Challenge in Medical Ultrasound,” 
%-- IEEE International Ultrasonics Symposium, IUS, vol. 2016-November, 2016.
%-------------------------------------------------------------------------%
clear
% tic
%% Set all the required paths
tic
addpath('functions');
addpath('lib');
%% Enable/Disable the required beamformers
% inVivo = 1;
% rfDataType = 'PICMUS';
DAS_enable = 0;
pDAS_enable = 0;
BMAS_enable = 1;
%% Load PICMUS dataset
dataPath = 'data\';
fileName = 'carotid_cross_expe_dataset_rf';

path_dataset = strcat(dataPath,fileName,'.hdf5');
rawData = us_dataset();
rawData.read_file(path_dataset);

sysPara.txType = 'plane_wave';
sysPara.txAngles = rawData.angles;
sysPara.maxRxSteerAng = max(sysPara.txAngles(:))*180/pi;
sysPara.rxStartTime = 0;
sysPara.nElements = size(rawData.data,2);         %rawDataset.channels;
sysPara.elem_pos = rawData.probe_geometry;        % Element positions  [m]
sysPara.fs = rawData.sampling_frequency;          % Sampling frequency [Hz]
sysPara.c = rawData.c0;                           % Sound speed [m/s]
sysPara.fc = 5.2e06;                              % center frequency
sysPara.lambda = sysPara.c/sysPara.fc;            % Wavelength [m]
sysPara.nFrames = rawData.frames;
sysPara.PRF = rawData.PRF;
sysPara.no_emission_types = rawData.firings;
sysPara.flow_em_per_frame = 1;                    % Number of flow emissions per frame
sysPara.focal_delay = zeros(size(rawData.angles));


rawDataset = rawData.data;

%% Select the beamformer parameters
txAngIdx_0 = sysPara.txAngles==0;
maxRxSteerAng = sysPara.maxRxSteerAng; %in degrees
%% Define the grid parameters
timeVector(:,1)=((0:size(rawData.data,1)-1)/sysPara.fs);%+sysPara.rxStartTime;
zAxisOrg = timeVector*sysPara.c/2;
% Choose the ROI --- x_axis and z_axis to beamform

clear gridParams
xAxisOrg = sysPara.elem_pos(:,1);
gridParams.x_axis = xAxisOrg;
gridParams.z_axis = zAxisOrg(zAxisOrg>=0.005 & zAxisOrg<=0.05);
gridParams.theta = 0; % Select the grid rotation angle for directional beamforming

%% Load the estimation points grid create the grid (grid for conventional BF)
zl=length(gridParams.z_axis);
xl=length(gridParams.x_axis);

[xGrid,zGrid] = meshgrid(gridParams.x_axis,gridParams.z_axis);
est_points = [xGrid(:), zeros(size(xGrid(:))), zGrid(:)];

% Generate the beamforming points
bf_points = est_points;
clear est_points
%% Focal delay
focal_delay = sysPara.focal_delay;
rfDataset = squeeze(rawDataset(:,:,txAngIdx_0));

probe_geometry = sysPara.elem_pos(:,1);
sysPara.elemPitch = sysPara.elem_pos(2,1)-sysPara.elem_pos(1,1);

%% Conventional DAS beamforming for 75 angles
if DAS_enable == 1
    beamformedData_DAS = single(zeros(size(zGrid,1),size(zGrid,2),numel(sysPara.txAngles)));
    for txAngIdx = 1:numel(sysPara.txAngles)
        thisDataCube=data_cube_extract(squeeze(rawDataset(:,:,txAngIdx)),timeVector,...
            focal_delay(txAngIdx),xGrid,zGrid,probe_geometry,sysPara.c,sysPara.txAngles(txAngIdx));
        beamformedData_DAS(:,:,txAngIdx) = DAS(thisDataCube,xGrid,zGrid,probe_geometry);
        clc;
        disp(strcat('DAS beamforming: ',num2str(txAngIdx),'/',num2str(numel(sysPara.txAngles))));
    end
    beamformedData_DAS(isnan(beamformedData_DAS))=0;
    disp('DAS beamforming completed for all angles')
else
    disp('DAS beamforming skipped')
end

%% BMAS beamforming
if BMAS_enable == 1
    set_BMAS_params_optimal();

    %%% BMAS beamforming for 75 angles
    beamApod_enable = 0;
    beamformedData_BMAS = single(zeros(2*size(zGrid,1),size(zGrid,2),numel(sysPara.txAngles)));
    wb = waitbar(0,'BMAS beamforming');
    for txAngIdx = 1:numel(sysPara.txAngles)
        waitbar(txAngIdx/numel(sysPara.txAngles),wb,sprintf('BMAS Beamforming: %d of %d',txAngIdx,numel(sysPara.txAngles)));
        beamformedData_BMAS(:,:,txAngIdx) = rfBMAS(squeeze(rawDataset(:,:,txAngIdx)),timeVector,...
            focal_delay(txAngIdx),xGrid,zGrid,probe_geometry,sysPara.txAngles(txAngIdx),sysPara,bfParams,vectors,beamApod_enable);
        clc;
        disp(strcat('BMAS beamforming: ',num2str(txAngIdx),'/',num2str(numel(sysPara.txAngles))));
    end
    close(wb);
    beamformedData_BMAS(isnan(beamformedData_BMAS))=0;
    disp('BMAS beamforming completed for all angles')
else
    disp('BMAS beamforming skipped')
end

%%
figure()
envData = abs(hilbert(squeeze(mean(beamformedData_BMAS,3))));
imagesc(gridParams.x_axis*100,gridParams.z_axis*100,20*log10(envData./max(envData(:))));
axis equal tight
xlabel('x [cm]')
xlabel('z [cm]')
colormap("gray")
clim([-55,0])




