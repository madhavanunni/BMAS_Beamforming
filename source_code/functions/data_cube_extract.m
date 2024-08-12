%-------------------------------------------------------------------------%
%-- Script for conventional geometry based delay compensation.
%-- Authors: Madhavanunni A N and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 05 - Aug - 2024
%-------------------------------------------------------------------------%

function dataCube=data_cube_extract(rawData,timeVector,focal_delay,x_grid,z_grid,probe_geometry,c,txAngle)

N_elements=length(probe_geometry);
x = x_grid(:);
z = z_grid(:);
% c=1540;

dataCube=zeros(size(x_grid,1)*size(x_grid,2),N_elements);
%-- transmit delay
% transmit_delay=z;
transmit_delay = z*cos(txAngle)+x*sin(txAngle);
for nrx=1:N_elements  
    %-- receive delay
    receive_delay = sqrt((probe_geometry(nrx)-x).^2+z.^2);
    %-- total delay
    delay = ((transmit_delay+receive_delay)/c)+focal_delay;  
    dataCube(:,nrx) = interp1(timeVector,rawData(:,nrx),delay,'spline',0);
end

dataCube=(reshape(dataCube,size(x_grid,1),size(x_grid,2),N_elements));