%-------------------------------------------------------------------------%
%-- Implements the conventional Delay And Sum (DAS) beamforming
%-- Dependencies:
%-- apodization.m from PICMUS evaluation code
%-------------------------------------------------------------------------%
function beamformedData=DAS(dataCube,x_grid,z_grid,probe_geometry)

N_elements=size(dataCube,3);
x = x_grid(:);
z = z_grid(:);
pixels=size(z_grid,1)*size(z_grid,2);

rx_f_number = 1.75;
rx_aperture = z/rx_f_number;

rx_aperture_distance = abs(x*ones(1,N_elements)-ones(pixels,1)*probe_geometry.');
receive_apodization = apodization(rx_aperture_distance,rx_aperture*ones(1,N_elements),'tukey25');%'hanning'

beamformedData = sum(receive_apodization.*reshape(dataCube,pixels,N_elements),2);
beamformedData(isnan(beamformedData))=0;
beamformedData=reshape(beamformedData,size(x_grid));