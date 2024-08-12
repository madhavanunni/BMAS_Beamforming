%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit-Receive angle estimation                                      %
%                                                                        %
% Estimates all possible transmit to recieve angles using Law of Cosines %
% Also estimates the alpha value pairs where alpha left and alpha right  %
% is equal.                                                              %
%                                                                        %
% Author: Madhavanunni A N                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alphas = getTxRxAngles(sysPara,bf_points,txElements,vectors)

n = size(sysPara.elem_pos,1);
bfPara.no_emission_types = size(txElements,1);
alphas = single(NaN(bfPara.no_emission_types,n,size(bf_points,1)));
etype = 1; %for zero angle plane wave
for bfp=1:size(bf_points,1)
    tx_idx = txElements(etype,bfp);
    for i=1:n
        vectors.tx_to_rx(etype,i) = abs(sysPara.elem_pos(tx_idx)-sysPara.elem_pos(i));
        if vectors.tx_to_rx(etype,i)==0 %%To avoid complex values of alpha
            nmrtr = 1;
            dnmtr = 1;
        else
            nmrtr = (vectors.rx_distance(i,bfp))^2+(vectors.tx_distance(etype,bfp))^2-(vectors.tx_to_rx(etype,i))^2;
            dnmtr = (2*vectors.rx_distance(i,bfp)*vectors.tx_distance(etype,bfp));
        end
        if nmrtr>dnmtr
            warning('Complex alpha values');
        end
        alphas(etype,i,bfp) = acos(nmrtr/dnmtr);
    end
    if rem(bfp,1000)==0
        clc
        disp(strcat(num2str(bfp),'/',num2str(size(bf_points,1))))
    end
end
