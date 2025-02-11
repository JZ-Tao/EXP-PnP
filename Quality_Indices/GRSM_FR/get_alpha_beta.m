%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Alpha and Beta for QNR and QNR-like protocols
% 
% Interface:
%           [alpha, beta] = get_alpha_beta(tag_qnr,tag_case,tag_sensor)
%
% Inputs:
%           tag_qnr    [string] = Protocol ('QNR','FQNR','HQNR','RQNR')
%           tag_case   [string] = Estimation case ('unitary','sensor-based','overall')
%           sensor     [string] = Sensor name ('GeoEye1','IKONOS','WV3') [required for the sensor-based case]
% 
% Outputs:
%           alpha      [scalar] = Exponential weight for spectral quality
%           beta       [scalar] = Exponential weight for spatial quality
%
% References:
%           [Vivone15]  G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”, 
%                       IEEE Transactions on Geoscience and Remote Sensing, vol. 53, no. 5, pp. 2565–2586, May 2015.
%           [Vivone21]  G. Vivone, M. Dalla Mura, A. Garzelli, R. Restaino, G. Scarpa, M. O. Ulfarsson,L. Alparone, and J. Chanussot, "A new benchmark based on recent advances
%                       in multispectral pansharpening: Revisiting pansharpening with classical and emerging pansharpening methods", IEEE Geoscience and Remote
%                       Sensing Magazine, 9(1):53-81, 2021.          
%           [Arienzo22] A. Arienzo, G. Vivone, A. Garzelli, L. Alparone and J. Chanussot, "Full Resolution Quality Assessment of Pansharpening: Theoretical and hands-on Approaches”, 
%                       IEEE Geoscience and Remote Sensing Magazine, 10(2):2-35, 2022.
%
% Notes: If the specified sensor is not present in the list,
%        by default the overall case is considered.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha,beta] = get_alpha_beta(tag_qnr,tag_case,tag_sensor)

% Check of inputs
qnr_like_list = {'QNR','FQNR','HQNR','RQNR'};
sensors_list  = {'GeoEye1','IKONOS','WV3'};
case_list     = {'unitary','sensor-based','overall'};

if ismember(tag_qnr,qnr_like_list)
else
    error('Protocol not found!! Please try one of the following: QNR, FQNR, HQNR, RQNR.');
end

if ismember(tag_case, case_list)
else
    error('Estimation case not found!! Please try one of the following: unitary, sensor-based, overall.');
end

if strcmp(tag_case,'sensor-based') && nargin == 2
    error('Missing input parameter. Sensor not specified!');        
end

if nargin == 3
    if ismember(tag_sensor,sensors_list) || strcmp(tag_case,'unitary')
    else
        tag_case = 'overall';
    end
else
end

% Retrieval of alpha and beta
if strcmp(tag_case,'unitary')
	alpha = 1;
	beta  = 1;
else 
    switch tag_qnr
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'QNR'
            switch tag_case
                case 'sensor-based'
                    switch tag_sensor                        
                        case 'GeoEye1'
                            alpha = 0.01;
                            beta  = 1.220;	
                        case 'IKONOS'
                            alpha = 0.000;
                            beta  = 1.030;	
                        case 'WV3' 
                            alpha=0.000;
                            beta=1.330;
                    end	
                case 'overall'			
                    alpha = 0.000;
                    beta  = 1.190;
            end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'HQNR'
            switch tag_case
                case 'sensor-based'
                    switch tag_sensor
                        case 'GeoEye1' 
                            alpha = 0.890;
                            beta  = 0.910;	
                        case 'IKONOS' 
                            alpha = 0.860;
                            beta  = 0.770;	
                        case 'WV3' 
                            alpha = 1.420;
                            beta  = 0.940;
                    end
                case 'overall'
                    alpha = 0.980;
                    beta  = 0.890;
            end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'FQNR'
            switch tag_case
                case 'sensor-based'
                    switch tag_sensor
                        case 'GeoEye1'
                            alpha = 1.200;
                            beta  = 1.120;	
                        case 'IKONOS'
                            alpha = 0.790;
                            beta  = 1.240;	
                        case 'WV3' 
                            alpha = 1.490;
                            beta  = 0.980;
                    end	
                case 'overall'			
                    alpha = 1.080;
                    beta  = 1.100;
            end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'RQNR'
            switch tag_case
                case 'sensor-based'
                    switch tag_sensor
                        case 'GeoEye1'
                            alpha = 1.900;
                            beta  = 0.430;	
                        case 'IKONOS'
                            alpha = 1.380;
                            beta  = 0.760;	
                        case 'WV3'
                            alpha = 2.450;
                            beta  = 1.310;
                    end	
                case 'overall'			
                    alpha = 1.810;
                    beta  = 0.900;
            end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end