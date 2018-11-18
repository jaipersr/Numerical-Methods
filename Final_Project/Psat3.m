% Ryan Jaipersaud
% ChE 352
% This functions will return the vapor pressure for isobutane at a 
% specficed T between 170 and 408 Kelvin
function [VaporPressure] = Psat3(T)
    if((T<408)&&(T>170))
     [VaporPressure] = exp(9.186 -(2128/(T-29.53)));
    else
        fprintf('Error Psat3: Not in valid range\n');
        [VaporPressure] = NaN;
        
    end
end
