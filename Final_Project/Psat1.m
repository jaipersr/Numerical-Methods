% Ryan Jaipersaud
% ChE 352
% This functions will return the vapor pressure for ethane at a 
% specficed T between 170 and 305.4 Kelvin

function [VaporPressure] = Psat1(T)
    if((T<305.4)&&(T>=170))
         [VaporPressure] = exp(9.568-(1706/(T-6.064)));
    else
        fprintf('Error Psat1: Not in valid range\n');
        [VaporPressure] = NaN;
    end
    

end