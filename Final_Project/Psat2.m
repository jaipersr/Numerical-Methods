% Ryan Jaipersaud
% ChE 352
% This functions will return the vapor pressure for vinyl chloride monomer at a 
% specficed T between 170 and 432 Kelvin


function [VaporPressure] = Psat2(T)
    if((T<432)&&(T>170))
         [VaporPressure] = exp(9.001 -(1946/(T-43.15)));
    else
        fprintf('Error Psat2: Not in valid range\n');
        [VaporPressure] = NaN;
    end
    

end