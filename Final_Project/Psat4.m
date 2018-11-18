% Ryan Jaipersaud
% ChE 352
% This functions will return the vapor pressure for 1-butene at a 
% specficed T between 170 and 418.6 Kelvin
function [VaporPressure] = Psat4(T)
if((T<418.6)&&(T>170))
[VaporPressure] = exp(9.331 -(2190/(T-32.18)));
else
    fprintf('Error Psat4: Not in valid range\n');
        [VaporPressure] = NaN;
end

end