% Ryan Jaipersaud
% ChE 352
% The following code solves a system of
% 20 equations in 25 unknowns.
% It will first check if the input pressure is within the dew and bubble
% points and then execute flash1. If the previous condition was not met it
% will check to see if the amount of liquid is between 0 and the feed and
% execute flash2. Checks if y1 is between 0 and 1 then executes flash3.
% The default is returning all outputs as NaN.

function [x, y, V, L, P] = RJflash(zF,T,P,L,y1)

% calculates F and z vector
F = norm(zF,1); 
z = zF/F;

% calculates Pbubble and Pdew
Pbubble = (z(1)*Psat1(T))+(z(2)*Psat2(T))+(z(3)*Psat3(T))+(z(4)*Psat4(T));
Pdew = 1/((z(1)/Psat1(T))+(z(2)/Psat2(T))+(z(3)/Psat3(T))+(z(4)/Psat4(T)));

if( (P<Pbubble) &&(P>Pdew))
    [x, y, V, L] = flash1(zF,T,P); % executes flash1 is P is within envelope
elseif((L>0)&&(L<F))% for when the prior condition is not met
   [x, y, V, P]= flash2(zF,T,L); % exectue flash2 if L is between 0 and F
elseif((y1>0)&&(y1<1))% for when neither of the first two conditions are met
   [x, y2, y3, y4, V, L, P]= flash3(zF,T,y1);% executes flash3
   y = [y1;y2;y3;y4];
else % default for when no conditions are met
    fprintf('Error: All inputs invalid');
    x = NaN;
    y = NaN;
    V = NaN;
    L = NaN;
    P = NaN;
end
end