% Ryan Jaipersaud
% ChE 352
% The following code solves a system of
% 20 equations in 25 unknowns.
% Newton's method for solving nonlinear equations 
% will be implemented to solve for y and P

function [x, y, V, P] = flash2(zF,T,L)
F = norm(zF,1); % returns feed

% checks if all Psats are valid and if a feed exists and if L is greater
% than 0 but less than F
if (sum(isnan([Psat1(T),Psat2(T),Psat3(T),Psat4(T)]))== 0 && F ~=0 && zF(1) >= 0 && zF(2) >= 0 && zF(3) >= 0 && zF(4) >= 0 && L<F && L>0)
z = zF/F; % returns z vector
TOL = 10^-11; % sets tolerance
V = F-L; % calculates vapor feed

j =[.25;.25;.25;.25;4]; % [y1,y2,y3,y4,P] initial guess

% start of newton's method
for N = 1:100

% analytical jacobian
J = [ -V-((j(5)*L)/Psat1(T)),0,0,0,(-j(1)*L)/Psat1(T);
      0,-V-((j(5)*L)/Psat2(T)),0,0,(-j(2)*L)/Psat2(T);
      0,0,-V-((j(5)*L)/Psat3(T)),0,(-j(3)*L)/Psat3(T);
      0,0,0,-V-((j(5)*L)/Psat4(T)),(-j(4)*L)/Psat4(T);
      -1,-1,-1,-1,0];

% Matrix of nonlinear equations 
f = [ (z(1)*F)-(j(1)*V)-((j(5)*j(1)*L)/Psat1(T)); (z(2)*F)-(j(2)*V)-((j(5)*j(2)*L)/Psat2(T)); (z(3)*F)-(j(3)*V)-((j(5)*j(3)*L)/Psat3(T));(z(4)*F)-(j(4)*V)-((j(5)*j(4)*L)/Psat4(T));1-j(1)-j(2)-j(3)-j(4)];

u =J\f; % u is my y superscript k as shown in the lecture notes 
j = j-u; % updates j

if(norm(u)<TOL) % breaks out of for loop if tolerance is met
    break;
end

end


% Assigns P and calculate bubble/dew pressures
P = j(5);
Pbubble = (z(1)*Psat1(T))+(z(2)*Psat2(T))+(z(3)*Psat3(T))+(z(4)*Psat4(T));
Pdew = 1/((z(1)/Psat1(T))+(z(2)/Psat2(T))+(z(3)/Psat3(T))+(z(4)/Psat4(T)));


%checks if calculated P is within the envelope
if( (P>Pdew) && (P<Pbubble) )
    % assignment of variables
    y = [j(1);j(2);j(3);j(4)];    
    k =[Psat1(T)/P; Psat2(T)/P; Psat3(T)/P; Psat4(T)/P]; 
    x = [y(1)/k(1); y(2)/k(2); y(3)/k(3); y(4)/k(4)];
    
else
    % assignment when not within envelope
    fprintf('Error flash2: P is not between Pbubble and Pdew\n');
    x = NaN;
    y = NaN;
    V = NaN;
    P = NaN;
end

else % for when at least one input is invalid
    fprintf('Error: At least one input was invalid');
    x = NaN;
    y = NaN;
    V = NaN;
    P = NaN;
end
end

