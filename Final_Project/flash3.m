% Ryan Jaipersaud
% ChE 352
% The following code solves a system of
% 20 equations in 25 unknowns.
% This code initial uses Newton's method for non linear equations to create
% a starting guess for steepest decent for nonlinear equation.
% A simple line search was used.

function [ x, y2, y3, y4, V, L, P] = flash3(zF,T,y1)
F = norm(zF,1); % this finds my feed

% checks if all Psats are valid and if a feed exists and if y1 is greater
% than z but less than 1
if (sum(isnan([Psat1(T),Psat2(T),Psat3(T),Psat4(T)])) == 0 && F~=0 && zF(1) >= 0 && zF(2) >= 0 && zF(3) >= 0 && zF(4) >= 0 && y1 > z(1) && y1<1)
z = zF/F; % this calculates my z vector
TOL = 10^-10; % sets Tolerance
N = 1000; % sets max number of iterations of outer for loop

Pbubble = (z(1)*Psat1(T))+(z(2)*Psat2(T))+(z(3)*Psat3(T))+(z(4)*Psat4(T)); % calculates bubble pressure
Pdew = 1/((z(1)/Psat1(T))+(z(2)/Psat2(T))+(z(3)/Psat3(T))+(z(4)/Psat4(T))); % calculates dew pressure
guess1 = (1-y1)/3; % initial guess for y2,y3,y4
guess = [guess1; guess1; guess1; F/2; F/2;(Pbubble + Pdew)/2]; %[y2,y3,y4,V,L,P] initial guess vector

% This is the sum of the square of functions
G = @ (j) (F*z(1)-y1*j(4)-(j(5)*y1*j(6)/Psat1(T)))^2 ...
         +(F*z(2)-j(1)*j(4)-(j(5)*j(1)*j(6)/Psat2(T)))^2 ...
         +(F*z(3)-j(2)*j(4)-(j(5)*j(2)*j(6)/Psat3(T)))^2 ...
         +(F*z(4)-j(3)*j(4)-(j(5)*j(3)*j(6)/Psat4(T)))^2 ...
         +(1-y1-j(1)-j(2)-j(3))^2 ...
         +(F-j(4)-j(5))^2;
% This is a column vector of the functions     
func = @(j)  [F*z(1)-y1*j(4)-(j(5)*y1*j(6)/Psat1(T));
              F*z(2)-j(1)*j(4)-(j(5)*j(1)*j(6)/Psat2(T));
              F*z(3)-j(2)*j(4)-(j(5)*j(2)*j(6)/Psat3(T));
              F*z(4)-j(3)*j(4)-(j(5)*j(3)*j(6)/Psat4(T));
              1-y1-j(1)-j(2)-j(3);
              F-j(4)-j(5)];
% this is the Jacobian of the function vector
J  = @(j) [0, 0 , 0,                          -y1  , -y1*j(6)/Psat1(T), -j(5)*y1/Psat1(T);
           -j(4)-(j(5)*j(6)/Psat2(T)), 0 , 0, -j(1), -j(1)*j(6)/Psat2(T), -j(5)*j(1)/Psat2(T);
           0, -j(4)-(j(5)*j(6)/Psat3(T)) , 0, -j(2), -j(2)*j(6)/Psat3(T), -j(5)*j(2)/Psat3(T);
           0, 0 , -j(4)-(j(5)*j(6)/Psat4(T)), -j(3), -j(3)*j(6)/Psat4(T), -j(5)*j(3)/Psat4(T);
           -1,-1,-1,0,0,0;
           0,0,0,-1,-1,0];
for iter = 1:3 % Using Newton's method to get a decent starting guess
u =J(guess)\func(guess); % u is my y superscript k as shown in the lecture notes 
guess = guess-u; % updates guess for 3 iterations
end

% this is the start of steepest decent
for i=1:N

grad = 2*(J(guess))'*func(guess); % calculates the gradient
old = G(guess); % calculates G at the current guess

for u = 1:100 % the purpose of this to find an appropriate step that decreases my G
    t = (0.5)^u; % calculates step size
    new = G(guess-t*grad); % calculates G in the new step direction
    if ( new < old) % checks if the G at the new step is less than the previous
        break; % breaks out of inner for loop
    end 
end

guess = guess - t*grad; % updates guess when new < old

if( t*(norm(-grad)/norm(guess)) < TOL) % checks guess
    break; % this will break out of outer for loop if TOL is met
end

end


P = guess(6); % assigns P
if ( P<Pbubble && P>Pdew) % checks if P is withing envelope
   
    % assignment of variables
    y = [y1; guess(1);guess(2);guess(3)];
    y2 = guess(1); 
    y3 = guess(2);
    y4 = guess(3);
    V = guess(4);
    L = guess(5);
    k = [Psat1(T)/P,Psat2(T)/P,Psat3(T)/P,Psat4(T)/P];
    
    
% calculates x
for o =1:4
    x(o,1) = y(o)/k(o);
end
    
else % This is for when P is outside of Pbubble and Pdew
    fprintf('Error flash3: P is not between Pbubble and Pdew\n');
    x = NaN;
    y2 = NaN; 
    y3 = NaN;
    y4 = NaN;
    V = NaN;
    L = NaN;
    P = NaN; 
end
else % for when at least one input is invalid
    fprintf('Error: At least one input was invalid');
    x = NaN;
    y2 = NaN; 
    y3 = NaN;
    y4 = NaN;
    V = NaN;
    L = NaN;
    P = NaN;   
end
end



