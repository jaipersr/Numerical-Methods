% Ryan Jaipersaud
% ChE 352
% The following code uses the rachford rice equation to solve a system of
% 20 equations in 25 unknowns.
% Root finding will be implemented to solve for psi which is defined as
% the vapor fraction. 

function [x, y, V, L] = flash1(zF, T, P)
F = norm(zF,1);% calculates Feed
if (sum(isnan([Psat1(T),Psat2(T),Psat3(T),Psat4(T)]))== 0 && F~=0 && zF(1) >= 0 && zF(2) >= 0 && zF(3) >= 0 && zF(4) >= 0)
% This calculates my z vector. 
z = zF/F;

% Calculate bubble point and dew point pressure
Pbubble = (z(1)*Psat1(T))+(z(2)*Psat2(T))+(z(3)*Psat3(T))+(z(4)*Psat4(T));
Pdew = 1/((z(1)/Psat1(T))+(z(2)/Psat2(T))+(z(3)/Psat3(T))+(z(4)/Psat4(T)));


%Checks if the input pressure is between the dew and bubble point pressures
if( (P>Pdew) && (P<Pbubble) )
    % Calculate k value for each species
    k = [Psat1(T)/P,Psat2(T)/P,Psat3(T)/P,Psat4(T)/P];
    
    % This is the Rachford Rice equation 
    f= @(psi)( ((z(1)*(1-k(1)))/(1+(psi*(k(1)-1)))) + ((z(2)*(1-k(2)))/(1+(psi*(k(2)-1)))) + ((z(3)*(1-k(3)))/(1+(psi*(k(3)-1)))) + ((z(4)*(1-k(4)))/(1+(psi*(k(4)-1)))) );
    psi = Bisectionnew(0,1,10^-11,f);% This implements root finding to find psi, the subroutine is at the bottom
    
    % Calculate V and L
    V = F*psi;
    L = F - V;

    % Updates the y and x vectors
    for i = 1:4
         y(i,1) = (z(i)*k(i))/(1+(psi*(k(i)-1)));
         x(i,1) = z(i)/(1+(psi*(k(i)-1)));
    end
       
% For when the pressure is not between the bubble and dew point
else
    fprintf('Error flash1: P is not between Pbubble and Pdew\n');
    x = NaN;
    y = NaN;
    V = NaN;
    L = NaN;
end
else % for when at least one input was invalid
    fprintf('Error: At least one input was invalid');
    x = NaN;
    y = NaN;
    V = NaN;
    L = NaN;
end
end

% This function will find the root of a function f.  A root will be
% provided once the relative error is below the tolerance. The function requires
% the user to specify a a lower and upper bound, the tolerance, and f. This
% uses the Bisection method to provide a solution.

function root = Bisectionnew(a,b,TOL,f)


N = 100;  % max iterations


root = []; % default root value

   if ( f(a)*f(b) < 0 )% First check if a root exists
    for i = 1:N       
    
        root = a + (b - a)/2;  % sets the value of root to the midpoint of a and b
               
        if ( abs((root-a)/root) <= TOL )  % relative error stop condition
            break; % Breaks out of for loop if tolerance is met
        elseif ( f(b)*f(root) < 0 )  % Checks if root is btwn b and root
            a = root;% sets root as the left bound for the next iteration   
        else
            b = root;% Set root as the right bound for the next iteration
        end
    
    end

    else% For when there is no root
    disp('Error: no root in interval') 
    end
end
