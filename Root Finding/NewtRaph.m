% Ryan Jaipersaud
% ChE 352
% Newton Raphson script

% This function will find the root of an unknown function f. It will return
% the number of iterations it took to find the root. A root will be
% provided once the relative error is below the tolerance. The function requires
% the user to specify a starting guess, the tolerance, f and fprime. This
% uses the NR method to provide a solution.

function [root,NumIter] = NewtRaph(p,TOL,f,fprime)


N = 50;     % max number of iterations
NumIter = 1;  % default value of count

for NumIter = 1:N % this will cause the code to iterate a maximum of N times

     oldp = p;
     p = p - (f(p)/fprime(p)); % this creates the next guess in the NR method
    
    if (abs((oldp-p)/p)< TOL) % this checks to see if the function at the root is within the tolerance
        root = p;
        fprintf('The program iterated %i times\n',NumIter); % prints the number of iterations it took to converge
        fprintf('The root is %4.15f',root); % prints the root
        break; % exits the for loop when the condition is met
   
        
    end
    
end

if ( NumIter >= N ) % Check if maximum number of iterations was reached
    root = 'Error no root found';
    fprintf('%s\n', root);
    fprintf('WARNING: Program reached the maximum number of iterations.\n') % error for when solution didn't converge within alloted number of iterations
end

end