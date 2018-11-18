% Ryan Jaipersaud
% ChE 352
% Bisection script

% This function will find the root of an unknown function f. It will return
% the number of iterations it took to find the root. A root will be
% provided once the relative error is below the tolerance. The function requires
% the user to specify a a lower and upper bound, the tolerance, and f. This
% uses the Bisection method to provide a solution.


function [root,NumIter] = Bisection(a,b,TOL,f)


N = 50;     % max number of iterations
root = a + (b-a)/2;      % sets the value of root to the midpoint of a and b
NumIter = 0;  % counter for iterations


if ( f(a)*f(b)< 0 ) % checks to see if there is a root on the interval
  
    while ( abs((root-a)/root) >= TOL ) % continues to execute code as long as the difference between the bounds is greater than the tolerance

        
        NumIter = NumIter + 1; % increments iterator
        
        if ( f(root)*f(b) < 0 ) % checks where the root is
            a = root;       % makes the root the left bound for the next iteration
        else
            b = root;       % default sets the root as the right bound for the next iteration
        end
        
        root = a + (b-a)/2;  % sets the value of root to the midpoint of a and b
        
        if (NumIter == N)  % checks if max iterations is reached
            fprintf('Error max iterations reached\n'); % error message
            break;       % breaks out of while loop
        end
        
    end
    fprintf('The root is %4.15f.\nThe program iterated %4i times.\n',root, NumIter); % prints root and iterations for convergence
else
    root = 'Error no root found'; % sets root to return a string
    fprintf('%s\n', root); % error message
    
end

end
    