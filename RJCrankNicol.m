% Ryan Jaipersaud
% ChE 352
% The following code implements the Crank Nicolson method for solving a
% partial differential equation that varies in time and in space. This
% method will solve the mass transfer equation in x for the concentration
% between the wall and the reactor.

a = (3.4*10^-5); % the alpha value squared
L = .05925; % This length returns a concentration of .8 at the wall after 10 seconds
time = 10; % sets the time I want the code to run
N = 100; % sets N
M = 20; % sets M

k = time/(M-1); % a step in time
h =L/N;% a step in distance
L = (k*(a))/(h^2); % this is lambda

% This creates my A for any N
A = zeros(N,N);
for i = 1:N
    if(i == N)
        A(i,i) = (L/2)+1; % the element in position N,N must change
    else
    A(i,i) = L+1; % the main diagonal elements
    A(i,i+1) = -L/2; % the off diagonal elements
    A(i+1,i) = -L/2; % the off diagonal elements
    end
end

% creates B for any N
B = zeros(N,N);
for i = 1:N
    if(i == N)
        B(i,i) = 1-(L/2);  % the element in position N,N must change
    else
    B(i,i) = 1-L; % the main diagonal elements
    B(i,i+1) = L/2; % the off diagonal elements
    B(i+1,i) = L/2; % the off diagonal elements
    end
end

% this is my d vector, only the first cell will be changed
d = zeros(N,1);
d(1) = L*18;
w = zeros(1,N); % creates first guess for w
 
for iter = 1:M-1 % iterates from time 1 to M-1
    
 w(iter+1,1:N) = A\(B*(w(iter,1:N)') + d);  % This is the iterative step for this method  
end

%w this is M by N
w' % this is N by M
  
  
  
  