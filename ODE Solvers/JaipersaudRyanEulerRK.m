clear all;
clc;
format long;

% The following code will implement Euler's method and RK4 to solve an IVP
% Plots with varying step sizes will be produced
% ODE45 is also used to solve the IVP

a = 0; % starting time in units of seconds
b = 1800; % ending time in units of seconds
N = 30; % number of iterations N = 1800 for 1 sec step size, N = 180 for 10 sec step size, N = 30 for 1 minutes step size
h = (b-a)/N; % step size
c = 3.06;   %initial concentration of A and B

Euler = zeros(1,N); % createa a 1 by N array for Euler
t = zeros(1,N);     % creates a 1 by N array for t
RK = zeros(1,N);    % creates a 1 by N array for  RK

Euler(1) = c;       % sets initial concentration value
RK(1) = c;          % sets initial concentration value
t(1) = a;           % sets initial time value

z = pi/120;

for i = 1:(N-1) 
    
  Euler(i+1) = Euler(i) + h*(-4000*(exp(-5000/(60*cos(z*t(i))+320)))*Euler(i)); % Iteration for Euler's method
    
  k1 = h*(-4000*(exp(-5000/(60*cos(z*t(i))+320)))*RK(i)); % first function evaluation for RK
  k2 = h*(-4000*(exp(-5000/(60*cos(z*(t(i)+(h/2)))+320)))*(RK(i)+(k1/2))); % second function evaluation for RK
  k3 = h*(-4000*(exp(-5000/(60*cos(z*(t(i)+(h/2)))+320)))*(RK(i)+(k2/2)));% third function evaluation for RK
  k4 = h*(-4000*(exp(-5000/(60*cos(z*(t(i)+h))+320)))*(RK(i)+k3));% fourth function evaluation for RK
  RK(i+1) = RK(i) + (1/6)*(k1+(2*k2)+(2*k3)+k4); % iterative step for RK
  
  t(i+1) = t(i) + h; % increments time 
end

% This will print out the concentration of B at the end of 30 minutes
fprintf('The amount of CB using a step size of %3.3i seconds according the Euler and RK4 are %4.6i, %4.6i, repsectively.\n',h, 3.06 + (3.06 - Euler(end)),3.06 + (3.06 - RK(end)));

% This will create a figure that will contain plots of Eulers method and
% RK4 solutions
figure(1)
plot(t,Euler,'g-','LineWidth',2)
hold on
plot(t,RK,'b-','LineWidth',2)
title('Concentration of A as a function of time using Euler and RK4');
ylabel('CA [mol/m^3]');
xlabel('Time [s]');
legend('Euler','RK4','Location', 'Northeast');


% This function will be plugged into ODE45 ( x = t, u = CA)
functhand = @(x, u) (-4000.*(exp(-5000/(60.*cos(z.*x)+320)))).*u;
[TOUT,YOUT] = ode45(functhand,[a b],c); % TOUT will store time values, and YOUT will store concentration of a values

%This will create a plot of ODE45 solution
fprintf('The amount of CB according to ODE45 is %4.6f.\n', 3.06 +(3.06 - YOUT(end)));
figure(2)
plot(TOUT,YOUT,'k-','LineWidth',2)
title('Concentration of A as a function of time using using ODE45');
ylabel('CA [mol/m^3]');
xlabel('Time [s]');

%[TOUT YOUT]
    
    
