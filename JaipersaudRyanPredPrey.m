%Ryan Jaipersaud
%ChE 352

% The following code implements RK4 for a system of differential equations
% The code also uses ode45 to solve the same system
% Questions 2&3
clear all;
clc;
format long;



a = 0; % starting time in units of seconds
b = 4; % ending time in units of seconds
N = 1800; % number of iterations N = 1800 for 1 sec step size, N = 180 for 10 sec step size, N = 30 for 1 minutes step size
h = (b-a)/N; % step size


t = zeros(1,N);     % creates a 1 by N array for t
w1 = zeros(1,N);    % creates a 1 by N array for prey
w2 = zeros(1,N);    % creates a 1 by N array for predator

%initial conditions
w1(1) = 1000;          
w2(1) = 500; 
t(1) = 0;  

% constants
k1 = 3;
k2 = .002;
k3 = .0006;
k4 = .5;


for i = 1:(N-1)
    
  % The second subscript for g indicates the species, the first indicates
  % the number of function evaluations
  g11 = h*(k1*w1(i)-k2*w1(i)*w2(i)); % first function evaluation for RK for prey
  g12 = h*(k3*w1(i)*w2(i)-k4*w2(i)); % first function evaluation for RK for predator
  
  g21 = h*(k1*(w1(i)+(g11/2))-k2*(w1(i)+(g11/2))*(w2(i)+(g12/2)));
  g22 = h*(k3*(w1(i)+g11/2)*(w2(i)+g12/2)-k4*(w2(i)+g12/2));
  
  g31 = h*(k1*(w1(i)+g21/2)-k2*(w1(i)+g21/2)*(w2(i)+g22/2));
  g32 = h*(k3*(w1(i)+g21/2)*(w2(i)+g22/2)-k4*(w2(i)+g22/2));
  
  g41 = h*(k1*(w1(i)+g31)-k2*(w1(i)+g31)*(w2(i)+g32));
  g42 = h*(k3*(w1(i)+g31)*(w2(i)+g32)-k4*(w2(i)+g32));
  
 
  w1(i+1) = w1(i) + (1/6)*(g11+(2*g21)+(2*g31)+g41); % iterative step for RK for prey
  w2(i+1) = w2(i) + (1/6)*(g12+(2*g22)+(2*g32)+g42); % iterative step for RK for predator
  
  t(i+1) = t(i) + h; % increments time 
end

%This will create a plot of w1(prey)& w2(pred) against t
figure(1)
plot(t,w1,'k-','LineWidth',2)
hold on
plot(t,w2,'g-','LineWidth',2)
title('Population Growth Over Time RK4')
xlabel('t')
ylabel('population')
legend('w1','w2','Location', 'Northeast');



% function handle for ode45 to solve
f=@(t,x) [k1*x(1)-k2*x(1)*x(2);k3*x(1)*x(2)-k4*x(2)];
[TOUT,YOUT] = ode45(f,[a b],[1000 500]);
%[TOUT,YOUT] = ode45(f,time interval,initial conditions);

%This will create a plot of ODE45 solution
figure(2)
plot(TOUT,YOUT(:,1),'r-','LineWidth',2)
hold on
plot(TOUT,YOUT(:,2),'k-','LineWidth',2)
title('Population Growth Over Time ODE45')
xlabel('t')
ylabel('population')
legend('x1','x2','Location', 'Northeast');




    