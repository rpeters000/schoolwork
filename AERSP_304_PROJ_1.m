%Ryan Peters
%AERSP 304
%Project 1
%Using direct cosine matrices to transfer between coordinate frames and
%plot angle, distance, angular momentum, and specific energy vs time
%%
clc
clear
%% Declaring constants
r0 = 7200;
rdot = 0;
theta = 0;
thetadot = .001084;
mu = 398600;
%% Setting tolerances
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
cond = [r0, rdot, theta, thetadot];
tspan = [0, 7121];
[t,y] = ode45(@fun, tspan, cond, options);
b = pi/2;

%Direct cosine matrix that maps a vector in the N-frame to a vector in the E-frame 
Aen = [sin(b).*cos(y(:,3)), -sin(b).*sin(y(:,3)), -cos(b); sin(y(:,3)), cos(y(:,3)), 0; cos(b).*cos(y(:,3)), -cos(b).*sin(y(:,3)), sin(b)];
%% Generating subplots on one figure
subplot(8,1,1)
plot(t,y(:,1))
title('Distance vs Time')

subplot(8,1,2)
plot(t,y(:,3))
title('Angle vs Time')

h=y(:,3).^2 .* y(:,4);
e = (y(:,2).^2)/2 - mu/y(:,1);


subplot(8,1,3)
plot(t,h)
title('Angular Momentum vs Time')

subplot(8,1,4)
plot(t,e)
title('Specific Energy vs Time')


options = odeset('RelTol', 1e-2, 'AbsTol', 1e-2);
cond = [r0, rdot, theta, thetadot];
tspan = [0, 7121];
[t,y] = ode45(@fun, tspan, cond, options);

subplot(8,1,5)
plot(t,y(:,1))
title('Distance2 vs Time')

subplot(8,1,6)
plot(t,y(:,3))
title('Angle2 vs Time')

h=y(:,3).^2 .* y(:,4);
e = (y(:,2).^2)/2 - mu/y(:,1);


subplot(8,1,7)
plot(t,h)
title('Angular Momentum2 vs Time')

subplot(8,1,8)
plot(t,e)
title('Specific Energy2 vs Time')


%%The graphs look the same from the previous even though there is a
%%different tolerance. Likely, the one with the higher tolerance is more
%%accurate than the lower tolerance. 



function dx=fun(t,x)
mu=398600;
dx(1)=x(2);
dx(3)=x(4);
dx(2)=x(1)*x(4)^2 - mu/(x(1)^2);
dx(4)=-(2*x(2)*x(4))/x(1);
dx=dx(:);
end



