%% AERSP 458 Project 1
%Ryan Peters and Dominico Vano
%3/4/2020
close all;
clear all;
clc;

%% Declaring constants
mu=4*pi^2;
t0=0;


%% Case 1
r0=[-2.4258, 0.8420, -2.2982]; 
v0=[-1.9857, -0.4478, 0.4309];
v0mag=(norm(v0));
r0mag=(norm(r0));
t1=2.1;
a=mu*r0mag/(2*mu-r0mag*v0mag^2);
sigma0=(dot(r0,v0)/sqrt(mu));
alpha=1/a;

%Applying functions
[X1]=newton_raphson(alpha,mu,t0,t1,sigma0,r0);
[U0_1,U1_1,U2_1]=top_down(alpha,X1);
[r1,v1]=lagrange_sol(r0,v0,U1_1,U2_1);

%Output r1 and v1
fprintf('Case 1\n r = [%.8f, ', r1(1,1));
fprintf(' %.8f, ', r1(1,2));
fprintf(' %.8f] LU\n', r1(1,3));
fprintf(' v = [%.8f, ', v1(1,1));
fprintf(' %.8f, ', v1(1,2));
fprintf(' %.8f] LU/TU\n', v1(1,3));



%% Case 2
r0=[1.3546, 0.4957, -0.8378];
v0=[-1.9497, -2.2778, -9.3040];
v0mag=(norm(v0));
r0mag=(norm(r0));
t1=0.5;
a=mu*r0mag/(2*mu-r0mag*v0mag^2);
sigma0=(dot(r0,v0)/sqrt(mu));
alpha=1/a;

%Applying functions
[X2]=newton_raphson(alpha,mu,t0,t1,sigma0,r0);
[U0_2,U1_2,U2_2]=top_down(alpha,X2);
[r2,v2]=lagrange_sol(r0,v0,U1_2,U2_2);

%Output r2 and v2
fprintf('Case 2\n r = [%.8f, ', r2(1,1));
fprintf(' %.8f, ', r2(1,2));
fprintf(' %.8f] LU\n', r2(1,3));
fprintf(' v = [%.8f, ', v2(1,1));
fprintf(' %.8f, ', v2(1,2));
fprintf(' %.8f] LU/TU\n', v2(1,3));



%% Functions
%Newton Raphson method to find chi(X)
%Universal kepler equation
%Calling top_down function
function [X]=newton_raphson(a,mu,t0,t1,sigma0,r0)

%Initial guess for X
X=pi;
%Set first tolerance check to value greater than tolerance but close enough
%to increase code performance
check=0.0001;

    while check>10^(-6)
        
    [U0,U1,U2]=top_down(a,X);
    
    f=a*sqrt(mu)*(t1-t0)+(1-a*norm(r0))*U1-a*sigma0*U2-X; %Universal time equation
    
    f_prime=(1-a*norm(r0))*U0-a*sigma0*U1-1; %Derivative of universal time equation
    
    X_n=X-(f/f_prime);
    
    check=abs(X_n-X);
    
    X=X_n;
    end
end

% Function to solve for Lagrange coefficient solutions and then r and v vectors
function [r,v]=lagrange_sol(r0,v0,U1,U2)
%Declaring values used in operations
r0mag=(norm(r0));

mu=4*pi^2;

sigma0=(dot(r0,v0)/sqrt(mu));

%Finding F and G coefficients to calculate r
F=1-(U2/r0mag);

G=(r0mag*U1)/sqrt(mu)+(sigma0*U2)/sqrt(mu);

r=F*r0+G*v0;

%Finding Ft and Gt coefficients to calculate v
rmag=(norm(r));

Ft=(-sqrt(mu)*U1)/(rmag*r0mag);
 
Gt=1-(U2/rmag);

v=Ft*r0+Gt*v0;
end

% Top down method to solve for continued fractions in universal functions U0, U1, and U2
function [U0,U1,U2] = top_down(alpha,X)
%Setting initial values
a0=0.5*X;
b0=1;
s0=a0/b0;
u0=a0/b0;
d0=1;

%Setting repeating values
an=alpha*(.5*X)^2;
bn=3;

%Incorporating inital values
bn_1=b0;
dn_1=d0;
un_1=u0;
sn_1=s0;
un=1;
n=1;

    %Creating while loop
    while (abs(un)>10^(-6))
        
        dn=1/(1-(an/bn_1/bn)*dn_1);
        un=un_1*(dn-1);
        sn=sn_1+un;
        
        %Defining loop characteristics
        n=n+1;
        bn=bn+2;
        bn_1=bn_1+2;
        un_1=un;
        sn_1=sn;
        dn_1=dn;
    end
    
u=sn;

%Outputting universal function values    
U0=(1-(alpha*u^2))/(1+(alpha*u^2));
U1=(2*u)/(1+(alpha*u^2));
U2=(2*u^2)/(1+(alpha*u^2));

end





