% Spencer King and Ryan Peters 
% Aersp 450 
% Proj 1
% 9/19/2020

%% Project statement

% Determine quaternions for given anuglar velocity and time to derive the
% Classical Rodrigues Parameters

%% Declaring variables and constants

t=linspace(0,120,12001); %Every .01 seconds up to 2 mins
w=20;                      %Constant
Wbn = [w*sin(0.01*t); .01+t-t; w*cos(0.01*t)];

%% 1.a

for i = 1:length(t)
    Quaternion = [ 1 Wbn(1,i) Wbn(2,i) Wbn(3,i) ];

end


%% 1.b

BWbn = [ 0 -Wbn(1) -Wbn(2) -Wbn(3) ;
           Wbn(1) 0 Wbn(3) -Wbn(2) ; 
           Wbn(2) -Wbn(3) 0 Wbn(1) ; 
           Wbn(3) Wbn(2) -Wbn(1) 0 ];
    
Phi = expm(-.5*bwBN*T);

%% 1.c

CRP = [Quaternion(2)/Quaternion(1) Quaternion(3)/Quaternion(1) Quaternion(4)/Quaternion(1) ];

 