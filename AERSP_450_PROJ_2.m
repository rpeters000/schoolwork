% Ryan Peters and Spencer King 
% Aersp 450 
% Project 2

% Through this project, we will consider the spacecraft attitude estimation
% problem while making use of the Davenport’s q method and the 
% Optimal Linear Attitude Estimator (OLAE).

clear all
close all

t=[0:5:3600];
Omega=1.84*(10^(-6));                      %Constant
w1 = Omega*sin(Omega*t);             %Angular velocity array
w2 = Omega*cos(Omega*t);
w3 = Omega+t-t;
wBN = [w1; w2; w3];             %Angular velocity array of sat

r1=[-0.1517;-0.96690;.2050];     %stars in newtonian frame
r2=[-0.83930;.44940;.3044];
r3=[-0.0886;-0.58560;.8000];
r4=[0.8814;-0.03030;.5202];
R = [r1 r2 r3 r4];

%% Calculating CBN,Rodriques, and Quaterions
%CRP2DCM lines 21 to 46
b0 = [ 1 0 0 0];            %Initial conditions of quaternions
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t,b] = ode45(@odefunction,t,b0,options);   %Integrating the Quaternions

Q=b';

for i = 1:length(t)     %Loop iterates through each time step
    CRP(:,i) = [Q(2,i)/Q(1,i) Q(3,i)/Q(1,i) Q(4,i)/Q(1,i) ];  %Used values form b cause they are same as a
    CRPf(i,:) = CRP(:,i)';
end

for i = 1:length(t) %Finding the DCM for each time t 
    
    CBN(1,1,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(1+CRP(1,i)*CRP(1,i)-CRP(2,i)*CRP(2,i)-CRP(3,i)*CRP(3,i));
    CBN(2,2,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(1-CRP(1,i)*CRP(1,i)+CRP(2,i)*CRP(2,i)-CRP(3,i)*CRP(3,i));
    CBN(3,3,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(1-CRP(1,i)*CRP(1,i)-CRP(2,i)*CRP(2,i)+CRP(3,i)*CRP(3,i));
    CBN(1,2,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(2*(CRP(1,i)*CRP(2,i)+CRP(3,i)));
    CBN(1,3,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(2*(CRP(1,i)*CRP(3,i)-CRP(2,i)));
    CBN(2,1,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(2*(CRP(2,i)*CRP(1,i)-CRP(3,i)));
    CBN(2,3,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(2*(CRP(2,i)*CRP(3,i)+CRP(1,i)));
    CBN(3,1,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(2*(CRP(3,i)*CRP(1,i)+CRP(2,i)));
    CBN(3,2,i)= 1/(1+CRP(:,i)'*CRP(:,i))*(2*(CRP(3,i)*CRP(2,i)-CRP(1,i)));
    
    CBNConstraint(:,:,i) = CBN(:,:,i)*inv(CBN(:,:,i)); %Finding the constraint of CBN
    CBNDeter(i) = det(CBNConstraint(:,:,i));
end

%figure(5)              %check if it works
%plot(t,CBNDeter')


%% Part a
f=65*10^-3; %Setting focal length parameter

%normalizing star vectors
for i=1:length(R)
    R(:,i)=R(:,i)/norm(R(:,i));
end

for i = 1:length(t)   %Loop through each postion vector for calculation
    for j = 1:4      %Loop through each star
        X(j,i)= -f*(CBN(1,1,i)*R(1,j)+CBN(1,2,i)*R(2,j)+CBN(1,3,i)*R(3,j))/(CBN(3,1,i)*R(1,j)+CBN(3,2,i)*R(2,j)+CBN(3,3,i)*R(3,j));
        Y(j,i)= -f*(CBN(2,1,i)*R(1,j)+CBN(2,2,i)*R(2,j)+CBN(2,3,i)*R(3,j))/(CBN(3,1,i)*R(1,j)+CBN(3,2,i)*R(2,j)+CBN(3,3,i)*R(3,j));
        
        X(j,i)= X(j,i)+normrnd(0,10^-6);  %Add "noise" to measurements
        Y(j,i)= Y(j,i)+normrnd(0,10^-6);  %normrnd is superior to randn
        
        %Calculates the plane projection
        B(1,j,i)= (1/sqrt(X(j,i)^2+Y(j,i)^2+f^2))*-X(j,i);
        B(2,j,i)= (1/sqrt(X(j,i)^2+Y(j,i)^2+f^2))*-Y(j,i);
        B(3,j,i)= (1/sqrt(X(j,i)^2+Y(j,i)^2+f^2))*f;
        
        %Finding error and normalizing it
        Projerror(:,j,i) = norm((B(:,j,i)-CBN(:,:,i)*R(:,j)));
    end  
    
end

%Convert the projection error array into 2D so it can be plotted
figure(1)
Proj = Projerror(:,1,:);
Proj = permute(Proj,[1 3 2]);
subplot(2,2,1);
plot(t,Proj);
title('Star 1 Projection Error')
Proj = Projerror(:,2,:);
Proj = permute(Proj,[1 3 2]);
subplot(2,2,2);
plot(t,Proj);
title('Star 2 Projection Error')
Proj = Projerror(:,3,:);
Proj = permute(Proj,[1 3 2]);
subplot(2,2,3);
plot(t,Proj);
title('Star 3 Projection Error')
Proj = Projerror(:,4,:);
Proj = permute(Proj,[1 3 2]);
subplot(2,2,4);
plot(t,Proj);
title('Star 4 Projection Error')

%% Part b
%CRP are the true rodriques parameters

%OLAE method
for i = 1:length(t) 
    Di(:,:,i)=B(:,:,i)-R(:,:);
    Si(:,:,i)=B(:,:,i)+R(:,:);
    CRPest(:,:,i)=pinv(Si(:,:,i)')*Di(:,:,i)';
    CRPestimated(1,i)=-CRPest(1,2,i);
    CRPestimated(3,i)=CRPest(1,3,i);
    CRPestimated(3,i)=-CRPest(2,3,i);
    
end


figure(2) 
subplot(2,1,1);
plot(t,CRP)   %True rodriques parameters
title('True rodriques parameters')
subplot(2,1,2); 
plot(t,CRPestimated)   %Estimated rodriques parameters
title('Estimated rodriques parameters')

%Plot of Quaterion to check if its right
%  figure(9)
%  subplot(2,2,1)
%  plot(t,Q(1,:))
%  subplot(2,2,2)
%  plot(t,Q(2,:))
%  subplot(2,2,3)
%  plot(t,Q(3,:))
%  subplot(2,2,4)
%  plot(t,Q(4,:))


%% Part c
%CBN calculated prior to part a using method from project 1
%DCM2CRP lines 140 to 144
for i = 1:length(t)
CBNerror(i)=det(CBN(:,:,i)*CBN(:,:,i)');
Roderror(i)=det(CRP(:,i)*CRP(:,i)');
RodEsterror(i)=det(CRPestimated(:,i)*CRPestimated(:,i)');
end

figure(3)
subplot(3,1,1);
plot(t,CBNerror)
title('CBN Error')
subplot(3,1,2);
plot(t,Roderror)
title('True Rodriques Error')
subplot(3,1,3);
plot(t,RodEsterror)
title('Estimated Rodriques Error')

%% Part d

Wk = 1; %We only have 1 sensor each star equal wieght
for i = 1:length(t)
    BTemp=zeros(3);
    for j = 1:length(R)
        BTemp(:,:,j) = Wk.*B(:,j,i)*R(:,j)'; %R is saved as a 3X4
    end
    %Solve for values in K matrix
    Bbar(:,:,i) = BTemp(:,:,1)+BTemp(:,:,2)+BTemp(:,:,3)+BTemp(:,:,4); 
    S(:,:,i) = Bbar(:,:,i)+Bbar(:,:,i)';
    r(i) = trace(Bbar(:,:,i));
    Z(:,i) = [Bbar(2,3,i)-Bbar(3,2,i) Bbar(3,1,i)-Bbar(1,3,i) Bbar(1,2,i)-Bbar(2,1,i)]';
    K(:,:,i) = [ r(i) Z(:,i)'; Z(1,i) S(1,:,i); Z(2,i) S(2,:,i); Z(3,i) S(3,:,i);];
    %Finding Eigenvalues
    [EigenVector,EigenValue] = eig(K(:,:,i));
    %Find vector corresponding to max eigenvalue
    [maxValue, maxIndex] = max(EigenValue);
    [maxValue, maxIndex] = max(maxValue);
    QuatParameter(:,i) = EigenVector(:,maxIndex);
    %This^ is the estimated quaterion
end
%Calculation for true quaterions
%DCM2quat lines 21 to 25


figure(4) 
subplot(2,1,1);
plot(t,Q)  
title('True Quaterions')
subplot(2,1,2); 
plot(t,QuatParameter)  
title('Estimated Quaterions')

%% Part e

for i = 1:length(t)
%Quat2DCM
    
 CBNest(:,:,i) = [Q(1,i)^2 + Q(2,i)^2 - Q(3,i)^2 - Q(4,i)^2, 2*(Q(2,i)*Q(3,i)+Q(1,i)*Q(4,i)), 2*(Q(2,i)*Q(4,i)-Q(1,i)*Q(3,i));
          2*(Q(2,i)*Q(3,i)-Q(1,i)*Q(4,i)), Q(1,i)^2 - Q(2,i)^2 + Q(3,i)^2 - Q(4,i)^2, 2*(Q(3,i)*Q(4,i)+Q(1,i)*Q(2,i));
          2*(Q(2,i)*Q(4,i)+Q(1,i)*Q(3,i)), 2*(Q(3,i)*Q(4,i)-Q(1,i)*Q(2,i)), Q(1,i)^2 - Q(2,i)^2 - Q(3,i)^2 + Q(4,i)^2;];
    
CBNesterror(i)=det(CBNest(:,:,i)*CBNest(:,:,i)');   
    
Q(:,i)=QuatParameter(:,i);
%Used DCM2Quat from part d to get quaterions
%Identity error method 
Truequaterror(i) = det(Q(:,i)\Q(:,i));
%error from quaterion multiplication
Truequaterrorm(i) = Q(1,i)^2+Q(2,i)^2+Q(3,i)^2+Q(4,i)^2;
end

figure(5)
subplot(3,1,1)
plot(t,Truequaterrorm)
title('true quat error multiplication method')
subplot(3,1,2)
plot(t,Truequaterror)
title('true quat error identity method')
subplot(3,1,3)
plot(t,CBNesterror)
title('estimated Cbn error')


%% Functions

function bdot = odefunction(t,b)  %Intergates angular velocity over given time span
    Omega=1.84*(10^(-6));
    w1 = Omega*sin(Omega*t);             %Angular velocity array
    w2 = Omega*cos(Omega*t);
    w3 = Omega+t-t;
    
    bdot = .5*[ 0 -w1 -w2 -w3 ; w1 +0 +w3 -w2 ; w2 -w3 +0 +w1 ; w3 +w2 -w1 +0 ]*b;

end
