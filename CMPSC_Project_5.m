%Ryan Peters
%CMPSCI 201, Sec 13
%Project 5
%December 7, 2017

%Problem 1: Program to analyze thermal conductivities of different
%formulations from different plants for data in a supplied file using
%arrays, use a nested loop and if statement to check if values are between 3
%and 10 inclusive, output error message if they are not

clear 
clc

Array=importdata('TC.mat'); %import data from file
arraysize=size(Array);

for x= 1:arraysize(1) 
    
    maxtc(x,1)=max(Array(x,:));
    mintc(x,1)=min(Array(x,:));
    avgtc(x,1)=(sum(Array(x,1))- mintc(x,1)-maxtc(x,1))/(arraysize(2)-2);
    vector(x,1)=x;
end

NosuspectTC=0; 

for x=1:arraysize(1)
    for z=1:arraysize(2)
    if (Array(x,z)>10 || Array(x,z)<3)
        NosuspectTC=NosuspectTC+1;
        fprintf('The value of %f from plant %d for formulation %d is suspect\n', Array(x,z),z,x);
    end
    end
end
fprintf('The number of suspect values is %d\n',NosuspectTC); %output total number of suspect values

Outputarray=[vector,mintc,maxtc,avgtc]; %define output array vectors

for x=1:arraysize(1) 
    fprintf('Formulation %d had a minimum thermal conductivity of %.1f a maximum thermal conductivity of %.1f, and an average thermal conductivity of %.2f \n',Outputarray(x,1),Outputarray(x,2), Outputarray(x,3), Outputarray(x,4));
end
%%
% Problem 2: program to plot bacterial growth as a function
% of time, using linspace function to create an array from 0-7 in .1 steps
% Equation y=(yi)*exp(K*T)

k=1.386;
yi = 1;
t = linspace(0, 7, 71);
y=yi*exp(k*t);

figure
subplot(2,2,1)
plot(t,y)
title('Time Vs. Bacteria Growth') %linearly scaled axes
ylabel('Number of bacteria')
xlabel('Time (hours)')

subplot(2,2,2)
semilogy(t,y)
title('Time Vs. log(Bacteria Growth)') %logarithmic scaled y axis
ylabel('log(Number of bacteria)')
xlabel('Time (hours)')

subplot(2,2,3)
semilogx(t,y)
title('log(Time) Vs. Bacteria Growth') %logarithmic scaledx axis
ylabel('Number of bacteria')
xlabel('log(Time) (hours)')

subplot(2,2,4)
loglog(t,y)
title('log(Time) Vs. log(Bacteria Growth)') %logarithmically scaled axes
ylabel('log(Number of bacteria)')
xlabel('log(Time) (hours)')


%Problem 3: program to create estimates for the constant pi using the monte
%carlo simulation and a circle with radius 1 in a 2x2 square both centered
%at 0. The program should use row vectors with 10000 x and y values and
%1000000 x and y values to estimate pi using the formula numpointscircle=N*pi/4
%The number points in the circle is found by x^2 + y^2 <= 1


ranx1=-1+(1+1)*rand(1,10000); 
rany1=-1+(1+1)*rand(1,10000);

Len1x=length(ranx1); %find length
Len1y=length(rany1);

rad1=ranx1.^2+rany1.^2; 
accept1=find(rad1<=1);
m=numel(accept1);

Pi1=m/(.25*Len1x); %close but less exact than larger array estimate

ranx2=-1+(1+1)*rand(1,1000000);
rany2=-1+(1+1)*rand(1,1000000);

Len2x=length(ranx2); %find vector length
Len2y=length(rany2);

rad2=ranx2.^2+rany2.^2; 
accept2=find(rad2<=1);
m=numel(accept2);

Pi2=m/(.25*Len2x); %more exact because of larger sample size

text1= 'PI from 10000 random values=';
disp(text1);
disp(Pi1);
text2= 'PI from 1000000 random values=';
disp(text2);
disp(Pi2);