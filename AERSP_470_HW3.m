%Ryan Peters%  %AERSP 470%  %HW 3%
%Problem 2%
%%
%Given a laminate with a ply stack of [0/+theta/-theta]s and its material properties, plot the
%axial strains vs theta with a range of 0<theta<pi/2 when the laminate is
%applied with a force of Nx=-100000 N/m and Ny=-200000 N/m .
%%
clc
%Defining variables and constants

%NM Matrix, thickness, stiffness matrix, transformation matrices for strains
%and stress
NM=[-100000 -200000 0 0 0 0].';
t=0.0001;
Q=[201.6 8.06 0; 8.06 40.3 0; 0 0 30];
%%
%Set step size for iterative loop and declare strain-curvature matrix
dx=pi/100;
strains(6,51)=zeros;
j=1;
%%
%Create loop to calculate strain-curvature vectors
for theta=0:dx:pi/2
    
    %Define ply stack information for [0/+theta/-theta]s
    ply_stack=[0 3*t 2*t; theta*180/pi 2*t t; -theta*180/pi t 0; 
               -theta*180/pi 0 -t; theta*180/pi -t -2*t; 0 -2*t -3*t];
    
    %Calculate A, D, and B matrices
    A=zeros(3);
    B=zeros(3);
    D=zeros(3);
    for i=1:6
        A=A+t*Q_thetai(Q,ply_stack(i,1));
    end
   
    for i=1:6
        D=D+(((ply_stack(i,2)^3)-(ply_stack(i,3)^3))/3)*Q_thetai(Q,ply_stack(i,1));
    end
    
    for i=1:6
        B=B+(((ply_stack(i,2)^2)-(ply_stack(i,3)^2))/2)*10^3*(Q_thetai(Q,ply_stack(i,1)));
    end
    
    %Define ABD matrix and calculate elongation and curvature
    ABD=[A B; B D];
    strains(:,j)=(ABD)\NM;
    j=j+1;
end
%%
%Plotting strains for x, y , and xy vs. theta
plot(0:dx:pi/2,strains(1,:),'g',...
     0:dx:pi/2,strains(2,:),'k',...
     0:dx:pi/2,strains(3,:),'--b')
title('\theta vs. Strain')
xlabel('\theta (rad)')
ylabel('Microstrain (10^{-6})')
legend('\epsilon_x','\epsilon_y','\gamma_{xy}')
%%
function output=Q_thetai(Q,A)
Tstress=[(cosd(A))^2 (sind(A))^2 2*cosd(A)*sind(A);
    (sind(A))^2 (cosd(A))^2 -2*cosd(A)*sind(A);
    -cosd(A)*sind(A) cosd(A)*sind(A) (cosd(A))^2-(sind(A))^2];

Tstrain=[(cosd(A))^2 (sind(A))^2 cosd(A)*sind(A);
    (sind(A))^2 (cosd(A))^2 -cosd(A)*sind(A);
    -2*cosd(A)*sind(A) 2*cosd(A)*sind(A) (cosd(A))^2-(sind(A))^2];

%Calculate stiffness matrix for a specific ply in laminate axes
output=Tstress^(-1)*Q*Tstrain;
end
%%
%Graph explanation:

%As shown on the graph, when the angle=0, the strain in the y direction is the largest while
%the x-direction is very close to zero.  As the plys are rotated and the angle=pi/2, the
%y-strain approaches -2 microstrain while the x-strain increases to around
%-1.5 microstrain.  The strains even out when angle=pi/2 since the plys are
%aligned with the applied forces and the y-strain is smaller due to the
%larger E(longitudinal) value of the material.

