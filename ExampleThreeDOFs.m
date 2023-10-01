%%
% Example in Section 5.2 of:
% @article{SONG2024116473,
% author = {Chongmin Song and Xiaoran Zhang},
% title = {High-order composite implicit time integration schemes based on rational approximations for elastodynamics},
% journal = {Computer Methods in Applied Mechanics and Engineering},
% volume = {418},
% pages = {116473},
% year = {2024},
% issn = {0045-7825},
% doi = {https://doi.org/10.1016/j.cma.2023.116473},
% url = {https://www.sciencedirect.com/science/article/pii/S0045782523005972},
% }
% Developed by Chongmin Song (c.song@unsw.edu.au)
%              Xiaoran Zhang (xiaoran.zhang1@student.unsw.edu.au)  
% Distributed under the MIT License (https://opensource.org/license/mit/)

%%
clear; close all;
dbstop if error;

addpath(".\src\");

%% Parameters for composite time integration
scheme = 'MP'; % 'M' or 'MP1'
nSubStep = 3; % 1 to 6 for M-scheme; 1, 2, 3, 5 for (M+1)-scheme
rhoInfty = 0; % for M-scheme only

%% Parameters of 3-DOF system
tmax = 5010;       % Simulation time
dt = 0.14;      % Time step size
% Spring constant
k1 = 10^7;
k2 = 1;
% Mass
m1 = 0;
m2 = 1;
m3 = 1;

np = 2;
M = [m2 0; 0, m3]; % Mass matrix
C = zeros(2); % Damping matrix
K = [k1+k2, -k2; -k2, k2]; % Stiffness matrix
F0 = [1; 0]*k1; % Force vector
u0 = [0; 0]; % Initial displacement
v0 = u0;  % Initial velocity
BC_Accl = [];

% Save results for particular degrees of freedom
pDOF = [1;2];

fHist = @(t) sin(1.2*t);      % excitation signal

ns = floor(tmax/dt)+1;
tp = (0:0.02:tmax);

[uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp);
R1Ref = k1*(fHist(tp)-uRef(1,:));

[dsp,vel,acc] = TimeSolverRExpn(scheme,nSubStep, rhoInfty, ...
                             ns,dt,fHist, K,M,C,F0,u0,v0,pDOF);
tn = (0:ns-1)*dt;
R1 = k1*((fHist(tn))'-dsp(:,1));

figure(1)
plot(tp,uRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,dsp(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Displacement')
title('Displacement of Mass 2');

figure(2)
plot(tp,uRef(2,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,dsp(:,2), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Displacement')
title('Displacement of Mass 3');

figure(3)
plot(tp,vRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,vel(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Velocity')
title('Velocity of Mass 2');

figure(4)
plot(tp,vRef(2,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,vel(:,2), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Velocity')
title('Velocity of Mass 3');

figure(5)
plot(tp,aRef(1,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,acc(:,1), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Acceleration')
title('Acceleration of Mass 2');

figure(6)
plot(tp,aRef(2,:), '-r',"DisplayName",'Reference')
hold on
plot(tn,acc(:,2), '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Acceleration')
title('Acceleration of Mass 3');

figure(7)
plot(tp,R1Ref, '-r',"DisplayName",'Reference')
hold on
plot(tn,R1, '--b',"DisplayName",'Present')
grid on
legend('show')
xlabel('Time');
ylabel('Force')
title('Rection Force');


function [uRef, vRef, aRef] = ThreeDOFsRefSln(M,K,F0,tp)
[Vec,D] = eig(full(K));
Mg = Vec'*M*Vec;
Kg = Vec'*K*Vec;
Fg = Vec'*F0;
o1 = sqrt(D(1,1));
o2 = sqrt(D(2,2));
Uex = @(o,t) 1/(o*o-1.2^2)*(sin(1.2*t) - 1.2/o*sin(o*t));
Uref = @(o,t) 1/(o*o-1.2^2)*(sin(1.2*t)       );
Vex = @(o,t) 1.2/(o*o-1.2^2)*(cos(1.2*t) - cos(o*t));
Vref = @(o,t) 1.2/(o*o-1.2^2)*(cos(1.2*t)     );
Aex = @(o,t) 1.2/(o*o-1.2^2)*(-1.2*sin(1.2*t) + o*sin(o*t));
Aref = @(o,t) 1.2/(o*o-1.2^2)*(-1.2*sin(1.2*t)     );
U  = [Fg(1)*Uex(o1,tp); Fg(2)*Uref(o2,tp)];
uRef = Vec*U;
V  = [Fg(1)*Vex(o1,tp); Fg(2)*Vref(o2,tp)];
vRef = Vec*V;
A  = [Fg(1)*Aex(o1,tp); Fg(2)*Aref(o2,tp)];
aRef = Vec*A;
end

