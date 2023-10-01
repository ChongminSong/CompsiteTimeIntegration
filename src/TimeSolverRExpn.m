function [dsp,vel,acc] = TimeSolverRExpn(scheme,p,rho,ns,dt,fHist,K,M,C,F,u0,v0,pDOF)
%
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

ndof = size(K,1);

K =  dt*dt*sparse(K);
C = dt*sparse(C);
M = sparse(M);
F =  dt*dt*F;
% Initial velocity --> normalized with dt
v0 = dt*v0;

ft = ForceSeries(p,fHist,ns-1,dt);


% initializing variables storing response history for output
dsp = zeros(ns, length(pDOF)); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

dM = decomposition(M);
Fb = (dM\F);

[r, prcoe, cf] = InitSchemeRExpn(scheme, p, rho );

Kd = sparse((r*r)*M + r*C + K);
dKd = decomposition(Kd);

f = zeros(2*ndof,p);
f(:,1)  = [Fb; zeros(ndof,1)];
for ii = 2:p
    f(:,ii) = SolverPadeAx(dM,K,C,f(:,ii-1));
end
Pf = f*cf';

z = zeros(2*ndof,p+1);

% Initial conditions
z(:,1) = [v0; u0];
for ii = 1:p
    z(:,ii+1) = r*z(:,ii) - SolverPadeAx(dM,K,C,z(:,ii));
end

% store responses for output
it = 1;
dsp(it,:) = u0(pDOF);
vel(it,:) = v0(pDOF);
fn = fHist(0); % at t=0
acc(it,:) = r*z(pDOF,1)-z(pDOF,2) + fn*Fb(pDOF);

for it = 2:ns
    is = it - 1;

    z(:,p+1) = (z(:,1:p+1))*prcoe' + Pf*ft(is,:)';
    for ip = p:-1:1
        z(1:ndof,ip) = dKd\(r*(M*z(1:ndof,ip+1)) - K*z(ndof+1:2*ndof,ip+1));
        z(ndof+1:2*ndof,ip) = (z(1:ndof,ip)+z(ndof+1:2*ndof,ip+1))/r;
    end

    % store responses for output
    vel(it,:) = z(pDOF,1);
    dsp(it,:) = z(ndof+pDOF,1);
    fn = fHist(is*dt);
    acc(it,:) = r*z(pDOF,1)-z(pDOF,2) + fn*Fb(pDOF);

end

vel = vel/dt;
acc = acc/(dt*dt);


end
