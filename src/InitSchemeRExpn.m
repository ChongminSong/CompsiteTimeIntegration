function [r, prcoe, cf] = InitSchemeRExpn(scheme, p, rho )
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

if contains(scheme,'MP1')
    r =  MP1SchemeRoot(p);
else
    r = MschemeRoot(p, rho);
end

pcoe = pCoefficients(p, r);

disp(['r = ', num2str(r)]);
disp(['Seclected rhoInfty = ', num2str(abs(pcoe(end)), '%.5f')]);

qrcoe = [zeros(1,p) 1];
qcoe = shiftPolycoe(qrcoe,r);

cf = TimeIntgCoeffForce(pcoe,qcoe);

prcoe = shiftPolycoe(pcoe,r);

end