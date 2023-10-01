function [r] = MschemeRoot(M, rhoInfty)
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

RHS = [1,  1, -1,  1, -1, -1]*rhoInfty;
ir  = [1,  2,  2,  2,  3,  3];
j = 0:M;
pMcoe = ((-1).^j).*factorial(M)./factorial(j)./(factorial(M-j).^2);
pMcoe(end) = pMcoe(end) - RHS(M);
rs = sort( roots(pMcoe) );
r  = rs(ir(M));
end
