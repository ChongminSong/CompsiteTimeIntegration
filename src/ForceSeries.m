function [ft] = ForceSeries(order,fHist,ns,dt)
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

pf = order + 1;
np = pf + 0;
s = forceSamplingPoints(np);
t = dt*reshape(repmat((0:ns-1), np,1)+repmat(s,1,ns),[],1);
f = reshape(fHist(t),np,[]);
T = transMtxPointsToPoly(s, pf);
ft=f'*T;

end