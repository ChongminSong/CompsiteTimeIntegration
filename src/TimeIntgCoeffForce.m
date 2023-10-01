function [C] = TimeIntgCoeffForce(p,q)
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

M = length(q) - 1;
tmp = p - q;
C = zeros(M+1,M);
C(1,:) = tmp(2:end);  %ltx order reduced by one {A}^{-1}\left({P}-{Q}\right) 
for k = 1:M
    tmp = ((-1/2)^k)*(p-((-1)^k)*q); %ltx order of $C_{k-1}$ is lower than $q$ by one
    tmp(1:M) = tmp(1:M) + k*C(k,:); 
    C(k+1,:) =  tmp(2:end);
end