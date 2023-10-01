function [T] = transMtxPointsToPoly(s, nf)
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

s  = reshape(s,[],1);
if length(s) >= nf
    a = (s-0.5).^(0:nf-1); %
    T = a/(a'*a);
else
    disp([' ******* The number of sampling points: ', num2str(length(s))]);
    disp(['         should be larger than the number of terms of polynomial: ' ...
        num2str(nf)]);
end