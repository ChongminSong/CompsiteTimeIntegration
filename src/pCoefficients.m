function pcoe = pCoefficients(M, r)
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

pcoe = zeros(1,M+1);
for ii = 0:M
    j = 0:ii;
    p = ((-1).^j).*factorial(M)./factorial(M-j)./factorial(j)./factorial(ii-j);
    pcoe(ii+1) = p*(r.^(M-j))';
end
end
