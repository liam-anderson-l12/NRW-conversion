L = 0.01;% length of sample
mu_r = 1;% relative permeability of sample
lambdac = 6.557E9;% WR-90 cutoff frequency Hz
[filename,folder] = uigetfile('*.s2p','Select the S-parameter file');
% token = strtok(filename, '.');
% result_folder = strcat(folder, token, '\Result');
% startfolder = cd(result_folder);
% datafilename = strcat(token, '.mat');

[freq, S, freq_noise, data_noise, Zo] = SXPParse(filename);
npoints = size(S,3);
eps_o = 8.854187E-12;
c = 2.99792458E8;
X = zeros(npoints,1);
Gamma = zeros(npoints,1);
eps_r = zeros(npoints,1);

for f = 1:npoints
    
    lambdao = c/freq(f);
    
    X(f) = ((S(1,1,f)^2) - (S(2,1,f)^2) + 1)/(2*S(1,1,f));
    Gamma1 = X(f) + sqrt((X(f)^2) - 1);
    Gamma2 = X(f) - sqrt((X(f)^2) - 1);
    if abs(Gamma1) <= 1
        Gamma = Gamma1;
    elseif abs(Gamma2) <=1
        Gamma = Gamma2;
    else
        Gamma = NaN;
    end
    
    T(f) = (S(1,1,f) + S(2,1,f) - Gamma)/(1-(S(1,1,f) + S(2,1,f))*Gamma);
    
    eps_r(f) = ((lambdao^2)/mu_r)*((1/lambdac^2)-((1/(2*pi()*L))*log(1/T(f)))^2);
    
    
%     eps_r(f) = mu_r*(((1-Gamma)^2)/(1+Gamma)^2)*(1-(lambdao^2)/(lambdac^2))+((lambdao^2)/(lambdac^2))*(1/mu_r);
    perm(f) = real(eps_r(f));
    cond(f) = imag(eps_r(f))*2*pi()*freq(f)*eps_o;
    
end



