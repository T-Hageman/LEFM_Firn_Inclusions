function [L] = LEFM_Density(v, D, H, rho_i, rho_sea, hw, rho_s, g, kcrit, x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

d_min = 10; %Initial Specified Notch
d = 10; %Initial crevasse depth for variable density case (notch depth)%
hs = x.*d; %Height of meltwater within crevasse - variable density case%
zs = H-d; % Height from base of ice shelf to bottom of crevasse - VARIABLE DENSITY%
n = 0.01; %Increment size%

pw = @(eta) 0.*(eta<(d-hs))+0.*(eta>d)+rho_sea.*g.*(hs-((H-eta)-zs)).*((eta>=(d-hs)) & (eta<=d)); % Meltwater pressure - VARIABLE DENSITY
phi = @(eta) sqrt(tan((pi.*d)./(2.*H)))./sqrt(1-(cos((pi.*d)./(2.*H))./cos((pi.*eta)./(2.*H))).^2);
f_1 = @(eta) 0.3.*(1-(eta./d).^1.25);
f_2 = 0.5.*(1-sin((pi.*d)./(2.*H))).*(2+sin((pi.*d)./(2.*H)));
M_D = @(eta) (2./(sqrt(2.*H))).*(1+f_1(eta).*f_2).*phi(eta); % LEFM weight function for double edge cracks formulation - VARIABLE DENSITY
sigmaxxrho = @(eta) v./(1-v).*(-rho_i.*g.*(eta)+(rho_i-rho_s).*D.*g.*(1-exp(-eta./D)))+v./(1-v).*rho_i.*g.*H./2 - rho_sea.*g.*hw.^2./2./H-v./(1-v).*(rho_i-rho_s).*D.*g+v./(1-v).*(rho_i-rho_s).*(D.^2).*g./H.*(1-exp(-H./D)) + pw(eta); % Far field longitudinal stress - VARIABLE DENSITY (including meltwater pressure)

G = @(eta) M_D(eta).*sigmaxxrho(eta); %Equation for SIF - VARIABLE DENSITY%

f = integral(G,0,d); %Evaluates SIF over entire depth - VARIABLE DENSITY%

while f > kcrit %LEFM Calculation to find stabilised crevasse depth - VARIABLE DENSITY 
  d = d+n; % Iterated crevasse depth - VARIABLE DENSITY
  hs = x.*d; % Iterated meltwater height - VARIABLE DENSITY
  zs = H-d; % Iterated height from base of ice shelf to bottom of crevasse - VARIABLE DENSITY%
  
  pw = @(eta) 0.*(eta<(d-hs))+0.*(eta>d)+rho_sea.*g.*(hs-((H-eta)-zs)).*((eta>=(d-hs)) & (eta<=d)); % Iterated meltwater pressure - VARIABLE DENSITY
  phi = @(eta) sqrt(tan((pi.*d)./(2.*H)))./sqrt(1-(cos((pi.*d)./(2.*H))./cos((pi.*eta)./(2.*H))).^2);
  f_1 = @(eta) 0.3.*(1-(eta./d).^1.25);
  f_2 = 0.5.*(1-sin((pi.*d)./(2.*H))).*(2+sin((pi.*d)./(2.*H)));
  M_D = @(eta) (2./(sqrt(2.*H))).*(1+f_1(eta).*f_2).*phi(eta); % LEFM weight function for double edge cracks formulation - VARIABLE DENSITY

  sigmaxxrho = @(eta) v./(1-v).*(-rho_i.*g.*(eta)+(rho_i-rho_s).*D.*g.*(1-exp(-eta./D)))+v./(1-v).*rho_i.*g.*H./2 - rho_sea.*g.*hw.^2./2./H-v./(1-v).*(rho_i-rho_s).*D.*g+v./(1-v).*(rho_i-rho_s).*(D.^2).*g./H.*(1-exp(-H./D)) + pw(eta); 

  G = @(eta) M_D(eta).*sigmaxxrho(eta); %Equation for SIF - VARIABLE DENSITY%

  f = integral(G,0,d); %Evaluates SIF over entire depth - VARIABLE DENSITY%
  
  if d>H
      break 
  end

end
y = max(d,d_min); %Compares the crevasse depth due to LEFM with the initial crevasse specified in the geometry% 
L = round(y/H,4); % Normalised crevasse depth - VARIABLE DENSITY

end