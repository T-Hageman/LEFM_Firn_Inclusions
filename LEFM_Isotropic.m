function [L1] = LEFM_Isotropic(v, H, rho_i, rho_sea, hw, g, kcrit, x)

d_min = 10; %Initial Specified Notch
d1 = 10; %Initial crevasse depth for isotropic ice case (notch depth)%
hs1 = x.*d1; %Height of meltwater within crevasse  - isotropic ice case%
zs1 = H-d1; % Height from base of ice shelf to bottom of crevasse  - ISOTROPIC ICE%
n = 0.01; %Increment size%

pw1 = @(eta) 0.*(eta<(d1-hs1))+0.*(eta>d1)+rho_sea.*g.*(hs1-((H-eta)-zs1)).*((eta>=(d1-hs1)) & (eta<=d1)); % Meltwater pressure for constant density
phiconst = @(eta) sqrt(tan((pi.*d1)./(2.*H)))./sqrt(1-(cos((pi.*d1)./(2.*H))./cos((pi.*eta)./(2.*H))).^2);
f_1const = @(eta) 0.3.*(1-(eta./d1).^1.25);
f_2const = 0.5.*(1-sin((pi.*d1)./(2.*H))).*(2+sin((pi.*d1)./(2.*H)));
M_Dconst = @(eta) (2./(sqrt(2.*H))).*(1+f_1const(eta).*f_2const).*phiconst(eta); % LEFM weight function for double edge cracks formulation - ISOTROPIC ICE

sigmaxx = @(eta) v./(1-v).*(rho_i.*g.*((H-eta)-H./2))-0.5.*rho_sea.*g.*hw.^2./H +pw1(eta); % Far field longitudinal stress - ISOTROPIC ICE (including meltwater pressure)

Gconst =  @(eta) M_Dconst(eta).*sigmaxx(eta); %Equation for SIF - ISOTROPIC ICE%
fconst = integral(Gconst,0,d1); %Evaluates SIF over entire depth - ISOTROPIC ICE%

while fconst > kcrit %LEFM Calculation to find stabilised crevasse depth - ISOTROPIC ICE
  d1 = d1+n; % Iterated crevasse depth - ISOTROPIC ICE
  hs1 = x.*d1; % Iterated meltwater height - ISOTROPIC ICE
  zs1 = H-d1; % Iterated height from base of ice shelf to bottom of crevasse - ISOTROPIC ICE%

  pw1 = @(eta) 0.*(eta<(d1-hs1))+0.*(eta>d1)+rho_sea.*g.*(hs1-((H-eta)-zs1)).*((eta>=(d1-hs1)) & (eta<=d1)); % Iterated meltwater pressure - ISOTROPIC ICE
  phiconst = @(eta) sqrt(tan((pi.*d1)./(2.*H)))./sqrt(1-(cos((pi.*d1)./(2.*H))./cos((pi.*eta)./(2.*H))).^2);
  f_1const = @(eta) 0.3.*(1-(eta./d1).^1.25);
  f_2const = 0.5.*(1-sin((pi.*d1)./(2.*H))).*(2+sin((pi.*d1)./(2.*H)));
  M_Dconst = @(eta) (2./(sqrt(2.*H))).*(1+f_1const(eta).*f_2const).*phiconst(eta); % LEFM weight function for double edge cracks formulation - ISOTROPIC ICE


  sigmaxx = @(eta) v./(1-v).*(rho_i.*g.*((H-eta)-H./2))-0.5.*rho_sea.*g.*hw.^2./H + pw1(eta); % Far field longitudinal stress for isotropic ice (including meltwater pressure)
  Gconst =  @(eta) M_Dconst(eta).*sigmaxx(eta); %Equation for SIF - ISOTROPIC ICE%
  fconst = integral(Gconst,0,d1); %Evaluates SIF over entire depth - ISOTROPIC ICE%

  if d1>H
      break 
  end
  
end

y1 = max(d1,d_min); %Compares the crevasse depth due to LEFM with the initial crevasse specified in the geometry% 
L1 = round(y1/H,4); % Normalised crevasse depth - VARIABLE DENSITY


end