function [L2] = LEFM_Modulus(E1, E0, D, v, H, rho_i, rho_sea, hw, g, kcrit, x)

d_min = 10; %Initial Specified Notch
d2 = 10; %Initial crevasse depth for isotropic ice case (notch depth)%
hs2 = x.*d2; %Height of meltwater within crevasse  - isotropic ice case%
zs2 = H-d2; % Height from base of ice shelf to bottom of crevasse  -  VARIABLE MODULUS%
n = 0.01; %Increment size%

pw2 = @(eta) 0.*(eta<(d2-hs2))+0.*(eta>d2)+rho_sea.*g.*(hs2-((H-eta)-zs2)).*((eta>=(d2-hs2)) & (eta<=d2)); % Meltwater pressure - VARIABLE MODULUS
phi_Modulus = @(eta) sqrt(tan((pi.*d2)./(2.*H)))./sqrt(1-(cos((pi.*d2)./(2.*H))./cos((pi.*eta)./(2.*H))).^2);
f_1_Modulus = @(eta) 0.3.*(1-(eta./d2).^1.25);
f_2_Modulus = 0.5.*(1-sin((pi.*d2)./(2.*H))).*(2+sin((pi.*d2)./(2.*H)));
M_D_Modulus = @(eta) (2./(sqrt(2.*H))).*(1+f_1_Modulus(eta).*f_2_Modulus).*phi_Modulus(eta); % LEFM weight function for double edge cracks formulation - VARIABLE MODULUS

sigmaxxE = @(eta) (E1.*exp(H./D)-E0.*exp((H-eta)./D)).*(rho_i.*g.*v.*(H-eta)-E0./E1.*rho_i.*g.*v.*H.*exp(-eta./D))./(E1.*exp(H./D)-E0.*exp((H-eta)./D)-v.*E1.*exp(H./D)+v.*E0.*exp((H-eta)./D)) -( v./(v-1).*E0./E1.*D.*rho_i.*g.*H.*(exp(-H./D)-1)+v./(v-1).*rho_i.*g.*H.^2./2-rho_sea.*g.*hw.^2./2)./(E0.*D.*(exp(H./D)-1)-E1.*H.*exp(H./D)).*(E1.*exp(H./D)-E0.*exp((H-eta)./D)) + pw2(eta); % Far field longitudinal stress - VARIABLE MODULUS (including meltwater pressure)

G_Modulus = @(eta) M_D_Modulus(eta).*sigmaxxE(eta); %Equation for SIF - VARIABLE MODULUS%
f_Modulus = integral(G_Modulus,0,d2); %Evaluates SIF over entire depth - VARIABLE MODULUS%

while f_Modulus > kcrit %LEFM Calculation to find stabilised crevasse depth - VARIABLE MODULUS 
  d2 = d2+n; % Iterated crevasse depth - VARIABLE MODULUS
  hs2 = x.*d2; % Iterated meltwater height - VARIABLE MODULUS
  zs2 = H-d2; % Iterated height from base of ice shelf to bottom of crevasse - VARIABLE MODULUS%
  
pw2 = @(eta) 0.*(eta<(d2-hs2))+0.*(eta>d2)+rho_sea.*g.*(hs2-((H-eta)-zs2)).*((eta>=(d2-hs2)) & (eta<=d2)); % Meltwater pressure for constant density
phi_Modulus = @(eta) sqrt(tan((pi.*d2)./(2.*H)))./sqrt(1-(cos((pi.*d2)./(2.*H))./cos((pi.*eta)./(2.*H))).^2);
f_1_Modulus = @(eta) 0.3.*(1-(eta./d2).^1.25);
f_2_Modulus = 0.5.*(1-sin((pi.*d2)./(2.*H))).*(2+sin((pi.*d2)./(2.*H)));
M_D_Modulus = @(eta) (2./(sqrt(2.*H))).*(1+f_1_Modulus(eta).*f_2_Modulus).*phi_Modulus(eta); % LEFM weight function for double edge cracks formulation - ISOTROPIC ICE

sigmaxxE = @(eta) (E1.*exp(H./D)-E0.*exp((H-eta)./D)).*(rho_i.*g.*v.*(H-eta)-E0./E1.*rho_i.*g.*v.*H.*exp(-eta./D))./(E1.*exp(H./D)-E0.*exp((H-eta)./D)-v.*E1.*exp(H./D)+v.*E0.*exp((H-eta)./D)) -( v./(v-1).*E0./E1.*D.*rho_i.*g.*H.*(exp(-H./D)-1)+v./(v-1).*rho_i.*g.*H.^2./2-rho_sea.*g.*hw.^2./2)./(E0.*D.*(exp(H./D)-1)-E1.*H.*exp(H./D)).*(E1.*exp(H./D)-E0.*exp((H-eta)./D)) + pw2(eta); 

G_Modulus = @(eta) M_D_Modulus(eta).*sigmaxxE(eta); %Equation for SIF - VARIABLE MODULUS%

f_Modulus = integral(G_Modulus,0,d2); %Evaluates SIF over entire depth - VARIABLE MODULUS%

  if d2>H
      break 
  end

end

y2 = max(d2,d_min); %Compares the crevasse depth due to LEFM with the initial crevasse specified in the geometry% 
L2 = round(y2/H,4); % Normalised crevasse depth - VARIABLE DENSITY

end