function [L3] = LEFM_Density_Modulus(E1, E0, D, v, H, rho_i, rho_sea, rho_s, hw, g, kcrit, x)

d_min = 10; %Initial Specified Notch
d3 = 10; %Initial crevasse depth for isotropic ice case (notch depth)%
hs3 = x.*d3; %Height of meltwater within crevasse  - isotropic ice case%
zs3 = H-d3; % Height from base of ice shelf to bottom of crevasse  -  VARIABLE MODULUS & DENSITY%
n = 0.01; %Increment size%

pw3 = @(eta) 0.*(eta<(d3-hs3))+0.*(eta>d3)+rho_sea.*g.*(hs3-((H-eta)-zs3)).*((eta>=(d3-hs3)) & (eta<=d3)); % Meltwater pressure -  VARIABLE MODULUS & DENSITY
phi_Modulus_Density = @(eta) sqrt(tan((pi.*d3)./(2.*H)))./sqrt(1-(cos((pi.*d3)./(2.*H))./cos((pi.*eta)./(2.*H))).^2);
f_1_Modulus_Density = @(eta) 0.3.*(1-(eta./d3).^1.25);
f_2_Modulus_Density = 0.5.*(1-sin((pi.*d3)./(2.*H))).*(2+sin((pi.*d3)./(2.*H)));
M_D_Modulus_Density = @(eta) (2./(sqrt(2.*H))).*(1+f_1_Modulus_Density(eta).*f_2_Modulus_Density).*phi_Modulus_Density(eta); % LEFM weight function for double edge cracks formulation -  VARIABLE MODULUS & DENSITY

sigmaxxErho = @(eta) (-(E1.*exp(H./D)-E0.*exp((H-eta)./D)).*(((v./(v-1)).*0.5.*rho_i.*g.*H.^2+0.5.*rho_sea.*g.*hw.^2+(v./(v-1)).*g.*H.*(rho_s-rho_i).*D+(v./(v-1)).*g.*(rho_i-rho_s).*(1-exp(-H./D)).*D.^2)./(E1.*H.*exp(H./D)-E0.*D.*(exp(H./D)-1)))+ (E1.*exp(H./D)-E0.*exp((H-eta)./D)).*(v./(v-1)).*g.*exp(-H./D).*(rho_s.*exp(H./D)+exp((H-eta)./D).*(rho_i-rho_s)+rho_i.*exp(H./D).*(H./D-1)-rho_i./D.*(H-eta).*exp(H./D))./((E1.*exp(H./D)-E0.*exp((H-eta)./D))./D)) + pw3(eta); % Far field longitudinal stress - VARIABLE MODULUS & DENSITY (including meltwater pressure)
G_Modulus_Density = @(eta) M_D_Modulus_Density(eta).*sigmaxxErho(eta); %Equation for SIF - VARIABLE MODULUS%
f_Modulus_Density = integral(G_Modulus_Density,0,d3); %Evaluates SIF over entire depth - VARIABLE MODULUS%

while f_Modulus_Density > kcrit %LEFM Calculation to find stabilised crevasse depth - VARIABLE MODULUS & DENSITY
  d3 = d3+n; % Iterated crevasse depth - VARIABLE MODULUS & DENSITY
  hs3 = x.*d3; % Iterated meltwater height - VARIABLE MODULUS & DENSITY
  zs3 = H-d3; % Iterated height from base of ice shelf to bottom of crevasse - VARIABLE MODULUS & DENSITY%
  
pw3 = @(eta) 0.*(eta<(d3-hs3))+0.*(eta>d3)+rho_sea.*g.*(hs3-((H-eta)-zs3)).*((eta>=(d3-hs3)) & (eta<=d3)); % Meltwater pressure - VARIABLE MODULUS & DENSITY
phi_Modulus_Density = @(eta) sqrt(tan((pi.*d3)./(2.*H)))./sqrt(1-(cos((pi.*d3)./(2.*H))./cos((pi.*eta)./(2.*H))).^2);
f_1_Modulus_Density = @(eta) 0.3.*(1-(eta./d3).^1.25);
f_2_Modulus_Density = 0.5.*(1-sin((pi.*d3)./(2.*H))).*(2+sin((pi.*d3)./(2.*H)));
M_D_Modulus_Density = @(eta) (2./(sqrt(2.*H))).*(1+f_1_Modulus_Density(eta).*f_2_Modulus_Density).*phi_Modulus_Density(eta); % LEFM weight function for double edge cracks formulation - VARIABLE MODULUS & DENSITY

sigmaxxErho = @(eta) (-(E1.*exp(H./D)-E0.*exp((H-eta)./D)).*(((v./(v-1)).*0.5.*rho_i.*g.*H.^2+0.5.*rho_sea.*g.*hw.^2+(v./(v-1)).*g.*H.*(rho_s-rho_i).*D+(v./(v-1)).*g.*(rho_i-rho_s).*(1-exp(-H./D)).*D.^2)./(E1.*H.*exp(H./D)-E0.*D.*(exp(H./D)-1)))+ (E1.*exp(H./D)-E0.*exp((H-eta)./D)).*(v./(v-1)).*g.*exp(-H./D).*(rho_s.*exp(H./D)+exp((H-eta)./D).*(rho_i-rho_s)+rho_i.*exp(H./D).*(H./D-1)-rho_i./D.*(H-eta).*exp(H./D))./((E1.*exp(H./D)-E0.*exp((H-eta)./D))./D)) + pw3(eta);

G_Modulus_Density = @(eta) M_D_Modulus_Density(eta).*sigmaxxErho(eta); %Equation for SIF - VARIABLE MODULUS & DENSITY%

f_Modulus_Density = integral(G_Modulus_Density,0,d3); %Evaluates SIF over entire depth - VARIABLE MODULUS & DENSITY%

  if d3>H
      break 
  end

end

y3 = max(d3,d_min); %Compares the crevasse depth due to LEFM with the initial crevasse specified in the geometry% 
L3 = round(y3/H,4); % Normalised crevasse depth - ISOTROPIC ICE

end