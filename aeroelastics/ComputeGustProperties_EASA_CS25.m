%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [Uds_grid, U_sig] = ComputeGustProperties_EASA_CS25(Fg, Altitude, H_grid_meters, VTAS, DynPress_Dive, Mach_Dive)



R = 287.053;     %[J/kg/deg]
Kappa = 1.4;
T_0 = 288.15;    %[K]
a_0 = 340.294;   %[m/s]
rho_0 = 1.225;   %[kg/m^3]

T = T_0 - 0.0065*Altitude;
a = sqrt(Kappa*R*T);

Mach = VTAS/a;

if Altitude < 11000 % Atmosphere
    sigma=(1-6.5e-3/288.15*(Altitude))^(9.80665/287.053/6.5e-3-1);
else
    sigma=exp(-9.80665/287.053*(Altitude-11000)/(288.15*0.75187))*0.29708;
end
        
rho = rho_0*sigma;

DynPress = rho*VTAS^2/2;

%% see CS 25.341 (a) "Discrete Gust Design Criteria"        
       
if Altitude < 15000*0.3048
    U_ref = (56-(56-44)/4572 * Altitude)*0.3048;                              % U_ref in m/s EAS
elseif and(Altitude>4572,Altitude<=18288) 
    % [1996] version of CS 25.341
    %U_ref = (44-(44-26)/(15240-4572) * (Altitude-4572))*0.3048;              % U_ref in m/s EAS
    % [2010] version of CS 25.341
    U_ref = (44-(44-20.86)/(18288-4572) * (Altitude-4572))*0.3048;            % U_ref in m/s EAS  
elseif Altitude>18288
    Uref = 6.36;
end
        
        
H_grid = H_grid_meters/0.3048;
        
% Conversion from EAS 2 TAS:
        
%At aeroplane speeds between V_B and V_C:
Uds_grid = U_ref*Fg*(H_grid/350).^(1/6)/sqrt(sigma);                          %TAS in m/s compute from FAR
        
% Compute reduced values for dive speed V_D:
if (Mach >= Mach_Dive) || (DynPress >= DynPress_Dive)  
    Uds_grid = Uds_grid*0.5;                                                  %TAS in m/s compute from FAR
end
                
%% see CS 25.341 (b) "Continuous Turbulence Design Criteria"  [2010]      
if Altitude < 7315
    U_sig_ref = (24.08-27.43)/7315*Altitude + 27.43;                          %m/s TAS
else
    U_sig_ref = 24.08;                                                        %m/s TAS
end
%At aeroplane speeds between V_B and V_C:   
U_sig = U_sig_ref * Fg;                                                       %m/s TAS
        
% Compute reduced values for dive speed V_D:
if (Mach >= Mach_Dive) || (DynPress >= DynPress_Dive)  
    U_sig = U_sig*0.5;                                                  %TAS in m/s compute from FAR
end
            
% Note: between V_C and V_D linear interpolation is allowed for U_sig.
% Not so for U_ds (according to CS-25)


end

