function AeroPara = Aero(alpha)


% [getcoeff is just getting the required aerodynamic coefficients and
% derivatives for pullup]

alpha = alpha*180/pi; % now alpha is in degree
a_o = aero_lookup_with_alpha_coefficients(alpha);
% all derivatives are in per radians

Cm0        = a_o(32);
Cm_alpha    = a_o(33)*180/pi;
Cm_q       = a_o(34)*180/pi;
Cm_del     = (a_o(35)+a_o(36))*180/pi;

CL0        = a_o(27);
CL_alpha    = a_o(28)*180/pi;
Cl_q        = a_o(29)*180/pi;
CL_del      = (a_o(30)+a_o(31))*180/pi;

CD0        = a_o(22);
CD_alpha    = a_o(23)*180/pi;
CD_q      = a_o(24)*180/pi;
CD_del    = (a_o(25)+a_o(26))*180/pi;

AeroPara = [Cm0,Cm_alpha,Cm_q,Cm_del,CL0,CL_alpha ,Cl_q,CL_del,CD0,CD_alpha,CD_q,CD_del]';
end
