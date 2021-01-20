%% Function I want to integrate
function f = integranda(a, z_M, z_MP, delta_z, x, kappa)
    
    f = exp(-1i*kappa*sqrt(a^2+(z_MP-x*delta_z-z_M)^2))/sqrt((a^2+(z_MP-x*delta_z-z_M)^2));
    %f = 1/(sqrt(a^2+(z_MP-x*delta_z-z_M)^2))-1i*kappa-kappa^2*sqrt(a^2+(z_MP-x*delta_z-z_M)^2)/2+1/6*1i*kappa^3*sqrt(a^2+(z_MP-x*delta_z-z_M)^2)^2+kappa^4*sqrt(a^2+(z_MP-x*delta_z-z_M)^2)^3/24
end
